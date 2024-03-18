from typing import List
from .consts import DB_URL, AUTH, gds, graph, USER, PASSWORD
from itertools import combinations, product
from .func import (
    rSubset,
    ProjectGraph,
    ProjectGraph_WithProperties,
    DirectProduct_Ordered,
    DirectProduct_Unordered,
    ProjectGraphWithPatients,
)
from graphdatascience import GraphDataScience
from tqdm.auto import tqdm

from p_tqdm import p_umap, p_map

# from .p_map import p_umap, p_map
import pandas as pd
import numpy as np
import math
import ast
import itertools
import multiprocessing as mp
import csv

# need to put here the link to synthea output files
link_to_output = "input.keepme/synthea/"
cohorts_ref = [
    "breast_cancer",
    "cerebral_palsy",
    "colorectal_cancer",
    "dialysis",
    "hypertension",
    "normal_sample",
    "prostate_cancer",
]


def Fetch_Output_Synthea(cohort_selection: str = ""):
    # Fetch the input from given link. If no cohorts, put nothing as variable
    link = link_to_output + cohort_selection + "/"
    input_procedures = pd.read_csv(link + "procedures.csv")
    input_conditions = pd.read_csv(link + "conditions.csv")
    input_observations = pd.read_csv(link + "observations.csv")
    input_patients_infos = pd.read_csv(link + "patients.csv")
    return input_procedures, input_conditions, input_observations, input_patients_infos


def Retrieve_Additional_Infos_Synthea(cohort_selection: str = ""):
    # We create a list of additional snomed code we want to have in the snomed encoding of the synthea patients
    # we retrieve the gender, the smoking status and the BMI range of the synthetic patient
    _, _, input_observations, input_patients_infos = Fetch_Output_Synthea(
        cohort_selection
    )
    ## We extract the informations ##
    patients_data_infos = input_patients_infos.loc[:, ("Id", "GENDER")]
    patients_observations = input_observations.loc[
        :, ("PATIENT", "DESCRIPTION", "VALUE")
    ]
    patients_observations.loc[:, "VALUE"] = patients_observations.loc[:, "VALUE"].apply(
        str
    )
    patients_unique_pid = patients_data_infos["Id"].drop_duplicates().tolist()
    patients_observation_unique_index = patients_observations.loc[
        :, ("PATIENT", "DESCRIPTION")
    ].values.tolist()
    patients_observation_unique_index = [
        pid + "_" + description
        for [pid, description] in patients_observation_unique_index
    ]
    patients_observations["Index"] = patients_observation_unique_index
    patients_data_infos = patients_data_infos.set_index("Id")
    patients_observations = patients_observations.set_index("Index")
    patients_tobacco_status = []
    ## We retrieve the wanted infos with a for loop ##
    for pid in patients_unique_pid:
        if patients_data_infos.loc[pid, "GENDER"] == "M":
            patients_data_infos.loc[pid, "GENDER"] = 248153007
        else:
            patients_data_infos.loc[pid, "GENDER"] = 248152002

        try:
            neversmoked = any(
                patients_observations.loc[str(pid) + "_Tobacco smoking status", "VALUE"]
                == "Never smoked tobacco (finding)"
            ) and not all(
                patients_observations.loc[str(pid) + "_Tobacco smoking status", "VALUE"]
                == "Ex-smoker (finding)"
            )
            didsmoke = any(
                patients_observations.loc[str(pid) + "_Tobacco smoking status", "VALUE"]
                == "Ex-smoker (finding)"
            )
        except:
            neversmoked = False
            didsmoke = False
        else:
            neversmoked = any(
                patients_observations.loc[str(pid) + "_Tobacco smoking status", "VALUE"]
                == "Never smoked tobacco (finding)"
            ) and not all(
                patients_observations.loc[str(pid) + "_Tobacco smoking status", "VALUE"]
                == "Ex-smoker (finding)"
            )
            didsmoke = any(
                patients_observations.loc[str(pid) + "_Tobacco smoking status", "VALUE"]
                == "Ex-smoker (finding)"
            )
        if neversmoked:
            patients_tobacco_status.append(266919005)
        elif didsmoke:
            patients_tobacco_status.append(8517006)
        else:
            patients_tobacco_status.append(119955)
        # Cannot put BMI range since the file is time-dependent. At some point the patient
        # had a low BMI index, and then a high BMI index.. cannot summarize this without including
        # a time component.
        # if any(
        #     patients_observations.loc[
        #         str(pid) + "_Body mass index (BMI) [Ratio]", "VALUE"
        #     ]
        #     < 20
        # ):
        #     patients_BMI_range.append(310252000)
        # elif any(
        #     patients_observations.loc[
        #         str(pid) + "_Body mass index (BMI) [Ratio]", "VALUE"
        #     ]
        #     in range(20, 25, 0.01)
        # ):
        #     patients_BMI_range.append(412768003)
        # elif any(
        #     patients_observations.loc[
        #         str(pid) + "_Body mass index (BMI) [Ratio]", "VALUE"
        #     ]
        #     in range(25, 30, 0.01)
        # ):
        #     patients_BMI_range.append(162863004)

    ## Extraction of additional data is done, but now we want it in a nice form for db injection ##
    patients_data_infos["TOBACCO"] = patients_tobacco_status
    patients_data_infos = patients_data_infos.reset_index()
    patient_additional_infos_gender = patients_data_infos.loc[:, ("Id", "GENDER")]
    patient_additional_infos_gender.rename(
        columns={"Id": "pid", "GENDER": "sctid"},
        inplace=True,
    )
    patient_additional_infos_tobacco = patients_data_infos.loc[:, ["Id", "TOBACCO"]]
    patient_additional_infos_tobacco.rename(
        columns={"Id": "pid", "TOBACCO": "sctid"},
        inplace=True,
    )
    patient_additional_infos = pd.concat(
        [patient_additional_infos_gender, patient_additional_infos_tobacco],
        ignore_index=True,
    )
    return patient_additional_infos


def Table_of_Patients_From_Synthea(
    cohort_selection: str = "", pid_rename: bool = False
):
    # Extract the data from synthea csv outputs and returns a panda table with patients ids and associated snomed codes
    input_procedures, input_conditions, _, _ = Fetch_Output_Synthea(cohort_selection)
    ## We extract the SNOMED codes in the procedures and condition files ##
    print(
        "--------------Synthea output detected, extraction of patient snomed codes--------------"
    )
    if cohort_selection != "":
        print(
            f"------------------cohort of patients selected : {cohort_selection}------------------"
        )
    patients_data_procedures = input_procedures.loc[:, ("PATIENT", "CODE")]
    patients_data_procedures.rename(
        columns={"PATIENT": "pid", "CODE": "sctid"},
        inplace=True,
    )
    patients_data_conditions = input_conditions.loc[:, ("PATIENT", "CODE")]
    patients_data_conditions.rename(
        columns={"PATIENT": "pid", "CODE": "sctid"},
        inplace=True,
    )
    ## We include the additional information wanted ##
    print(
        "--------------Retrieving additional informations from synthea files...--------------"
    )
    patient_additional_infos = Retrieve_Additional_Infos_Synthea(cohort_selection)
    ## We make a unique table of all snomed codes ##
    patients_data = pd.concat(
        [patients_data_conditions, patients_data_procedures, patient_additional_infos],
        ignore_index=True,
    )
    ## We sort and reset index to have a pretty table ##
    patients_data = (
        patients_data.drop_duplicates().sort_values(by=["pid"]).reset_index(drop=True)
    )
    # synthea makes the mistake of generating sometimes the same pid between cohorts.
    # we add the first letters of the cohort name to the pid to avoid problems.
    if pid_rename and cohort_selection != "":
        patients_data["pid"] = cohort_selection[:2] + "-" + patients_data["pid"]
    ## We export the table created
    patients_data["sctid"] = patients_data["sctid"].astype(
        np.int64
    )  # to make sure we have integers
    patients_data.to_csv(
        "output.keepme/Output_data_synthetic_patients_" + cohort_selection + ".csv",
        index=False,
    )
    return patients_data


def EHR_files_creation(relation):
    gds = GraphDataScience(DB_URL, auth=(USER, PASSWORD))
    gds.run_cypher(
        f"""
            MATCH (n:Patient)
            WHERE n.pid = '{relation[0]}'
            MATCH (m:ObjectConcept)
            WHERE m.sctid= '{relation[1]}'
            MERGE (n)-[:EHR_CONTAINS]->(m)
            RETURN n
            """
    )


def Injection_Patients_Neo4j_db(cohort_selection: str = "", pid_rename: bool = False):
    # Retrieve using Table_of_Patients_From_Synthea() a clean table of patients data
    # (id and associated snomed codes) and inject it in the neo4j db
    gds.run_cypher(
        """
        CREATE CONSTRAINT personIdConstraint IF NOT EXISTS FOR (patient:Patient) REQUIRE patient.pid IS UNIQUE
        """
    )  # we need to create a constraint on the db to ensure no duplicata of patient is created

    ## We retrieve the snomed codes from the output of Synthea, using above built functions
    patients_data = Table_of_Patients_From_Synthea(cohort_selection, pid_rename)
    patients_unique_pid = patients_data["pid"].drop_duplicates().tolist()
    number_of_patients = len(patients_unique_pid)
    patients_data_list = patients_data.values.tolist()
    ## We now create new patients from the retrieved tables ##
    # note that the code allows to be run again, with same synthea output, and does not create problems
    # since in the cypher queries MERGE commands have been used instead of CREATE.
    print(
        f"------------------{number_of_patients} patients ready to be created, interacting now with neo4j db------------------"
    )
    if cohort_selection == "":
        for pid in tqdm(patients_unique_pid, desc="creation of patient nodes"):
            gds.run_cypher(
                f"""
                MERGE (p:Patient {{pid:'{pid}'}})
                RETURN p LIMIT 10
                """
            )
    # if patient is member of a cohort, we put this information in the node
    else:
        for pid in tqdm(patients_unique_pid, desc="creation of patient nodes"):
            gds.run_cypher(
                f"""
                MERGE (p:Patient {{pid:'{pid}',cohort:'{cohort_selection}'}})
                RETURN p LIMIT 10
                """
            )
    print("creation of the EHR files of each patient")
    p_umap(EHR_files_creation, patients_data_list)
    print(
        f"------------------{number_of_patients} patients and corresponding EHR successfully created------------------"
    )
