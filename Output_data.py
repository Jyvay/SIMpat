from typing import List
from .consts import gds
from .Metrics import (
    RetrieveAllSnomedConcepts,
    InformationContent,
    DescendantsNb,    
    InformationContentAll_WRITECOST,
    InformationContentRelationships_WRITECOST,
    InformationContentRoleGroup_WRITECOST,
    ShortestPathInterPatients_WithPrecedentData,
    ShortestPathInterPatients_ICWeights,
    All_Metrics_optimized,
)
from .ShortestPathPatient import (
    FindPatient_ids,
    ShortestPath_Between2Nodes_alt_PathOfNodes,
    PrepareNodes,
)
from .func import ListOfReference_Patients_Nodes_pid,check_if_IC_property_exists
from p_tqdm import p_umap, p_map
import pandas as pd
import ast
import os


def Output_data_ListOfReference_Patients_Nodes_pid():
    # return a csv file with patient id from the Neo4j db, their pid, and the list of nodes id attached to each patient
    List_of_patients = FindPatient_ids()
    output = ListOfReference_Patients_Nodes_pid(List_of_patients)
    output.to_csv(
        "output.keepme/Output_data_ListOfReference_Patients_Nodes_pid.csv", index=False
    )
    return output


def Output_data_ShortestPath_patients():
    # From a list of patients, take all pairs of nodes and find shortest path between them using the Snomed Graph
    List_of_patients = FindPatient_ids()
    List_pair_of_nodes = PrepareNodes(List_of_patients)
    List_pair_of_nodes = List_pair_of_nodes["pairs_of_nodes"].to_list()
    All_pathes = p_map(ShortestPath_Between2Nodes_alt_PathOfNodes, List_pair_of_nodes)
    output = pd.DataFrame({"pairs_of_node": List_pair_of_nodes, "path": All_pathes})
    output["distance"] = output["path"].apply(len)
    output.to_csv("output.keepme/Output_data_ShortestPath_patients", index=False)
    return output


def Output_data_Nodes_IC():
    # return all unique nodes of patients along with their Information Content
    query = gds.run_cypher(
        "MATCH (n:Patient)-[r]-(m) RETURN DISTINCT id(m) as id_nodes"
    )
    List_of_nodes = query["id_nodes"].values.tolist()
    Information_Content = p_map(InformationContent, List_of_nodes)

    output = pd.DataFrame({"id_nodes": List_of_nodes, "IC": Information_Content})
    output.sort_values(by="id_nodes", inplace=True)
    output.to_csv("output.keepme/Output_data_Nodes_IC.csv", index=False)
    return output


def Output_data_Nodes_DescendantsNb():
    List_of_nodes = RetrieveAllSnomedConcepts()
    Information_Content = p_map(DescendantsNb, List_of_nodes)

    output = pd.DataFrame({"node_id": List_of_nodes, "IC": Information_Content})
    output.to_csv("output.keepme/Output_data_Nodes_DescendantsNb", index=False)
    return output


def Output_data_ShortestPath_Node_to_Node():
    List_of_patients = FindPatient_ids()
    output = ShortestPathInterPatients_WithPrecedentData(List_of_patients)
    output.to_csv(
        "output.keepme/Output_data_ShortestPath_Node_to_Node.csv",
        index=False,
    )

    return output


def Output_data_ShortestPath_Node_to_Node_ICWeights():
    List_of_patients = FindPatient_ids()
    output = ShortestPathInterPatients_ICWeights(List_of_patients)
    output.to_csv(
        "output.keepme/Output_data_ShortestPath_Node_to_Node_ICWeights.csv", index=False
    )
    return output


def Output_data_AllPatientsWithFourMetrics():
    List_of_patients = FindPatient_ids()
    Output_data_ListOfReference_Patients_Nodes_pid()
    Output_data_Nodes_IC()
    Patient_pid = pd.read_csv(
        "output.keepme/Output_data_ListOfReference_Patients_Nodes_pid.csv"
    )
    Patient_pid = Patient_pid["pid"].values.tolist()

    (
        Metric_BagsOfFindings,
        Metric_AverageLinks,
        Metric_AverageLinks_ICWeights,
        Metric_ICPath,
    ) = All_Metrics_optimized(List_of_patients)

    print("--------Labeling of the data--------")
    Metric_BagsOfFindings = pd.DataFrame(
        data=Metric_BagsOfFindings, index=Patient_pid, columns=Patient_pid
    )
    Metric_AverageLinks = pd.DataFrame(
        data=Metric_AverageLinks, index=Patient_pid, columns=Patient_pid
    )
    Metric_AverageLinks_ICWeights = pd.DataFrame(
        data=Metric_AverageLinks_ICWeights, index=Patient_pid, columns=Patient_pid
    )
    Metric_ICPath = pd.DataFrame(
        data=Metric_ICPath, index=Patient_pid, columns=Patient_pid
    )
    print("--------Exporting now the Metrics--------")
    Metric_BagsOfFindings.to_csv(
        "output.keepme/Metric_BagsOfFindings.csv", index_label="pid"
    )
    Metric_AverageLinks.to_csv(
        "output.keepme/Metric_AverageLinks.csv", index_label="pid"
    )
    Metric_AverageLinks_ICWeights.to_csv(
        "output.keepme/Metric_AverageLinks_ICWeights.csv", index_label="pid"
    )
    Metric_ICPath.to_csv("output.keepme/Metric_ICPath.csv", index_label="pid")


def Metric_Similarity_Assistant():
    # from a neo4j db with patients, run all needed functions to have as output the metrics.
    # has built-in prompt that ask user if he wants to continue, and from where.
    functions_to_run = [
        "Output_data_ShortestPath_Node_to_Node",
        "Output_data_ShortestPath_Node_to_Node_ICWeights",
        "Output_data_AllPatientsWithFourMetrics",
    ]
    yes_choices = ["yes", "y"]
    no_choices = ["no", "n"]
    print("Welcome to the similarity metrics generator helper")
    user_input = input(
        "Is this the first time do you launch this program (yes/no)?\nYour answer: "
    )
    if user_input.lower() in yes_choices:
        print(
            "This assistant will guide you through the steps needed to obtain a csv file with metrics of similarity."
        )
        print(
            "Please keep in mind to look at the readme file first before running this program."
        )
        print(
            "To generate the metrics, preliminary files need first to be generated. These files do take quite a long time\nto generate. At the end of each step, this program will ask you if you want to continue with the next step or terminate.\nIf you terminate the program, please note at which step you currently are. The next time you run this program, you can\nskip this introductory information and directly select the next step."
        )

    print("The functions that will be run by this assistant are the following:")
    iteration = 1
    for function_name in functions_to_run:
        print(f"{iteration} : " + function_name + "()")
        iteration += 1
    user_input = int(
        input(
            "Please select the step you want to resume (1-3). If you want to launch from the beginning, answer 1.\nYour answer: "
        )
    )
    if user_input == 1:
        print(f"Launching module for step {user_input}...")
        ## We check if the file already exists ##
        FileExists = os.path.isfile(
            "output.keepme/" + f"{functions_to_run[user_input-1]}" + ".csv"
        )
        if FileExists:
            user_input = input(
                "The file "
                + f"{functions_to_run[user_input-1]}"
                + ".csv is already existing and will be overwritten. Do you wish to continue (yes/no)?\nYour answer: "
            )
        else:
            user_input = "yes"
        if user_input.lower() in yes_choices:
            # User chose to run this step
            Output_data_ShortestPath_Node_to_Node()
        user_input = 1
        user_input = input(
            "The file "
            + f"{functions_to_run[user_input-1]}"
            + ".csv is sucessfully generated. Do you wish to continue with step 2/3 (yes/no)?\nYour answer: "
        )
        if user_input.lower() in yes_choices:
            user_input = 2
        else:
            user_input = 1  # we put again 1 so the 'else' part at the ends knows at wich step we stopped.
    if user_input == 2:
        print(f"Launching module for step {user_input}...")
        FileExists = os.path.isfile(
            "output.keepme/" + f"{functions_to_run[user_input-1]}" + ".csv"
        )
        PreviousFileExists = os.path.isfile(
            "output.keepme/" + f"{functions_to_run[user_input-2]}" + ".csv"
        )
        if FileExists:
            user_input = input(
                "The file "
                + f"{functions_to_run[user_input-1]}"
                + ".csv is already existing and will be overwritten. Do you wish to continue (yes/no)?\nYour answer: "
            )
        else:
            user_input = "yes"
        if (not PreviousFileExists) and (user_input != "no"):
            user_input = 2
            user_input = input(
                "The file "
                + f"{functions_to_run[user_input-2]}"
                + ".csv from previous step is missing. Do you wish to continue (yes/no)?\nYour answer: "
            )
        check_value = check_if_IC_property_exists()
        if check_value:
            Output_data_ShortestPath_Node_to_Node_ICWeights()
        else:
            user_input = input(
                "The assistant detected that the present neo4j database has no property cost_InverseIC\nwhich is needed for this step. Do you want to solve this issue (yes/no)?\nYour answer:"
            )
            if user_input.lower() in yes_choices:
                print("Launching the Information Content module...")
                Module_Information_Content()
                print(
                    "The information content of all nodes in the neo4j database have been generated. The assistant will now proceed with step 4."
                )
                Output_data_ShortestPath_Node_to_Node_ICWeights()
            else:
                print(
                    "The property cost_InverseIC is needed and missing. This assistant will quit now..."
                )
                quit()
        user_input = 3
        user_input = input(
            "The file "
            + f"{functions_to_run[user_input-1]}"
            + ".csv is sucessfully generated. Do you wish to continue with the final step (yes/no)?\nYour answer: "
        )
        if user_input.lower() in yes_choices:
            user_input = 3
        else:
            user_input = 2

    if user_input == 3:
        print(f"Launching module for step {user_input}...")
        FileExists = os.path.isfile(
            "output.keepme/" + f"{functions_to_run[user_input-1]}" + ".csv"
        )
        PreviousFileExists = os.path.isfile(
            "output.keepme/" + f"{functions_to_run[user_input-2]}" + ".csv"
        )
        if FileExists:
            user_input = input(
                "The file "
                + f"{functions_to_run[user_input-1]}"
                + ".csv is already existing and will be overwritten. Do you wish to continue (yes/no)?\nYour answer: "
            )
        else:
            user_input = "yes"
        if (not PreviousFileExists) and (user_input != "no"):
            user_input = 3
            user_input = input(
                "The file "
                + f"{functions_to_run[user_input-2]}"
                + ".csv from previous step is missing. Do you wish to continue (yes/no)?\nYour answer: "
            )
        if user_input.lower() in yes_choices:
            Output_data_AllPatientsWithFourMetrics()
        user_input = 3
        print(
            "The file "
            + f"{functions_to_run[user_input-1]}"
            + ".csv is sucessfully generated. Similarity metrics are now available, yeah!"
        )
    else:
        print(
            "You completed step "
            + f"{user_input}, succesfully generated "
            + f"{functions_to_run[user_input-1]}"
            + ".csv and chose to stop.\nThis assistant will quit now..."
        )
        quit()

def Module_Information_Content():
    print("Actually computing nodes information content from neo4j database")
    Output_data_Nodes_DescendantsNb()
    print(
        "Output_data_Nodes_DescendantsNb.csv has been generated. Accessing the neo4j database for property writing..."
    )
    InformationContentAll_WRITECOST()
    InformationContentRelationships_WRITECOST()
    InformationContentRoleGroup_WRITECOST()
    print("Finished.")
