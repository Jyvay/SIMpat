from typing import List
from .consts import DB_URL, AUTH, USER, PASSWORD, graph, gds
from itertools import combinations, product
from .func import (
    rSubset,
    ProjectGraph,
    ProjectGraph_WithProperties,
    DirectProduct_Ordered,
    DirectProduct_Unordered,
    ProjectGraphWithPatients,
    ListOfNodesFromOnePatient,
    ListOfNodesFromAllPatients,
    Put_ShortestPath_table_into_good_form,
)
from neo4j import GraphDatabase, RoutingControl
from .ShortestPathPatient import ShortestPath_All, FindPatient_ids
from graphdatascience import GraphDataScience
from tqdm.auto import tqdm
from functools import partial
from p_tqdm import p_umap, p_map

# from .p_map import p_umap, p_map
import pandas as pd
import numpy as np
import math
import ast
import itertools
import multiprocessing as mp
import csv
import os


def InformationContent(id_node: int):
    query = gds.run_cypher(
        """
        MATCH (startNode)<-[:ISA*]-(descendant)
        WHERE id(startNode)=$node
        WITH DISTINCT descendant
        RETURN count(descendant) as NumberOfChildren
        """,
        {"node": id_node},
    )  # we retrieve the number of descendant
    nb_idesc = query["NumberOfChildren"][0] + 1
    nb_nodes = 453738  # computed via snomed query
    IC = nb_idesc / nb_nodes
    IC = -math.log(IC)
    return IC


def DescendantsNb(id_node: int):
    query = gds.run_cypher(
        """
        MATCH (startNode)<-[:ISA*]-(descendant)
        WHERE id(startNode)=$node
        WITH DISTINCT descendant
        RETURN count(descendant) as NumberOfChildren
        """,
        {"node": id_node},
    )  # we retrieve the number of descendant
    nb_idesc = query["NumberOfChildren"][0] + 1
    nb_nodes = 453738  # computed via snomed query
    IC = nb_idesc / nb_nodes
    return IC


def RetrieveAllSnomedConcepts():
    query = gds.run_cypher(
        """
        MATCH (source)-[r:ISA]-(target)
        WHERE NOT source:Patient
        RETURN DISTINCT id(source) AS id_nodes LIMIT 25
        """
    )  # Delete LIMIT 25 in above code to run on ALL snomed concepts. Let it if just want to test the function.
    id_nodes = query["id_nodes"].tolist()
    return id_nodes


def check_if_IC_property_exists():
    query = gds.run_cypher(
        """
        MATCH (n)-[ISA]->(m) WHERE n.cost_InverseIC IS NOT NULL RETURN count(n) as check_int
        """
    )
    if query["check_int"][0] > 0:
        return True
    else:
        return False


def InformationContent_WRITECOST(id_node: int):
    gds = GraphDataScience(DB_URL, auth=(USER, PASSWORD))
    data = pd.read_csv("output.keepme/Output_data_Nodes_DescendantsNb.csv")
    IC = data[data["node_id"] == id_node].reset_index()
    IC = IC["IC"][0]
    gds.run_cypher(
        """
        MATCH (Node)
        WHERE id(Node) = $node
        SET Node.cost_InverseIC=$IC
        RETURN Node
        """,
        {"node": id_node, "IC": IC},
    )


def InformationContentAll_WRITECOST():
    # Compute the information content
    List_of_nodes = gds.run_cypher(
        """
        MATCH (source)-[r]->(target)
        WHERE r:ISA
        WITH DISTINCT id(source) as id_nodes
        RETURN id_nodes
        """
    )
    List_of_nodes = List_of_nodes["id_nodes"].values.tolist()
    p_umap(InformationContent_WRITECOST, List_of_nodes)


def InformationContentRelationships_WRITECOST():
    # we define the information content of a relationship as the average of the information content of the two connected nodes
    gds.run_cypher(
        """
        MATCH (n)-[r]-(m)
        WHERE n.cost_InverseIC IS NOT NULL AND m.cost_InverseIC IS NOT NULL
        WITH n.cost_InverseIC AS cost_1, m.cost_InverseIC AS cost_2,n,m,r
        SET r.cost_AverageInverseIC = (cost_1 + cost_2)/2
        RETURN n LIMIT 25
        """
    )


def InformationContentRoleGroup_WRITECOST():
    gds.run_cypher(
        """
        MATCH (n:RoleGroup)-[r]-(m)
        WITH sum(m.cost_InverseIC) AS InverseIC,n
        SET n.cost_InverseIC = InverseIC
        RETURN InverseIC LIMIT 25
        """
    )


def InformationContent_RetrieveFromData(id_node: int):
    # given a node id, look into the data where all computed id nodes are to find its associated IC_value
    data = pd.read_csv("output.keepme/Output_data_Nodes_DescendantsNb.csv")
    Inverse_IC = data[data["node_id"] == id_node].reset_index(drop=True)
    Inverse_IC = Inverse_IC["IC"][0]
    return Inverse_IC


#### FOR THE ROLEGROUPS
# Rolegroups node are not included initially when executing the above code, but can be included very easily with a simple cypher query
# since Rolegroups IC are simply the sum of the IC of the nodes connected to it.
# the cypher query to execute is the following:
# MATCH (n:RoleGroup)-[r]-(m)
# WITH sum(m.cost_InverseIC) AS InverseIC,n
# SET n.cost_InverseIC = InverseIC
# RETURN InverseIC LIMIT 25


def Retrieve_Snomed_Labels():
    query = gds.run_cypher(
        """
        MATCH (source)-[r]->(target)
        WHERE NOT source:Patient
        RETURN DISTINCT labels(source) As Node_labels
        """
    )
    List_nodes = query["Node_labels"]
    List_nodes = List_nodes.tolist()
    Label_nodes = []
    for category in List_nodes:
        for label in category:
            Label_nodes.append(label)
    Label_nodes = list(dict.fromkeys(Label_nodes))
    return Label_nodes


def ShortestPath_Between2Nodes_WithICWeights(
    pair_of_nodes: list[int],
):  # Return a path of relations Ids given a projected graph, a source id and a target id
    gds = GraphDataScience(DB_URL, auth=(USER, PASSWORD))
    G = gds.graph.get(
        "noPatientUndirected_WithProperties"
    )  # we specify here the projected graph (if we change it need also to change the function ProjectGraph() in main() )
    result = gds.shortestPath.dijkstra.stream(
        G,
        sourceNode=pair_of_nodes[0],
        targetNode=pair_of_nodes[1],
        relationshipWeightProperty="cost_AverageInverseIC",
    )
    nodeIds = result["nodeIds"][0]
    query = gds.run_cypher(
        """
            MATCH (n)
            WHERE id(n) in $list
            WITH sum(n.cost_InverseIC) as IC_path
            RETURN IC_path
        """,
        {"list": nodeIds},
    )
    IC_path = query["IC_path"][0]
    return [nodeIds, IC_path]


def ShortestPath_All_WithICWeights(List_pair_of_nodes):
    # Return a path of relations ids given a projected graph, a source id and a target id.
    # The resulting path are in the order of the given array of pairs of nodes.
    All_pathes, IC_path = zip(
        *p_map(ShortestPath_Between2Nodes_WithICWeights, List_pair_of_nodes)
    )
    return All_pathes, IC_path


def ShortestPathInterPatients_WithPrecedentData(List_of_patients: list[int]):
    # Same function as the ones defined in ShortestPathPatient.py, put together and with precendent paths
    # computed in the first part of the project taken into account.
    print(
        f"--------------{len(List_of_patients)} Patient files detected, accessing now the nodes----------------"
    )
    List_of_nodes_all = ListOfNodesFromAllPatients(
        List_of_patients
    )  # we retrieve all uniques node from patients
    List_of_nodes = List_of_nodes_all["nodes_id"].values
    ProjectionExists = gds.graph.exists("noPatientUndirected")
    if ProjectionExists[
        "exists"
    ]:  # we check if the other kind of projected graph exist. If yes, we delete it to free some RAMs.
        G = gds.graph.get("noPatientUndirected")
        gds.graph.drop(G)

    ProjectionExists = gds.graph.exists(
        "noPatientUndirected_WithProperties"
    )  # We check if the projection exists
    if not ProjectionExists["exists"]:
        print(
            "----------------------------Projecting the Snomed graph---------------------------"
        )
        ProjectGraph_WithProperties()
    print(
        "----------------------------------Pairing the nodes-------------------------------"
    )
    List_of_pairs = rSubset(List_of_nodes, 2)
    List_of_pairs = List_of_pairs.tolist()
    FileExists = os.path.isfile(
        "output.keepme/Output_data_ShortestPath_Node_to_Node.csv"
    )
    # we check if the file from the first part of the project is present.
    # if yes, we can spare a lot of computation by using it.
    if FileExists:
        print(
            "---------Retrieving and comparison of already computed pairs of nodes-------------"
        )
        ShortestPath = pd.read_csv(
            "output.keepme/Output_data_ShortestPath_Node_to_Node.csv"
        )
        List_of_pairs_already_treated = ShortestPath["pairs_of_nodes"].values.tolist()
        List_of_pairs_already_treated = [
            list(ast.literal_eval(pairs)) for pairs in List_of_pairs_already_treated
        ]
        List_of_paths_already_treated = ShortestPath["paths"].values.tolist()
        List_of_paths_already_treated = [
            list(ast.literal_eval(paths)) for paths in List_of_paths_already_treated
        ]
        List_of_pairs_min = [
            pairs
            for pairs in List_of_pairs
            if pairs not in List_of_pairs_already_treated
        ]  # for pairs already treated, no need to include them.
        print(
            f"----------------{len(List_of_pairs_min)} pairs to be computed instead of {len(List_of_pairs)}----------------------"
        )
    else:
        # if file is not present, no worries but we have to compute all pairs.
        List_of_pairs_min = List_of_pairs
    print(
        f"--------------------Computing Shortest Paths of {len(List_of_pairs_min)} pairs of nodes-----------------"
    )
    paths = ShortestPath_All(
        List_of_pairs_min
    )  # we compute the shortest path for all pairs of unique concepts
    result = pd.DataFrame({"pairs_of_nodes": List_of_pairs_min, "paths": paths})
    print(
        "--------------------Operation succeeded, all paths have been found----------------"
    )
    if FileExists:
        data = pd.DataFrame(
            {
                "pairs_of_nodes": List_of_pairs_already_treated,
                "paths": List_of_paths_already_treated,
            }
        )
        output = pd.concat([data, result], ignore_index=True)
    else:
        output = result
    output["distance"] = output["paths"].apply(len)
    return output


def MinDistance_Nodes_to_Patients(table_shortest_path, List_of_patients: list[int]):
    # given a table of shortest path between nodes and a table of patients with their nodes
    # return a table with the shortest path from one patient to another
    MinDistance_Nodes_Patients = gds.run_cypher(
        "MATCH (n:Patient)-[r]-(m) RETURN DISTINCT id(m) AS id_nodes"
    )
    MinDistance_Nodes_Patients.sort_values(by="id_nodes", inplace=True)
    MinDistance_Nodes_Patients.reset_index(drop=True, inplace=True)
    dataframe_to_append = [MinDistance_Nodes_Patients]
    for patient in tqdm(
        List_of_patients, "Finding the minimum distance from nodes to patients"
    ):
        dist_node_patient = table_shortest_path
        query = gds.run_cypher(
            f"MATCH (n:Patient)--(m) WHERE id(n) = {patient} RETURN DISTINCT id(m) AS nodes"
        )
        nodes_patient = query["nodes"].values.tolist()
        mask = dist_node_patient["second_node"].isin(nodes_patient)
        dist_node_patient = dist_node_patient[mask]
        dist_node_patient = (
            dist_node_patient.loc[:, "distance"].groupby(level="first_node").min()
        )
        dist_node_patient = dist_node_patient.to_frame(name=str(patient))
        dist_node_patient.reset_index(drop=True, inplace=True)
        dataframe_to_append.append(dist_node_patient)
    MinDistance_Nodes_Patients = pd.concat(dataframe_to_append, axis=1)
    MinDistance_Nodes_Patients.set_index("id_nodes", inplace=True)

    return MinDistance_Nodes_Patients


def Table_for_Metric_Bags_of_Findings(table_shortest_path, List_of_patients: list[int]):
    # Give an indicator matrix with rows : nodes_id, cols: patients_id, entries : 1 if nodes is
    # in the patient file, 0 if the node isn't.
    # We use a (put in the good form) shortest path table as a framework for our indicator matrix
    # we do not use its values.
    table_shortest_path.reset_index(inplace=True)
    table_shortest_path["distance"] = (
        table_shortest_path["first_node"] == table_shortest_path["second_node"]
    )
    table_shortest_path.set_index("first_node", inplace=True)
    Does_Nodes_belong_to_Patients = gds.run_cypher(
        "MATCH (n:Patient)-[r]-(m) RETURN DISTINCT id(m) AS id_nodes"
    )
    Does_Nodes_belong_to_Patients.sort_values(by="id_nodes", inplace=True)
    Does_Nodes_belong_to_Patients.reset_index(drop=True, inplace=True)
    dataframe_to_append = [Does_Nodes_belong_to_Patients]
    for patient in tqdm(
        List_of_patients,
        "Constructing an indicator matrix for the metric Bags of Findings",
    ):
        dist_node_patient = table_shortest_path
        nodes_patient = pd.DataFrame(
            columns=["List_of_nodes", "Appears_in_patient"],
            data=ListOfNodesFromOnePatient(patient),
        )
        nodes_patient = nodes_patient["List_of_nodes"].values.tolist()
        mask = dist_node_patient["second_node"].isin(nodes_patient)
        dist_node_patient = dist_node_patient[mask]
        dist_node_patient = (
            dist_node_patient.loc[:, "distance"].groupby(level="first_node").sum()
        )
        dist_node_patient = dist_node_patient.to_frame(name=str(patient))
        dist_node_patient.reset_index(drop=True, inplace=True)
        dataframe_to_append.append(dist_node_patient)
    Does_Nodes_belong_to_Patients = pd.concat(dataframe_to_append, axis=1)
    Does_Nodes_belong_to_Patients.set_index("id_nodes", inplace=True)
    print(Does_Nodes_belong_to_Patients)
    return Does_Nodes_belong_to_Patients

def All_Metrics_optimized(List_of_patients: list[int]):
    # we first retrieve our data
    shortest_path_data = pd.read_csv(
        "output.keepme/Output_data_ShortestPath_Node_to_Node.csv"
    )
    shortest_path_data = Put_ShortestPath_table_into_good_form(shortest_path_data)
    shortest_path_data_ICWeights = pd.read_csv(
        "output.keepme/Output_data_ShortestPath_Node_to_Node_ICWeights.csv"
    )
    shortest_path_data_ICWeights = Put_ShortestPath_table_into_good_form(
        shortest_path_data_ICWeights
    )

    # we retrieve the information content of the nodes
    IC_nodes = pd.read_csv("output.keepme/Output_data_Nodes_IC.csv")

    # we then find the min. distance of each node to the patients set of nodes..
    print("-------Computing the min distances between nodes and patients-------")

    Table_Bags_of_Findings = Table_for_Metric_Bags_of_Findings(
        shortest_path_data, List_of_patients
    )

    Dist_node_to_patient = MinDistance_Nodes_to_Patients(
        shortest_path_data, List_of_patients
    )
    Dist_node_to_patient_ICPath = MinDistance_Nodes_to_Patients(
        shortest_path_data_ICWeights, List_of_patients
    )

    print("Exporting the tables...")
    Table_Bags_of_Findings.to_csv("output.keepme/Table_Bags_of_Findings.csv")
    Dist_node_to_patient.to_csv("output.keepme/Dist_node_to_patient.csv")
    Dist_node_to_patient_ICPath.to_csv("output.keepme/Dist_node_to_patient_ICPath.csv")

    # Transforming pandas arrays to numpy ones
    ones = Table_Bags_of_Findings >= 0
    ones = ones.astype(int)
    ones = np.transpose(ones.to_numpy())
    Table_Bags_of_Findings = Table_Bags_of_Findings.to_numpy()
    Dist_node_to_patient = Dist_node_to_patient.to_numpy()
    Dist_node_to_patient_ICPath = Dist_node_to_patient_ICPath.to_numpy()

    print("--------We now compute the metric 'Bags of Findings'--------")
    Table_Bags_of_Findings_tr = np.transpose(Table_Bags_of_Findings)
    Intersection_Matrix = np.dot(Table_Bags_of_Findings_tr, Table_Bags_of_Findings)

    card = np.dot(ones, Table_Bags_of_Findings)
    Union_Matrix = np.transpose(card) + card - Intersection_Matrix
    Metric_BagsOfFindings = 1 - np.divide(Intersection_Matrix, Union_Matrix)

    # We compute the metrics AverageLinks, AverageLinks with IC Weights, IC Path
    IC_list = IC_nodes["IC"].values.tolist()
    CoeffICWeights_diag_matrix = np.diag(IC_list)
    Sum_of_Children = gds.run_cypher(
        "MATCH (n:Patient)-[r]->(m) WITH id(n) as patient_id, count(m)^(-1) as children_nb RETURN patient_id, children_nb ORDER BY patient_id"
    )
    Sum_of_Children_list = Sum_of_Children["children_nb"].values.tolist()
    SumCoeff_diag_matrix = np.diag(Sum_of_Children_list)
    Sum_of_ICWeights = gds.run_cypher(
        "MATCH (n:Patient)-[r]->(m) WITH id(n) as patient_id, 1/sum(-log(m.cost_InverseIC)) as sum_IC RETURN patient_id, sum_IC ORDER BY patient_id"
    )
    Sum_of_ICWeights_list = Sum_of_ICWeights["sum_IC"].values.tolist()
    SumCoeffICWeights_diag_matrix = np.diag(Sum_of_ICWeights_list)
    Indicator_Matrix = Table_Bags_of_Findings_tr
    print("--------We now compute the metric 'Average Links'--------")
    Metric_AverageLinks = np.dot(SumCoeff_diag_matrix, Indicator_Matrix)
    Metric_AverageLinks = np.dot(Metric_AverageLinks, Dist_node_to_patient)
    print("------We now compute the metric 'Average Links' with IC coefficients-----")
    Coeff = np.dot(SumCoeffICWeights_diag_matrix, Indicator_Matrix)
    Coeff = np.dot(Coeff, CoeffICWeights_diag_matrix)

    Metric_AverageLinks_ICWeights = np.dot(Coeff, Dist_node_to_patient)

    print("------We now compute the metric 'IC Path' with IC coefficients-----")
    Metric_ICPath = np.dot(Coeff, Dist_node_to_patient_ICPath)
    print("--------------------We now symmetrize each metric------------------")

    Metric_BagsOfFindings = (
        Metric_BagsOfFindings + np.transpose(Metric_BagsOfFindings)
    ) * 0.5
    Metric_AverageLinks = (
        Metric_AverageLinks + np.transpose(Metric_AverageLinks)
    ) * 0.5
    Metric_AverageLinks_ICWeights = (
        Metric_AverageLinks_ICWeights + np.transpose(Metric_AverageLinks_ICWeights)
    ) * 0.5
    Metric_ICPath = (Metric_ICPath + np.transpose(Metric_ICPath)) * 0.5
    print("-----All metrics have been sucessfully generated------")
    return (
        Metric_BagsOfFindings,
        Metric_AverageLinks,
        Metric_AverageLinks_ICWeights,
        Metric_ICPath,
    )


def ShortestPathInterPatients_ICWeights(List_of_patients: list[int]):
    # Same function as the ones defined in ShortestPathPatient.py, put together and with precendent paths
    # computed in the first part of the project taken into account.
    print(
        f"--------------{len(List_of_patients)} Patient files detected, accessing now the nodes----------------"
    )
    List_of_nodes_all = ListOfNodesFromAllPatients(
        List_of_patients
    )  # we retrieve all uniques node from patients
    List_of_nodes = List_of_nodes_all["nodes_id"].values
    ProjectionExists = gds.graph.exists("noPatientUndirected")
    if ProjectionExists[
        "exists"
    ]:  # we check if the other kind of projected graph exist. If yes, we delete it to free some RAMs.
        G = gds.graph.get("noPatientUndirected")
        gds.graph.drop(G)

    ProjectionExists = gds.graph.exists(
        "noPatientUndirected_WithProperties"
    )  # We check if the projection exists
    if not ProjectionExists["exists"]:
        print(
            "----------------------------Projecting the Snomed graph---------------------------"
        )
        ProjectGraph_WithProperties()
    print(
        "----------------------------------Pairing the nodes-------------------------------"
    )
    List_of_pairs = rSubset(List_of_nodes, 2)
    List_of_pairs = List_of_pairs.tolist()
    FileExists = os.path.isfile(
        "output.keepme/Output_data_ShortestPath_Node_to_Node_ICWeights.csv"
    )
    # we check if the file from the first part of the project is present.
    # if yes, we can spare a lot of computation by using it.
    if FileExists:
        print(
            "---------Retrieving and comparison of already computed pairs of nodes-------------"
        )
        ShortestPath_data = pd.read_csv(
            "output.keepme/Output_data_ShortestPath_Node_to_Node_ICWeights.csv"
        )
        List_of_pairs_already_treated = ShortestPath_data[
            "pairs_of_nodes"
        ].values.tolist()
        List_of_pairs_already_treated = [
            list(ast.literal_eval(pairs)) for pairs in List_of_pairs_already_treated
        ]
        List_of_paths_already_treated = ShortestPath_data["paths"].values.tolist()
        List_of_paths_already_treated = [
            list(ast.literal_eval(paths)) for paths in List_of_paths_already_treated
        ]
        List_of_pairs_min = [
            pairs
            for pairs in List_of_pairs
            if pairs not in List_of_pairs_already_treated
        ]  # for pairs already treated, no need to include them.
        print(
            f"----------------{len(List_of_pairs_min)} pairs to be computed instead of {len(List_of_pairs)}----------------------"
        )
    else:
        # if file is not present, no worries but we have to compute all pairs.
        List_of_pairs_min = List_of_pairs
    print(
        f"--------------------Computing Shortest Paths of {len(List_of_pairs_min)} pairs of nodes-----------------"
    )
    # we compute the shortest path for all pairs of unique concepts
    paths, distance_ic = ShortestPath_All_WithICWeights(List_of_pairs_min)
    result = pd.DataFrame(
        {"pairs_of_nodes": List_of_pairs_min, "paths": paths, "distance": distance_ic}
    )
    print(
        "--------------------Operation succeeded, all paths have been found----------------"
    )
    if FileExists:
        output = pd.concat([ShortestPath_data, result], ignore_index=True)
    else:
        output = result
    return output