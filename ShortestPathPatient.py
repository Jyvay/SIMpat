from .consts import DB_URL, AUTH, USER, PASSWORD, graph, gds

from .func import (
    rSubset_iter,
    ProjectGraph,
    DeleteAllCounter,
    ProjectGraphWithPatients,
    FindPatient_ids,
)
from itertools import combinations
from neo4j import GraphDatabase, RoutingControl
from graphdatascience import GraphDataScience
from p_tqdm import p_umap, p_map

# from .p_map import p_umap, p_map
import pandas as pd
import numpy as np
import time
import multiprocessing as mp


def RetrievePairofNodes(Patient_id):
    # given a patient id, return its pairs of nodes along with the patient id.
    gds = GraphDataScience(DB_URL, auth=(USER, PASSWORD))
    query1 = gds.run_cypher(
        "MATCH (n:Patient)--(m) WHERE id(n)=$id RETURN id(m)", {"id": Patient_id}
    )
    id_nodes = query1["id(m)"]
    Pairs_of_nodes = rSubset_iter(id_nodes, 2)
    output = pd.DataFrame({"pairs_of_nodes": Pairs_of_nodes})
    output["pairs_of_nodes"] = output["pairs_of_nodes"].apply(
        lambda x: tuple(sorted(x))
    )  # we sort the pairs of nodes (minimum of the two first). We do NOT want (node_1,node_2) AND (node_1,node_2) in our array, only (node_1,node_2) (with counter updated accordingly). Else ShortestPath_All() will compute the path twice.
    output["counter"] = 1
    output["appears_in_patient_id"] = output.apply(lambda x: [Patient_id], axis=1)
    output = output.values.tolist()
    return output


def PrepareNodes(List_of_patients):
    # From a list of patients id returns all pairs of unique nodes. Use multiprocessing for a faster treatment.
    All_pairs_messy = p_map(RetrievePairofNodes, List_of_patients)
    All_pairs_flat = [
        item for sublist in All_pairs_messy for item in sublist
    ]  # we unwind the resulting array (it is nested)
    All_pairs = pd.DataFrame(
        data=All_pairs_flat,
        columns=["pairs_of_nodes", "counter", "appears_in_patient_id"],
    )
    Unique_pairs = All_pairs.groupby(
        "pairs_of_nodes"
    ).sum()  # we want distinct pairs of nodes
    Unique_pairs = Unique_pairs.reset_index(
        names=["pairs_of_nodes"]
    )  # we reset the index (groupby has set the index to "pairs_of_nodes")
    return Unique_pairs


def ShortestPath_Between2Nodes(
    pair_of_nodes,
):  # Return a path of relations Ids given a projected graph, a source id and a target id
    # use the property "cost" to avoid certain nodes.
    gds = GraphDataScience(DB_URL, auth=(USER, PASSWORD))
    ProjectionExists = gds.graph.exists("noPatientUndirected_WithProperties")
    # We check if the user projected a graph with the cost property. If yes, we will use it for the gds module shortest path
    if ProjectionExists["exists"]:
        G = gds.graph.get("noPatientUndirected_WithProperties")
        result = gds.shortestPath.dijkstra.stream(
            G, sourceNode=pair_of_nodes[0], targetNode=pair_of_nodes[1]
        )
    else:
        # Else, we assume the user called the standard projection coming with this project.
        G = gds.graph.get("noPatientUndirected")
        result = gds.shortestPath.dijkstra.stream(
            G,
            sourceNode=pair_of_nodes[0],
            targetNode=pair_of_nodes[1],
        )
    nodeIds = result["nodeIds"][
        0
    ]  # we select in the panda dataframe the entry corresponding to the path of nodes
    query = gds.run_cypher(
        """
                MATCH (n)-[r]->(m)
                WHERE id(n) IN $List_id AND id(m) IN $List_id
                RETURN DISTINCT id(r)
                """,
        {"List_id": nodeIds},
    )  # we want a list of relations id, not of nodes.
    id_relations = query["id(r)"]
    id_relations = id_relations.values.tolist()
    return id_relations


def ShortestPath_All(
    List_pair_of_nodes,
):  # Return a path of relations ids given a projected graph, a source id and a target id. The resulting path are in the order of the given array of pairs of nodes.
    All_pathes = p_map(ShortestPath_Between2Nodes, List_pair_of_nodes)
    return All_pathes


def ShortestPath_Between2Nodes_alt_PathOfNodes(
    pair_of_nodes: list[int],
):  # Return a path of nodes Ids given a projected graph, a source id and a target id
    ProjectionExists = gds.graph.exists("noPatientUndirected_WithProperties")
    if ProjectionExists[
        "exists"
    ]:  # we check if the other kind of projected graph exist. If yes, we delete it to free some RAMs.
        G = gds.graph.get("noPatientUndirected_WithProperties")
        gds.graph.drop(G)
    ProjectionExists = gds.graph.exists(
        "noPatientUndirected"
    )  # We check if the projection exists
    if not ProjectionExists["exists"]:
        print(
            "----------------------------Projecting the Snomed graph--------------------------------"
        )
        ProjectGraph()
    G = gds.graph.get("noPatientUndirected")
    result = gds.shortestPath.dijkstra.stream(
        G, sourceNode=pair_of_nodes[0], targetNode=pair_of_nodes[1]
    )
    nodeIds = result["nodeIds"][
        0
    ]  # we select in the panda dataframe the entry corresponding to the path of nodes
    # query = gds.run_cypher(
    #     """
    #             MATCH (n)-[r]->(m)
    #             WHERE id(n) IN $List_id AND id(m) IN $List_id
    #             RETURN DISTINCT id(r)
    #             """,
    #     {"List_id": nodeIds},
    # )  # we want a list of relations id, not of nodes.
    # id_relations = query["id(r)"]
    # id_relations = id_relations.values.tolist()
    return nodeIds


def PrepareRelations(
    Global_array,
):  # takes as input the array AFTER treatment by ShortestPath_all with columns "counter", "appears_in_patient_id" and "paths". Return as output an array with "relations_id", "counter" and "appears_in_patient_id". It is this array that will get passed as input in the ManageNeo4jCounter.
    unnested_data = []
    for _, row in Global_array.iterrows():
        counter = row["counter"]
        appears_in_patient_id = row["appears_in_patient_id"]
        paths = row["paths"]
        for relation_id in paths:
            unnested_data.append(
                {
                    "counter": counter,
                    "appears_in_patient_id": appears_in_patient_id,
                    "relation_id": relation_id,
                }
            )

    unnested_df = pd.DataFrame(unnested_data)
    unnested_df = unnested_df.groupby("relation_id").sum()
    unnested_df["appears_in_patient_id"] = unnested_df["appears_in_patient_id"].apply(
        lambda x: list(set(x))
    )  # we eliminate in the column "appears_in_patient_id" all values that appear multiple times
    List_of_relations = unnested_df.reset_index(names=["relation_id"])
    return List_of_relations


def ManageNeo4jCounter(
    Relations_to_update,
):  # Needs an array in input with relations to add a counter to, the counter number to add and a list of patient where this relations occurs
    # we assume that the array in input has the format 'relation_id','counter','appears_in_patient_id' as columns names
    List_of_relations_id = Relations_to_update["relation_id"]
    case_statement_counter = "\n".join(
        [
            f"WHEN {row['relation_id']} THEN {row['counter']}"
            for _, row in Relations_to_update.iterrows()
        ]
    )  # create a string for the cypher query. Used with the 'CASE' in the query, to associate the correct counter to its associated relation
    case_statement_patient = "\n".join(
        [
            f"WHEN {row['relation_id']} THEN {row['appears_in_patient_id']}"
            for _, row in Relations_to_update.iterrows()
        ]
    )
    gds.run_cypher(
        f"""
            MATCH (n)-[r]->(m)
            WHERE id(r) IN $Relations_id
            WITH n, m, r, id(r) AS Id_r //need to use WITH or else don't work
            SET r.counter =
                CASE Id_r
                    {case_statement_counter}
                END
            SET r.patient =
                CASE Id_r
                    {case_statement_patient}
                END
            RETURN r
            """,
        {"Relations_id": List_of_relations_id},
    )


def ShortestPath_BetweenPatients_And_Counter_Update():
    # Use all the functions defined above in a sequential order to compute the Shortest Path
    # between patients and update a counter to each relations used.
    ProjectionExists = gds.graph.exists("noPatientUndirected_WithProperties")
    if ProjectionExists[
        "exists"
    ]:  # we check if the other kind of projected graph exist. If yes, we delete it to free some RAMs.
        G = gds.graph.get("noPatientUndirected_WithProperties")
        gds.graph.drop(G)
    ProjectionExists = gds.graph.exists(
        "noPatientUndirected"
    )  # We check if the projection exists
    if not ProjectionExists["exists"]:
        print(
            "----------------------------Projecting the Snomed graph---------------------------"
        )
        ProjectGraph()
    List_of_patients = FindPatient_ids()  # We retrieve the list of all patients
    print(
        f"----------------{len(List_of_patients)} Patient files detected, accessing now the nodes----------------"
    )
    Global_array_unique_nodes = PrepareNodes(
        List_of_patients
    )  # We retrieve the pairs of unique nodes along with the number of times they occured and the patients they appeared in
    Pairs_of_nodes = Global_array_unique_nodes["pairs_of_nodes"].values.tolist()
    print(
        f"-----------------------Preparing to treat {len(Pairs_of_nodes)} pairs of nodes-----------------------"
    )
    paths = ShortestPath_All(
        Pairs_of_nodes
    )  # We apply Djikstra to all distinct pairs of nodes. Return an array of paths of relationships id.
    Global_array_unique_nodes[
        "paths"
    ] = paths  # We add the paths to their corresponding pair of nodes.
    print(
        "-----------------------Preparing the array of Relations...-----------------------"
    )
    List_of_relations = PrepareRelations(
        Global_array_unique_nodes
    )  # We prepare the final array with distinct relations along with their corresponding counter and the patients they appeared in
    print(
        f"-------------We are now updating {len(List_of_relations)} relations in the Neo4j database--------------"
    )
    ManageNeo4jCounter(List_of_relations)  # We update the Neo4j database
    print(
        "--------------Operation succeeded, all relations have been updated---------------"
    )
