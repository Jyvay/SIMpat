import itertools as iters
from math import factorial
from neo4j import GraphDatabase, RoutingControl
from graphdatascience import GraphDataScience
from .consts import DB_URL, AUTH, gds, graph
from operator import itemgetter
from p_tqdm import p_umap, p_map
import pandas as pd
import numpy as np
import csv
import ast
from tqdm.auto import tqdm


def FindPatient_ids():  # return all patients id
    source = gds.run_cypher("MATCH (n:Patient) RETURN id(n)")
    result = source["id(n)"]
    return result


def rSubset(array: np.ndarray, r: int = 3) -> np.ndarray:
    n_items = array.shape[0]
    num_combinations = factorial(n_items) // (factorial(n_items - r) * factorial(r))
    combination_idx = np.fromiter(
        iters.chain.from_iterable(
            iters.combinations(np.arange(n_items, dtype=np.int64), r=r)
        ),
        dtype=np.int64,
        count=num_combinations * r,
    ).reshape(-1, r)
    return array[combination_idx]


def rSubset_iter(arr, r):
    # return list of all subsets of length r
    # to deal with duplicate subsets use
    # set(list(combinations(arr, r)))
    print()
    return list(iters.combinations(arr, r))


def DirectProduct_Ordered(
    List1: list, List2: list
):  # Direct product between two sets, intended for nodes
    list_of_pairs = list(iters.product(List1, List2))
    list_of_pairs = list(dict.fromkeys(list_of_pairs))
    list_of_pairs = [sorted(list(pairs)) for pairs in list_of_pairs]
    output = sorted(list_of_pairs, key=itemgetter(0))
    return output


def DirectProduct_Unordered(List1: list, List2: list):
    list_of_pairs = list(iters.product(List1, List2))
    return list_of_pairs


def RetrieveListRelationTypes(threshold):
    query = gds.run_cypher(
        """
            MATCH (source)-[r]->()
            WHERE NOT source:Patient
            WITH type(r) AS Rel_type, count(*) AS count
            WHERE count > $number // Adjust the threshold as needed
            RETURN Rel_type, count
            ORDER BY count DESC;
            """,
        {"number": threshold},
    )
    Rel_type = query["Rel_type"]
    Rel_type = Rel_type.tolist()
    return Rel_type


def SubgraphConnectivityByType_Descend():
    # Compute all combination of relationship type that yield a connected graph, output is stored in a CSV "Connected_Subgraph_ByType.csv"
    # In a descending order, i.e. it will start from the maximum combinations and work its way to only unique relation types
    RelationshipType = RetrieveListRelationTypes(10000)
    Not_Connected = (
        []
    )  # we will store all relation type who together yield disconnected graphs
    K = len(RelationshipType)
    for i in range(K + 1, 0, -1):
        Connected = []
        Tuples_of_Type = rSubset(RelationshipType, i)
        for Selection_Types in Tuples_of_Type:
            Selection_Types = list(Selection_Types)
            print(f"studying now {Selection_Types}")
            if all(
                Selection_Types not in Tuple for Tuple in Not_Connected
            ):  # check if we have already determined that these types yield a disconnected graph
                Nb_comp = gds.run_cypher(
                    f"""
                CALL gds.wcc.stats('noPatientUndirected',{{relationshipTypes:{Selection_Types}}})
                YIELD componentCount
                RETURN componentCount
                """
                )
                Nb_comp = Nb_comp["componentCount"][0]
                print(f"Number of connected component of Combination is{Nb_comp}")
                if Nb_comp > 1:
                    print(
                        "\n--------------------------------This Combination is not connected--------------------------------\n"
                    )
                    Not_Connected.append(Selection_Types)
                else:
                    print(
                        "\n----------------------------------This Combination is connected----------------------------------\n"
                    )
                    Connected.append(Selection_Types)
            else:
                print("This case has already been treated.")
        with open("Connected_Subgraph_ByType", "w", encoding="UTF-8") as f:
            # using csv.writer method from CSV package
            write = csv.writer(f)
            write.writerow(Connected)


def Exists_subset(Set_of_Subset, Testing_Set):
    for Subset in Set_of_Subset:
        if all(element in Testing_Set for element in Subset):
            return True
    return False


Threshold = 20000
defaultUpTo = len(RetrieveListRelationTypes(Threshold))


def SubgraphConnectivityByType_Ascend(
    From, UpTo=defaultUpTo
):  # Compute all combination of relationship type that yield a connected graph, output is stored in a CSV "Connected_Subgraph_ByType.csv"
    # In an ascending order, i.e. it will start from the unique relation types and work its way up to full combinations of types
    # From determines the size of tuples we begin with (e.g. From = 2 we will consider list of size 2) and UpTo determines up to which size we run our algorithm.
    RelationshipType = RetrieveListRelationTypes(Threshold)
    # Not_Connected = []#we will store all relation type who together yield disconnected graphs
    Connected = [
        ["ISA", "HAS_ROLE_GROUP"]
    ]  # I have already determined this yielded a connected graph
    for i in range(From, UpTo):
        Tuples_of_Type = rSubset(RelationshipType, i)
        for Selection_Types in Tuples_of_Type:
            Selection_Types = list(Selection_Types)
            print(f"studying now {Selection_Types}")
            if Exists_subset(
                Connected, Selection_Types
            ):  # check if we have already determined that these types yield a disconnected graph
                # Check if there exists in Connected a tuple which is a subset of the selected list of types.
                print(
                    "\n---------------------This Combination has already been treated and is connected.---------------------\n"
                )
            else:
                Nb_comp = gds.run_cypher(
                    """
                    CALL gds.wcc.stats('noPatientUndirected',{relationshipTypes:%s})
                    YIELD componentCount
                    RETURN componentCount
                """
                    % Selection_Types
                )
                Nb_comp = Nb_comp["componentCount"][0]
                print(f"Number of connected component of this Combination is {Nb_comp}")
                if Nb_comp > 1:
                    print(
                        "\n--------------------------------This Combination is not connected--------------------------------\n"
                    )
                    ## Not_Connected.append(Selection_Types) #we do not need in this case to remember disconnected subgraphs
                else:
                    print(
                        "\n----------------------------------This Combination is connected----------------------------------\n"
                    )
                    Connected.append(Selection_Types)

                    with open("Connected_Subgraph_ByType", "w", encoding="UTF-8") as f:
                        # using csv.writer method from CSV package
                        write = csv.writer(f)
                        write.writerow(Selection_Types)


def stats(
    ListOfPatients,
):  # computes the number of nodes max connected to a patient cluster and the number of nodes in average
    nb_of_patients = len(ListOfPatients)
    array_temp = np.zeros(nb_of_patients)
    counting = 0
    for Patient_id in ListOfPatients:
        nodes_id, _ = PrepareNodes(Patient_id)
        nb_of_nodes = len(nodes_id)
        array_temp[counting] = nb_of_nodes
        counting = counting + 1
    max_nb_nodes = np.max(array_temp)
    max_nb_nodes = max_nb_nodes.tolist()  # to obtain a float type
    max_nb_nodes = round(max_nb_nodes, 1)
    average = np.average(array_temp)
    average = average.tolist()
    average = round(average, 1)
    print(
        "The maximum number of nodes is "
        + repr(max_nb_nodes)
        + " and the average number of nodes per patient is "
        + repr(average)
    )


def ListOfNodesFromOnePatient(Patient_id: int):
    # given a patient id, return its nodes along with the patient id.
    gds = GraphDataScience(DB_URL, auth=AUTH)
    query1 = gds.run_cypher(
        "MATCH (n:Patient)--(m) WHERE id(n)=$id RETURN id(m)", {"id": Patient_id}
    )
    id_nodes = query1["id(m)"]
    output = pd.DataFrame({"List_of_nodes": id_nodes})
    output["appears_in_patient_id"] = output.apply(lambda x: [Patient_id], axis=1)
    output = output.values.tolist()
    return output


def ListOfNodesFromAllPatients(List_of_patients: list[int]):
    # From a list of patients id returns all unique nodes. Use multiprocessing for a faster treatment.
    All_nodes_messy = p_map(ListOfNodesFromOnePatient, List_of_patients)
    # return All_nodes_messy
    All_nodes_flat = [
        item for sublist in All_nodes_messy for item in sublist
    ]  # we unwind the resulting array (it is nested)
    All_nodes = pd.DataFrame(
        data=All_nodes_flat, columns=["nodes_id", "appears_in_patient_id"]
    )
    Unique_nodes = All_nodes.groupby("nodes_id").sum()  # we want distinct nodes
    Unique_nodes = Unique_nodes.reset_index(
        names=["nodes_id"]
    )  # we reset the index (groupby has set the index to "pairs_of_nodes")
    return Unique_nodes


def ListOfReference_Patients_Nodes_pid_single(Patient_id: int):
    # given a patient id, return its nodes along with the patient id.
    gds = GraphDataScience(DB_URL, auth=AUTH)
    query1 = gds.run_cypher(
        "MATCH (n:Patient)--(m) WHERE id(n)=$id RETURN id(m)", {"id": Patient_id}
    )
    query2 = gds.run_cypher(
        "MATCH (n:Patient) WHERE id(n)=$id RETURN n.pid AS pid", {"id": Patient_id}
    )
    id_nodes = query1["id(m)"].values.tolist()
    patient_pid = query2["pid"][0]
    output = [Patient_id, patient_pid, id_nodes]
    return output


def ListOfReference_Patients_Nodes_pid(List_of_patients: list[int]):
    # From a list of patients id returns all unique nodes. Use multiprocessing for a faster treatment.
    All_nodes = p_map(ListOfReference_Patients_Nodes_pid_single, List_of_patients)

    All_nodes = pd.DataFrame(
        data=All_nodes, columns=["patient_db_id", "pid", "list_of_nodes"]
    )
    return All_nodes


def Put_ShortestPath_table_into_good_form(Shortest_path_data):
    table_node_to_node = Shortest_path_data
    List_of_pairs_computed = table_node_to_node["pairs_of_nodes"].values.tolist()
    List_of_pairs_computed = [
        list(ast.literal_eval(pairs)) for pairs in List_of_pairs_computed
    ]
    List_of_patients = FindPatient_ids()
    List_of_nodes = ListOfNodesFromAllPatients(List_of_patients)
    List_of_nodes = List_of_nodes["nodes_id"].values.tolist()

    table_node_to_node["pairs_of_nodes"] = List_of_pairs_computed
    table_node_to_node = table_node_to_node[["pairs_of_nodes", "distance"]]
    # for an unknown reason, the .apply(list(reversed)) didn't work. Had to use a workaroung a bit ugly but functional
    tables_tobe_reversed = table_node_to_node.loc[:, "pairs_of_nodes"].values.tolist()
    distances = table_node_to_node.loc[:, "distance"].values.tolist()
    first_node_reversed = [pair_of_nodes[1] for pair_of_nodes in tables_tobe_reversed]
    second_node_reversed = [pair_of_nodes[0] for pair_of_nodes in tables_tobe_reversed]
    table_2 = pd.DataFrame(
        data={
            "first_node": first_node_reversed,
            "second_node": second_node_reversed,
            "distance": distances,
        }
    )
    # we add all identical pairs with distance of course 0
    table_3 = [[[node, node], 0] for node in List_of_nodes]
    table_3 = pd.DataFrame(data=table_3, columns=["pairs_of_nodes", "distance"])
    table_node_to_node = pd.concat([table_node_to_node, table_3], ignore_index=True)
    table_node_to_node[["first_node", "second_node"]] = pd.DataFrame(
        table_node_to_node["pairs_of_nodes"].tolist()
    )
    table_node_to_node = table_node_to_node.drop(columns=["pairs_of_nodes"])
    table_node_to_node = table_node_to_node.iloc[:, [1, 2, 0]]
    table_node_to_node = pd.concat([table_node_to_node, table_2])
    table_node_to_node.set_index(["first_node"], inplace=True)
    return table_node_to_node


def DeleteNodes(Nodes_to_delete):  # In entry an array of nodes_id
    records = gds.run_cypher(
        """
            MATCH (n)
            WHERE n.sctid IN $id
            DETACH DELETE n
            RETURN COUNT(n) AS AMOUNT
            """,
        {"id": Nodes_to_delete},
    )
    return records


def DeleteAllCounter():  # Delete all relation counters
    gds.run_cypher(
        """
        MATCH (n)-[r]-(m)
        WHERE r.counter IS NOT NULL
        REMOVE r.counter
        REMOVE r.patient
        RETURN r
        """
    )


def RetrieveStatsRelationTypes():  # stats about the types of relations
    global_result = pd.DataFrame(columns=["Relation_Type", "example"])
    query = gds.run_cypher(
        """
        MATCH (n)-[r]-(m)
        WHERE NOT n:Patient
        RETURN DISTINCT type(r) AS Relation_Type
        """
    )
    List_of_Types = query["Relation_Type"].values.tolist()
    print(List_of_Types)
    for specific_type in List_of_Types:
        result = gds.run_cypher(
            f"""
        MATCH (n)-[r]-(m)
        WHERE NOT n:Patient AND r:{specific_type} AND NOT (n.FSN = "None")
        RETURN DISTINCT type(r) AS Relation_Type, [n.FSN, m.FSN, m.sctid] AS example LIMIT 5
        """
        )
        print(result)
        global_result = pd.concat([global_result, result])
    global_result.to_csv("output.keppme/output.csv", index=False)


def PrepareNodes(
    Patient_id,
):  # Given a patient id, (returns one node per cluster) and all relations of the clusters
    query1 = gds.run_cypher(
        "MATCH (n:Patient)--(m) WHERE id(n)=$id RETURN id(m)", {"id": Patient_id}
    )
    id_nodes = query1["id(m)"]
    query2 = gds.run_cypher(
        "MATCH (n)-[r]-(m) WHERE id(n) IN $id AND id(m) IN $id RETURN DISTINCT id(r)",
        {"id": id_nodes},
    )
    id_relations = query2["id(r)"]
    return id_nodes, id_relations


def ProjectGraph():
    gds.run_cypher(
        """
        MATCH (source)-[r]->(target)
        WHERE NOT source:Patient
        WITH gds.graph.project(
            'noPatientUndirected',
            source,
            target,
            {},
            {undirectedRelationshipTypes: ['*']}) AS graph
        RETURN graph.nodeCount AS nodeCount,
            graph.relationshipCount AS relationshipCount
    """
    )
    result = gds.run_cypher("CALL gds.graph.list()")
    print(result)


def ProjectGraph_WithProperties():
    gds.run_cypher(
        """
        MATCH (source)-[r]->(target)
        WHERE NOT source:Patient AND NOT source:QualifierValue AND NOT target:QualifierValue
        WITH gds.graph.project(
            'noPatientUndirected_WithProperties',
            source,
            target,
            {relationshipProperties: {cost_AverageInverseIC: r.cost_AverageInverseIC,cost: r.cost}},
            {undirectedRelationshipTypes: ['*']}) AS graph
        RETURN graph.nodeCount AS nodeCount,
            graph.relationshipCount AS relationshipCount
        """
    )


def ProjectGraphWithPatients():
    gds.run_cypher(
        """
        MATCH (source)-[r]->(target)
        WITH gds.graph.project(
            'withPatientUndirected',
            source,
            target,
            {},
            {undirectedRelationshipTypes: ['*']}) AS graph
        RETURN graph.nodeCount AS nodeCount,
            graph.relationshipCount AS relationshipCount
        """
    )
    result = gds.run_cypher("CALL gds.graph.list()")
    print(result)


def ProjectGraphWithTypes():  # Project the graph with the different types of relaiton. Use a lot of RAM (~13gib)
    query = gds.run_cypher(
        """
        MATCH (source)-[r]->(target)
        WHERE NOT source:Patient
        RETURN DISTINCT type(r) As Rel_type
        """
    )
    Rel_type = query["Rel_type"]
    Rel_type = Rel_type.tolist()
    query = gds.run_cypher(
        """
        MATCH (source)-[r]->(target)
        WHERE NOT source:Patient
        RETURN DISTINCT labels(source) As Node_labels
        """
    )
    Label_nodes = []
    List = query["Node_labels"]
    List = List.tolist()
    for category in List:
        for label in category:
            Label_nodes.append(label)
    Label_nodes = list(dict.fromkeys(Label_nodes))
    query = gds.run_cypher(
        f"""
        CALL gds.graph.project('Snomed',{Label_nodes},{Rel_type})
        YIELD
            graphName AS graph,
            relationshipProjection AS knowsProjection,
            nodeCount AS nodes,
            relationshipCount AS rels
        RETURN nodes
        """
    )
    result = gds.run_cypher("CALL gds.graph.list()")
    print(result)


def ProjectGraphWithTypes_WithProperties():  # Project the graph with the different types of relaiton. Use a lot of RAM (~13gib)
    query = gds.run_cypher(
        """
            MATCH (source)-[r]->(target)
            RETURN DISTINCT type(r) As Rel_type
            """
    )
    Rel_list = query["Rel_type"]
    Rel_list = Rel_list.tolist()
    Rel_type = "{"
    Last_item = Rel_list[len(Rel_list) - 1]
    for specific_type in Rel_list:
        if Rel_list != Last_item:
            Rel_type = Rel_type + specific_type + ": {orientation: 'UNDIRECTED'}, "
    Rel_type = Rel_type + Last_item + ": {orientation: 'UNDIRECTED'}}"

    query = gds.run_cypher(
        """
            MATCH (source)-[r]->(target)
            RETURN DISTINCT labels(source) As Node_labels
            """
    )
    Label_nodes = []
    List = query["Node_labels"]
    List = List.tolist()
    for category in List:
        for label in category:
            Label_nodes.append(label)
    Label_nodes = list(dict.fromkeys(Label_nodes))
    gds.run_cypher(
        f"""
        CALL gds.graph.project('NoPatientUndirected',{Label_nodes},{Rel_type})
        YIELD
            graphName AS graph,
            relationshipProjection AS knowsProjection,
            nodeCount AS nodes,
            relationshipCount AS rels
        RETURN nodes
        """
    )
    result = gds.run_cypher("CALL gds.graph.list()")
    print(result)


def ProjectGraphWithTypes_NoPatients_Undirected_WithCost():  # Project the graph with the different types of relaiton. Use a lot of RAM (~13gib)
    query = gds.run_cypher(
        """
            MATCH (source)-[r]->(target)
            WHERE NOT source:Patient
            RETURN DISTINCT type(r) As Rel_type
            """
    )
    Rel_list = query["Rel_type"]
    Rel_list = Rel_list.tolist()
    Rel_type = "{"
    Last_item = Rel_list[len(Rel_list) - 1]
    for specific_type in Rel_list:
        if Rel_list != Last_item:
            Rel_type = (
                Rel_type
                + specific_type
                + ": {orientation: 'UNDIRECTED', properties: 'cost_AverageInverseIC'}, "
            )
    Rel_type = (
        Rel_type
        + Last_item
        + ": {orientation: 'UNDIRECTED', properties: 'cost_AverageInverseIC'}}"
    )

    query = gds.run_cypher(
        """
            MATCH (source)-[r]->(target)
            WHERE NOT source:Patient
            RETURN DISTINCT labels(source) As Node_labels
            """
    )
    Label_nodes = []
    List = query["Node_labels"]
    List = List.tolist()
    for category in List:
        for label in category:
            Label_nodes.append(label)
    Label_nodes = list(dict.fromkeys(Label_nodes))
    gds.run_cypher(
        f"""
        CALL gds.graph.project('noPatientUndirected_WithProperties',{Label_nodes},{Rel_type})
        YIELD
            graphName AS graph,
            relationshipProjection AS knowsProjection,
            nodeCount AS nodes,
            relationshipCount AS rels
        RETURN nodes
        """
    )
    result = gds.run_cypher("CALL gds.graph.list()")
    print(result)


def ProjectGraphJustISA():  # Project the graph with the different types of relaiton. Use a lot of RAM (~13gib)
    query = gds.run_cypher(
        """
        MATCH (source)-[r]->(target)
        WHERE r:ISA
        RETURN DISTINCT type(r) As Rel_type
        """
    )
    Rel_type = query["Rel_type"]
    Rel_type = Rel_type.tolist()
    query = gds.run_cypher(
        """
        MATCH (source)-[r]->(target)
        WHERE r:ISA
        RETURN DISTINCT labels(source) As Node_labels
        """
    )
    Label_nodes = []
    List = query["Node_labels"]
    List = List.tolist()
    for category in List:
        for label in category:
            Label_nodes.append(label)
    Label_nodes = list(dict.fromkeys(Label_nodes))
    gds.run_cypher(
        f"""
        CALL gds.graph.project('Snomed',{Label_nodes},{Rel_type})
        YIELD
            graphName AS graph,
            relationshipProjection AS knowsProjection,
            nodeCount AS nodes,
            relationshipCount AS rels
        RETURN nodes
        """
    )
    result = gds.run_cypher("CALL gds.graph.list()")
    print(result)


def AddLabel_Totreat(List_of_nodes, label):
    gds.run_cypher(
        f"""
        MATCH (n)
        WHERE id(n) IN $List_id
        SET n:{label}
        RETURN labels(n)
        """,
        {"List_id": List_of_nodes},
    )


def RemoveLabel_Totreat(List_of_nodes, label):
    gds.run_cypher(
        f"""
        MATCH (n)
        WHERE id(n) IN $List_id
        REMOVE n:{label}
        RETURN labels(n)
        """,
        {"List_id": List_of_nodes},
    )


def CheckLabel_Totreat(List_of_nodes, label):
    query = gds.run_cypher(
        f"""
        MATCH (n)
        WHERE id(n) IN $List_id
        WITH exists((n:{label})-[]-()) as existence
        RETURN DISTINCT existence
        """,
        {"List_id": List_of_nodes},
    )
    if query["existence"][0] == "true":
        return True
    else:
        return False


def ShortestPath_Between2Nodes_Filtering(
    pair_of_nodes, Snomed_Concept_Labels
):  # Return a path of relations Ids given a projected graph, a source id and a target id. We supposed that the filtering Id was in the second entry of the pair of nodes.
    G = gds.graph.get("Snomed")  # we specify here the projected graph
    Patient_id_selected = pair_of_nodes[1]
    labels = ["Patient_" + str(Patient_id_selected)] + Snomed_Concept_Labels
    result = gds.shortestPath.dijkstra.stream(
        G, sourceNode=pair_of_nodes[0], targetNode=pair_of_nodes[1], nodeLabels=labels
    )
    nodeIds = result["nodeIds"][
        0
    ]  # we select in the panda dataframe the entry corresponding to the path of nodes
    # query = gds.run_cypher(
    #             """
    #             MATCH (n)-[r]->(m)
    #             WHERE id(n) IN $List_id AND id(m) IN $List_id
    #             RETURN DISTINCT id(r)
    #             """,
    #             {"List_id":nodeIds}
    #             ) #we want a list of relations id, not of nodes.
    # id_relations = query["id(r)"
    # id_relations = id_relations.values.tolist()
    return nodeIds


def ShortestPath_All_FromLabel(
    label,
):  # Use the gds function all paths to compute all paths between nodes with a given label
    output = gds.run_cypher(
        f"""
        CALL gds.alpha.allShortestPaths.stream(
            'Snomed',
            {{nodeLabels: ['{label}']}})
        YIELD sourceNodeId,targetNodeId,distance
        RETURN sourceNodeId,targetNodeId,distance
        """
    )
    return output


def Set_Cost_Qualifier_Value_single(id_node: int):
    gds = GraphDataScience(DB_URL, AUTH)
    gds.run_cypher(f"MATCH (n) WHERE id(n)={id_node} SET n.cost=1 RETURN n")

def check_if_IC_property_exists():
    gds = GraphDataScience(DB_URL, AUTH)
    query = gds.run_cypher("MATCH (n:ObjectConcept) WHERE NOT n.cost_InverseIC IS NOT NULL RETURN n")
    query = query.values.tolist()
    check = not query
    return check

def Set_Cost_Qualifier_Value_all():
    id_nodes = gds.run_cypher(
        "MATCH (n) WHERE NOT n:QualifierValue RETURN DISTINCT id(n) AS id"
    )
    id_nodes = id_nodes.values.tolist()
    id_nodes = [id_single for id_bloc in id_nodes for id_single in id_bloc]
    p_umap(Set_Cost_Qualifier_Value_single, id_nodes)


def export_csv_metrics():
    list_metric = [
        "Metric_BagsOfFindings",
        "Metric_AverageLinks",
        "Metric_AverageLinks_ICWeights",
        "Metric_ICPath",
    ]
    Patient_pid = pd.read_csv(
        "output.keepme/Output_data_ListOfReference_Patients_Nodes_pid.csv"
    )
    Patient_pid = Patient_pid["pid"].values.tolist()
    matched_patients = pd.read_csv("input.keepme/matched_patients.csv", sep=";")
    matched_patients_id = (
        matched_patients["cohort"].str[:2] + "-" + matched_patients["old_id"]
    )
    matched_patients_id = matched_patients_id.values.tolist()
    for metric_name in tqdm(list_metric, "compression of the ultimate metrics"):
        metric = pd.read_csv(
            "output.keepme/" + metric_name + "_RAW.csv", index_col="Unnamed: 0"
        )
        # metric.columns = Patient_pid
        # metric.index = Patient_pid
        metric.index.name = "pid"
        print(
            f"---------------------computing the average metric of {metric_name}---------------------"
        )
        metric = (metric.T + metric) * 0.5
        print(metric)
        # print(
        #     f"---------------------compression of the metric {metric_name}---------------------"
        # )
        # metric.to_csv(
        #     "output.keepme/" + metric_name + ".bz2",
        #     index_label="pid",
        #     compression="bz2",
        # )
        # metric.to_csv("output.keepme/" + metric_name + ".csv", index_label="pid")
        # print("-------------Filtering the metric from matched selection-------------")
        # metric_matched = metric.loc[matched_patients_id, matched_patients_id]
        # metric.index.name = "pid"
        # print(
        #     f"---------------------compression of the matched metric {metric_name}---------------------"
        # )
        # metric_matched.to_csv(
        #     "output.keepme/" + metric_name + "_MatchedSelection.csv", index_label="pid"
        # )
