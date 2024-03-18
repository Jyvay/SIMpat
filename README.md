# Graph Project - JV Voegeli

Here is the code done by Jean-Virgile Voegeli during his interniship at Simed.

Link for the default Neo4j db used HERE : https://drive.google.com/file/d/1NUY1UT8bxvSfsrANeudnV2CN0t8FPFEr/view?usp=sharing

## Shortest path between concepts of patient synthetic files

For the first part, the task was to create a Neo4j database including Snomed CT browser and a generated synthetic patient dataset using Synthea. Then, the goal was to build functions to retrieve all distinct pair of nodes included in the synthetic patient files and to find the shortest path between them through Snomed CT, using Neo4j. At the same time, each relation used in the shortest path process had a counter added and updated each time that relation was part of a shortest path. For more information concerning this aspect of the project, please refer to section 3 of the Graph Project main report https://www.overleaf.com/read/hxgrctfytscq#251ad4.

The code for this part is available in the file ShortestPathPatient.py and a csv output can be obtained using the function Output_data_ShortestPath_patients() from Output_data.py. This will create a csv file called 'Output_data_ShortestPath_patients.csv' in 'output.keepme' repository.

The function ShortestPath_BetweenPatients_And_Counter_Update() in the file ShortestPathPatient.py is the main algorithm of the first part of this project. Running it will retrieve all patient files in the Neo4j database, prepare distinct pairs of nodes, run a shortest path algorithm, and finally will create and update according counters.

## Metrics for comparing patient synthetic EHR graphs

The second part of the project was to build metrics to compare two patient dataset, imported as graphs in the Neo4j database. The metrics are based on the paper 'Inter-patient distance metrics using snomed ct defining relationships.' by Melton and al. (2006). More informations about metrics can be found on section 4 of the Graph Project main report https://www.overleaf.com/read/hxgrctfytscq#251ad4.

To obtain the Metrics, one needs to run in a sequential order the functions<br />

Output_data_ShortestPath_Node_to_Node(),<br />
Output_data_ShortestPath_Node_to_Node_ICWeights(),<br /> <-- Please read remark below
Output_data_AllPatientsWithFourMetrics()<br />

To make the generation of metric similarity simpler, an assistant has been created. to launch it, just run the function Metric_Similarity_Assistant() in the file 'Output_data.py'. The assistant will check if needed input are available before running a particular step mentionned above.

### detailled explanations of the functions used

Output_data_ShortestPath_Node_to_Node() will retrieve all distinct pair of nodes, use an already computed file Output_data_ShortestPath_Node_to_Node.csv if existing to spare some running time, and compute the remaining shortest path not alreay computed. The difference with the pairs of nodes used in the first part of the project is the following: in the first part, we paired concomittent nodes of synthetic EHR whereas here we paired all nodes that appears in distinct synthetic EHR. For instance, if we have a patient 1 with nodes [1,2,3] and a patient 2 with nodes [1,5,6], in the first part we ran shortest path on [[1,2],[1,3],[2,3],[1,5],[1,6],[5,6]] whereas here we find the shortest paths for [[1,5],[1,6],[2,5],[2,6],[3,5],[3,6]] which is a different set of pairs of nodes. All shortest paths will be stored in a csv file called 'Output_data_ShortestPath_Node_to_Node.csv' in 'output.keepme'

Output_data_ShortestPath_Node_to_Node_ICWeights() will do the same thing as Output_data_ShortestPath_Node_to_Node() but instead of finding the shortest path minimizing the number of edges contained in the path, it will find the shortest path that maximizes the sum of the information content of the nodes included in the path. Information content is a semantic notion representing how much information a node in an onthology contains. For a simple example, if we take the tree of lifes representing the biology standard classification, "horse" will contain more information than "mammal" (if I say 'I own a mammal' instead of saying 'I own a horse', a lot of information is lost). All shortest distance will be stored in a csv file called 'Output_data_ShortestPath_Node_to_Node_ICWeights.csv' in 'output.keepme'

!! IMPORTANT NOTE !!<br />
Output_data_ShortestPath_Node_to_Node_ICWeights() supposes that all Snomed nodes in the database has a property called 'cost_InverseIC' which can be retrieved using from Ouput_data.py the function Output_data_Nodes_DescendantsNb() and then running from Metrics.py InformationContentAll_WRITECOST(), and InformationContentRelationships_WRITECOST. Note that if RouleGroups are present in the Db (which they are if you are using the by default db) you need to do a supplementary step:

Rolegroups node are not included initially when executing the above code. If Rolegroups nodes are present in the neo4j db, you also need to run the function InformationContentRoleGroup_WRITECOST().

Note that if you are using the standard db coming with this project, all nodes (including RoleGroups) have already the property 'cost_AverageInverseIC' computed. This needs to be done ONLY when starting fresh with a new importation of Snomed CT and synthetic patient files.

Output_data_AllPatientsWithFourMetrics() computes all metrics from the files Output_data_ShortestPath_Node_to_Node_.csv and Output_data_ShortestPath_Node_to_Node_ICWeights, supposedly found in output.keepme file, and return the four metrics in NxN matrices csv files, where N is the number of patients.
