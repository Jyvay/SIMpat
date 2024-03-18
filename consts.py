from graphdatascience import GraphDataScience

USER = "neo4j"
PASSWORD = "password"
DB_URL = "bolt://db:7687"
AUTH = ("neo4j", "password")


gds = GraphDataScience(DB_URL, auth=(USER, PASSWORD))
