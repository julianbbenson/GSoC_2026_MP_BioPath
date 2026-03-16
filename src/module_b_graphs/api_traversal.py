import requests
import networkx as nx
import json
from typing import List, Dict, Any, Set

class ReactomeMetabolicEngine:
    def __init__(self):
        self.base_url = "https://reactome.org/ContentService/data"
        self.headers = {'accept': 'application/json'}
        self.G = nx.MultiDiGraph() # MultiDiGraph allows multiple edge types (containment vs. flow)
        self.visited_pathways: Set[str] = set()

    def fetch_participants(self, reaction_id: str):
        """Retrieves physical inputs/outputs for a specific reaction."""
        endpoint = f"{self.base_url}/participants/{reaction_id}/physicalEntities"
        try:
            response = requests.get(endpoint, headers=self.headers)
            if response.status_code == 200:
                return response.json()
        except:
            return []
        return []

    def recursive_stitch(self, pathway_id: str, parent_node: str = None):
        """Recursive DFS that pulls hierarchy AND metabolic flow."""
        if pathway_id in self.visited_pathways:
            return
        self.visited_pathways.add(pathway_id)

        endpoint = f"{self.base_url}/pathway/{pathway_id}/containedEvents"
        print(f"[*] Processing: {pathway_id}...")
        
        try:
            events = requests.get(endpoint, headers=self.headers).json()
            for event in events:
                if not isinstance(event, dict): continue
                stId = event.get('stId')
                schema_class = event.get('schemaClass')
                
                self.G.add_node(stId, name=event.get('displayName'), type=schema_class)
                
                # Edge Type 1: Containment (Hierarchy)
                if parent_node:
                    self.G.add_edge(parent_node, stId, relationship="contains")

                if schema_class == 'Pathway':
                    self.recursive_stitch(stId, parent_node=stId)
                elif schema_class in ['Reaction', 'BlackBoxEvent']:
                    # Edge Type 2: Metabolic Flow
                    participants = self.fetch_participants(stId)
                    for p in participants:
                        p_id = p.get('peDbId') # Physical Entity ID
                        role = p.get('schemaClass') # Input/Output/Catalyst
                        # Logic to link reactions via shared participants would happen here
                        # For now, we tag nodes with their participants for the solver
                        nx.set_node_attributes(self.G, {stId: {'participants': participants}})

        except Exception as e:
            print(f"[!] Error: {e}")

    def export_solver_json(self, filename: str):
        data = nx.node_link_data(self.G)
        with open(filename, 'w') as f:
            json.dump(data, f, indent=4)
        print(f"\n[✔] Enriched JSON exported to: {filename}")

if __name__ == "__main__":
    engine = ReactomeMetabolicEngine()
    engine.recursive_stitch("R-HSA-5673001")
    engine.export_solver_json("data/enriched_pathway.json")