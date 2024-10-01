import re
import networkx as nx
from tqdm import tqdm


# Extract GO terms and their associated branches
def extract_go_terms_and_branches(file_path):
    with open(file_path, 'r') as file:
        content = file.read()
        stanzas = re.findall(r'\[Term\][\s\S]*?(?=\n\[|$)', content)

    branch_abbr = {'biological_process': 'BPO', 'cellular_component': 'CCO', 'molecular_function': 'MFO'}
    go_terms_dict = {}

    for stanza in stanzas:
        go_id = re.search(r'^id: (GO:\d+)', stanza, re.MULTILINE)
        namespace = re.search(r'^namespace: (\w+)', stanza, re.MULTILINE)

        if go_id and namespace:
            go_terms_dict[go_id.group(1)] = branch_abbr[namespace.group(1)]

    return go_terms_dict


# Class to manage protein predictions
class ProteinPredictions:
    def __init__(self):
        self.predictions = {}

    def add_prediction(self, protein, go_term, score, branch, bonus=0.1):
        if protein not in self.predictions:
            self.predictions[protein] = {'CCO': {}, 'MFO': {}, 'BPO': {}}

        score = min(float(score) + bonus, 1)
        if go_term in self.predictions[protein][branch]:
            self.predictions[protein][branch][go_term] = max(self.predictions[protein][branch][go_term], score)
        else:
            self.predictions[protein][branch][go_term] = score

    def get_predictions(self, output_file='predictions.tsv', top=30):
        with open(output_file, 'w') as f:
            for protein, branches in self.predictions.items():
                for branch, go_terms in branches.items():
                    top_go_terms = sorted(go_terms.items(), key=lambda x: x[1], reverse=True)[:top]
                    for go_term, score in top_go_terms:
                        f.write(f"{protein}\t{go_term}\t{score:.3f}\n")


# Read GO OBO file and build directed graph
def read_go_obo(file_path):
    graph = nx.DiGraph()
    with open(file_path, 'r') as file:
        current_term = None
        for line in file:
            line = line.strip()
            if line == "[Term]":
                current_term = {}
            elif not line and current_term:
                graph.add_node(current_term['id'], **current_term)
                current_term = None
            elif current_term:
                key, value = line.split(": ", 1)
                current_term[key] = value
                if key == "is_a":
                    parent_id = value.split(" ! ")[0]
                    graph.add_edge(parent_id, current_term['id'])
                elif key == "relationship" and "part_of" in value:
                    parent_id = value.split("part_of ")[1].split(" ! ")[0]
                    graph.add_edge(parent_id, current_term['id'])
    return graph


# Get all ancestor nodes for a given term
def get_ancestors(graph, start_term):
    ancestors = set()
    stack = [start_term]
    while stack:
        current_term = stack.pop()
        ancestors.add(current_term)
        stack.extend(graph.predecessors(current_term))
    return ancestors


# Add ancestor terms and their scores
def add_ancestors_to_predictions(graph, predictions):
    sample_dict = {}
    for sample_id, term, score in predictions:
        sample_dict.setdefault(sample_id, []).append((term, score))

    new_predictions = []
    for sample_id, terms_scores in tqdm(sample_dict.items(), desc="Processing samples"):
        for term, score in terms_scores:
            ancestors = get_ancestors(graph, term)
            for ancestor in ancestors:
                ancestor_score = min(score + (terms_scores.count((term, score)) - 1) * 0.05, 1.0)
                new_predictions.append((sample_id, ancestor, ancestor_score))

    return new_predictions


# File paths and GO term extraction
file_path = 'path_to_go-basic.obo'
go_terms_dict = extract_go_terms_and_branches(file_path)
protein_predictions = ProteinPredictions()

# Process protein predictions
for l in tqdm(open('path_to_quickgo.tsv')):
    temp_id, go = l.split('\t')[1:3]
    go = go.strip()
    if go in go_terms_dict:
        protein_predictions.add_prediction(temp_id, go, 1, go_terms_dict[go])

# Export prediction results
protein_predictions.get_predictions()

# Build graph and add ancestor terms
graph = read_go_obo(file_path)
new_predictions = add_ancestors_to_predictions(graph, protein_predictions.predictions)