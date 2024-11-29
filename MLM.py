# -*- coding: utf-8 -*-
"""
MLM Pipeline - Multi-Level-Migration Analysis
Author: Naclist
Version: 0.1
Date: 2024/11/29

This script performs Multi-Level-Migration (MLM) analysis using phylogenetic trees
and migration data from TreeTime to infer ancestral node statuses and weight nodes
in a phylogenetic tree according to their evolutionary significance.
"""

import sys
from io import StringIO
from Bio import Phylo
import pandas as pd
import logging

# Set up logging configuration for error tracking and info logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')


# Function to parse a Newick format tree string into a tree object
def parse_newick_tree(newick):
    """
    Parse a Newick formatted phylogenetic tree into a Bio.Phylo tree object.

    Args:
    newick (str): Newick string representing the phylogenetic tree.

    Returns:
    tree (Bio.Phylo.Tree): Phylogenetic tree object.
    """
    try:
        tree = Phylo.read(StringIO(newick), "newick")
        logging.info("Successfully parsed Newick tree.")
        return tree
    except Exception as e:
        logging.error(f"Error parsing Newick tree: {e}")
        raise


# Function to read a confidence matrix from a file
def read_confidence(confidence_mat):
    """
    Read the confidence matrix from a tab-separated file.

    Args:
    confidence_mat (str): Path to the confidence matrix file (CSV format).

    Returns:
    confi_df (pandas.DataFrame): DataFrame containing confidence values.
    """
    try:
        confi_df = pd.read_csv(confidence_mat, sep='\t')
        logging.info("Successfully read the confidence matrix.")
        return confi_df
    except Exception as e:
        logging.error(f"Error reading confidence matrix: {e}")
        raise


# Function to get the confidence value for a specific node and region
def get_confidence_by_name_and_region(confi_df, name, region):
    """
    Retrieve the confidence value for a given node and region from the confidence matrix.

    Args:
    confi_df (pandas.DataFrame): DataFrame containing confidence values.
    name (str): Name of the node to look up.
    region (str): Region to search for in the columns.

    Returns:
    float: The confidence value, defaulting to 1 if not found.
    """
    if region not in confi_df.columns:
        logging.warning(f"Region {region} not found in the confidence matrix. Returning default confidence value of 1.")
        return 1

    row = confi_df[confi_df[0] == name]

    if not row.empty and region in row.columns:
        return row.iloc[0][region]
    else:
        logging.warning(f"Node {name} not found for region {region}. Returning default confidence value of 1.")
        return 1


# Function to create a dictionary representing the nodes, their parents, and children in the tree
def create_node_dict(tree):
    """
    Create a dictionary representing the nodes in the phylogenetic tree.

    Args:
    tree (Bio.Phylo.Tree): The phylogenetic tree object to parse.

    Returns:
    node_dict (dict): Dictionary where keys are node names and values are dictionaries with 'parent', 'children', and 'weight'.
    """
    node_dict = {}
    for clade in tree.find_clades():
        node_dict[clade.name] = {
            'children': [],
            'parent': None,
            'weight': 0,
        }
    for clade in tree.find_clades():
        for child in clade.clades:
            node_dict[child.name]['parent'] = clade.name
            node_dict[clade.name]['children'].append(child.name)
    return node_dict


# Function to assign weights to the nodes based on the confidence matrix and tree structure
def calculate_weights(node_dict, confi_df):
    """
    Calculate the weights for each node based on its confidence value and tree structure.

    Args:
    node_dict (dict): Dictionary representing nodes in the phylogenetic tree.
    confi_df (pandas.DataFrame): DataFrame containing confidence values for each node.
    """

    def assign_weight(node):
        region = node_dict[node].get('region', 'Unknown')
        confidence = get_confidence_by_name_and_region(confi_df, node, region)
        if node_dict[node]['children']:
            weight = sum(assign_weight(child) * confidence for child in node_dict[node]['children'])
        else:
            weight = 1 * confidence
        weights[node] = weight
        return weight

    root = next(node for node, data in node_dict.items() if data['parent'] is None)
    weights = {node: 0 for node in node_dict}
    assign_weight(root)
    for node in node_dict:
        node_dict[node]['weight'] = weights[node]


# Function to assign geographical regions to each node
def assign_geography(node_dict, annotations):
    """
    Assign geographical region annotations to nodes in the tree.

    Args:
    node_dict (dict): Dictionary representing nodes in the phylogenetic tree.
    annotations (dict): Dictionary where keys are node names and values are regions.
    """
    for node, region in annotations.items():
        if node in node_dict:
            node_dict[node]['region'] = region
        else:
            logging.warning(f"Node {node} is not in the node dictionary.")


# Function to calculate the spread probabilities for each region
def calculate_spread_probabilities(node_dict, annotations):
    """
    Calculate spread probabilities for each region based on node weights and parent-child relationships.

    Args:
    node_dict (dict): Dictionary representing nodes in the phylogenetic tree.
    annotations (dict): Dictionary where keys are node names and values are regions.

    Returns:
    regions_weight (dict): Dictionary containing spread probabilities for each region.
    """
    regions_weight = {region: {'total': 0, 'sources': {}, 'local': 0} for region in set(annotations.values())}
    for node, data in node_dict.items():
        region = data.get('region', 'Unknown')
        regions_weight[region]['total'] += data['weight']
        parent_region = node_dict[data['parent']]['region'] if data['parent'] else None

        if parent_region == region:
            regions_weight[region]['local'] += data['weight']
        elif parent_region and parent_region != region:
            regions_weight[region]['sources'].setdefault(parent_region, 0)
            regions_weight[region]['sources'][parent_region] += data['weight']

    return regions_weight


# Main function to run the entire workflow
def main():
    if len(sys.argv) != 4:
        logging.error("Usage: python MLM.py [Tree] [Trait] [Confidence.mat]")
        sys.exit(1)

    newick_file = sys.argv[1]  # Path to the Newick tree file
    trait_file = sys.argv[2]  # Path to the Trait file (meta file generated from TreeTime .nex)
    confidence_file = sys.argv[3]  # Path to the Confidence matrix file

    # Read the Newick tree, trait file, and confidence matrix
    with open(newick_file, 'r') as f:
        newick = f.read()

    tree = parse_newick_tree(newick)
    confi_df = read_confidence(confidence_file)

    # Parse the trait annotations (assuming it is a dictionary of node -> region)
    annotations = {}
    with open(trait_file, 'r') as f:
        for line in f:
            node, region = line.strip().split('\t')  # assuming tab-separated values
            annotations[node] = region

    # Create the node dictionary for the phylogenetic tree
    node_dict = create_node_dict(tree)

    # Assign geographical regions to the nodes based on the annotations
    assign_geography(node_dict, annotations)

    # Calculate the weights for the nodes based on the tree structure and confidence matrix
    calculate_weights(node_dict, confi_df)

    # Calculate the spread probabilities for each region
    spread_probabilities = calculate_spread_probabilities(node_dict, annotations)

    # Print the results for spread probabilities
    logging.info("Spread probabilities by region:")
    for region, data in spread_probabilities.items():
        logging.info(
            f"Region: {region}, Total Weight: {data['total']}, Local Weight: {data['local']}, Sources: {data['sources']}")


if __name__ == "__main__":
    main()

