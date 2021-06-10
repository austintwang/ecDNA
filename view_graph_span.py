import pickle
import os
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import seaborn as sns

CHROM_MAP_PATH = "GRCh38_RefSeq2UCSC.txt"
CHROM_MAP = {}
with open(CHROM_MAP_PATH) as f:
    for l in f:
        seq, chrom = l.strip().split()
        CHROM_MAP[seq] = chrom

def plot_component(graph, comp, result_path, result_path_nl):
    G = graph.subgraph(comp)
    pos = nx.kamada_kawai_layout(G)
    # node_len_dict = nx.get_node_attributes(G, "len")
    # node_sizes = [node_len_dict[i] / 1e4 for i in G.nodes()]
    edge_weight_dict = nx.get_edge_attributes(G, "min_cn")
    edge_weights = [edge_weight_dict[i] for i in G.edges()]
    # node_labels = nx.get_node_attributes(G, "range")

    cmap = plt.cm.plasma

    nodes = nx.draw_networkx_nodes(
        G, 
        pos, 
        # node_size=node_sizes, 
        node_color="indigo"
    )
    edges = nx.draw_networkx_edges(
        G,
        pos,
        # node_size=node_sizes,
        arrowstyle="->",
        arrowsize=10,
        edge_color=edge_weights,
        edge_cmap=cmap,
        width=2,
        edge_vmin=0,
        edge_vmax=5
    )
    plt.savefig(result_path_nl, bbox_inches='tight')

    # labels = nx.draw_networkx_labels(
    #     G, 
    #     pos, 
    #     labels=node_labels, 
    #     font_size=8
    # )
    labels = nx.draw_networkx_edge_labels(
        G, 
        pos, 
        labels=edge_weights, 
        font_size=8
    )
    plt.savefig(result_path, bbox_inches='tight')
    plt.clf()
    plt.close()

def view_graph_span(in_path, results_dir):
    with open(in_path, "rb") as in_file:
        node_data, edge_data = pickle.load(in_file)

    graph = nx.DiGraph()
    graph.add_nodes_from(node_data)
    graph.add_edges_from(edge_data)
    # print(graph) ####

    components = nx.weakly_connected_components(graph)
    freqs = {}
    
    plot_dir = os.path.join(results_dir, "graph_span")
    os.makedirs(plot_dir, exist_ok=True)
    for ind, c in enumerate(components):
        freqs.setdefault(len(c), 0) ####
        freqs[len(c)] += 1 ####
        if len(c) <= 2:
            continue
        print(len(c)) ####
        result_path = os.path.join(plot_dir, f"span_{ind:02d}.svg")
        result_path_nl = os.path.join(plot_dir, f"span_{ind:02d}_nl.svg")
        # plot_component(graph, c, result_path, result_path_nl)

    for k in sorted(freqs.keys()):
        print(k, freqs[k]) ####

if __name__ == "__main__":
    data_dir = "/oak/stanford/groups/wjg/atwang/ecdna/data"

    results_dir_base = "/oak/stanford/groups/wjg/atwang/ecdna/results"

    in_path = os.path.join(data_dir, "COLO320DM_gDNA_nanopore_guppy_4.4_span.pickle")
    results_dir = os.path.join(results_dir_base, "COLO320DM_gDNA_nanopore_guppy_4.4")
    os.makedirs(results_dir, exist_ok=True)
    view_graph_span(in_path, results_dir)

    in_path = os.path.join(data_dir, "PC3_gDNA_combined_span.pickle")
    results_dir = os.path.join(results_dir_base, "PC3_gDNA_combined")
    os.makedirs(results_dir, exist_ok=True)
    view_graph_span(in_path, results_dir)