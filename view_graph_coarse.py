import pickle
import os
import numpy as np
import networkx as nx
import matplotlib as mpl
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
    node_len_dict = nx.get_node_attributes(G, "len")
    node_sizes = [node_len_dict[i] / 1e4 for i in G.nodes()]
    edge_weight_dict = nx.get_edge_attributes(G, "freq")
    edge_colors = [edge_weight_dict[i] for i in G.edges()]
    node_labels = nx.get_node_attributes(G, "range")

    cmap = sns.color_palette("flare", as_cmap=True)

    nodes = nx.draw_networkx_nodes(
        G, 
        pos, 
        node_size=node_sizes, 
        node_color="indigo"
    )
    edges = nx.draw_networkx_edges(
        G,
        pos,
        node_size=node_sizes,
        arrowstyle="->",
        arrowsize=10,
        edge_color=edge_colors,
        edge_cmap=cmap,
        width=2,
        edge_vmin=0,
        edge_vmax=5,
        alpha=0.6
    )
    pc = mpl.collections.PatchCollection(edges, cmap=cmap)
    pc.set_array(edge_colors)
    plt.colorbar(pc)

    ax = plt.gca()
    ax.set_axis_off()

    plt.savefig(result_path_nl, bbox_inches='tight')

    labels = nx.draw_networkx_labels(
        G, 
        pos, 
        labels=node_labels, 
        font_size=8
    )
    plt.savefig(result_path, bbox_inches='tight')
    plt.clf()
    plt.close()

def view_graph_coarse(in_path, results_dir):
    with open(in_path, "rb") as in_file:
        pt_data, node_data, edge_data, coarse_edges = pickle.load(in_file)

    blocks, blocks_cn, pt_block_map = node_data
    nodelist = []
    for k, v in blocks.items():
        start = v[0]
        end = v[-1]
        chrom = CHROM_MAP.get(start[0], start[0])
        strand = start[1]
        block_len = abs(start[2] - end[2])
        block_range = f"{chrom}{strand} {abs(int(start[2]))}:{abs(int(end[2]))}"
        data = (k, {"len": block_len, "range": block_range})
        nodelist.append(data)

    edgelist = []
    for k, v in coarse_edges.items():
        data = (k[0], k[1], {"freq": v})
        edgelist.append(data)

    # print(nodelist) ####
    # print(edgelist) ####

    graph = nx.DiGraph()
    graph.add_nodes_from(nodelist)
    graph.add_edges_from(edgelist)
    # print(graph) ####

    components = nx.weakly_connected_components(graph)
    
    plot_dir = os.path.join(results_dir, "graph2")
    os.makedirs(plot_dir, exist_ok=True)
    for ind, c in enumerate(components):
        if len(c) <= 2:
            continue
        print(len(c)) ####
        result_path = os.path.join(plot_dir, f"coarse_{ind:02d}.svg")
        result_path_nl = os.path.join(plot_dir, f"coarse_{ind:02d}_nl.svg")
        plot_component(graph, c, result_path, result_path_nl)

if __name__ == "__main__":
    data_dir = "/oak/stanford/groups/wjg/atwang/ecdna/data"

    results_dir_base = "/oak/stanford/groups/wjg/atwang/ecdna/results"

    in_path = os.path.join(data_dir, "COLO320DM_gDNA_nanopore_guppy_4.4_graph.pickle")
    results_dir = os.path.join(results_dir_base, "COLO320DM_gDNA_nanopore_guppy_4.4")
    os.makedirs(results_dir, exist_ok=True)
    view_graph_coarse(in_path, results_dir)

    in_path = os.path.join(data_dir, "PC3_gDNA_combined_graph.pickle")
    results_dir = os.path.join(results_dir_base, "PC3_gDNA_combined")
    os.makedirs(results_dir, exist_ok=True)
    view_graph_coarse(in_path, results_dir)