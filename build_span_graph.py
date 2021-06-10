import pickle
import os
import numpy as np

def build_graph(pt_data, node_data, edge_data):
    points, ranks, copy_num, seqs_break = pt_data
    blocks, blocks_cn, pt_block_map = node_data
    block_break_fwd, block_break_bwd, break_block_fwd, break_block_bwd = edge_data

    node_data = []
    for k in break_block_fwd.keys():
        pts = seqs_break[k]
        start = pts[0]
        end = pts[-1]
        data = (k, {"start_pt": start, "end_pt": end})
        node_data.append(data)

    edge_data = []
    for k, v in break_block_fwd.items():
        block, pt = v
        nexts = [(read, pt_n) for read, pt_n in block_break_fwd.get[block] if (pt_n > pt and read in break_block_fwd)]
        [k, read, {} for read, pt_n in nexts]
        for read, pt_n in nexts:
            r1 = k
            r2 = read
            pt_i1 = ranks[pt]
            pt_i2 = ranks[pt_n]
            min_cn = min(copy_num[pt_i1:pt_i2+1])
            data = (r1, r2, {"min_cn": min_cn})
            edge_data.append(data)

    return node_data, edge_data

def build_span_graph(in_path, out_path):
    with open(in_path, "rb") as in_file:
        pt_data, node_data, edge_data, coarse_edges = pickle.load(in_file)

    node_data, edge_data = build_graph(pt_data, node_data, edge_data)
    print(len(node_data), len(edge_data)) ####

    res = (node_data, edge_data)
    with open(out_path, "wb") as out_file:
        pickle.dump(res, out_file)

if __name__ == "__main__":
    data_dir = "/oak/stanford/groups/wjg/atwang/ecdna/data"

    in_path = os.path.join(data_dir, "COLO320DM_gDNA_nanopore_guppy_4.4_graph.pickle")
    out_path = os.path.join(data_dir, "COLO320DM_gDNA_nanopore_guppy_4.4_span.pickle")

    build_span_graph(in_path, out_path)

    in_path = os.path.join(data_dir, "PC3_gDNA_combined_graph.pickle")
    out_path = os.path.join(data_dir, "PC3_gDNA_combined_span.pickle")

    build_span_graph(in_path, out_path)