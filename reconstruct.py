import pickle
import os
import numpy as np

def parse_seq(seq, resolution):
    # print(seq) ####
    # try: ####
    strand = seq[2]
    # except Exception as e: ####
    #     print(seq) ####
    #     raise e ####
    chrom = seq[3]
    start = seq[4]
    end = seq[5]
    strand_sgn = 1 if strand == "+" else -1
    start_bucket = int(start // resolution) * resolution * strand_sgn
    end_bucket = int(end // resolution) * resolution * strand_sgn
    start_pt = (chrom, strand, start_bucket)
    end_pt = (chrom, strand, end_bucket)
    strand_rc = "-" if strand == "+" else "+"
    start_pt_rc = (chrom, strand, end_bucket * -1)
    end_pt_rc = (chrom, strand, start_bucket * -1)
    return start_pt, end_pt, start_pt_rc, end_pt_rc

def parse_points(seqs, resolution):
    points_set = set()
    # adj = set()
    for k, v in seqs.items():
        for i in v:
            start_pt, end_pt, start_pt_rc, end_pt_rc = parse_seq(i, resolution)
            points_set.add(start_pt)
            points_set.add(end_pt)
            points_set.add(start_pt_rc)
            points_set.add(end_pt_rc)
            # adj.add((start_pt, end_pt),)

    points = sorted(points_set)
    ranks = {val: ind for ind, val in enumerate(points)}

    return points, ranks

def merge_seqs(seqs, points, ranks, resolution):
    copy_num = np.zeros(len(points), dtype=int)
    seqs_break = {}
    seqs_break_start = {}
    seqs_break_end = {}
    for k, v in seqs.items():
        if len(v) == 1:
            i, = v
            start_pt, end_pt, start_pt_rc, end_pt_rc = parse_seq(i, resolution)
            start_idx = ranks[start_pt]
            end_idx = ranks[end_pt]
            copy_num[start_idx:end_idx] += 1
            start_idx_rc = ranks[start_pt_rc]
            end_idx_rc = ranks[end_pt_rc]
            copy_num[start_idx_rc:end_idx_rc] += 1
            continue

        seq_pts = []
        for i in v:
            start_pt, end_pt = parse_seq(i, resolution)
            start_idx = ranks[start_pt]
            end_idx = ranks[end_pt]
            seq_pts.extend(points[start_idx:end_idx])
            copy_num[start_idx] += 1
            copy_num[end_idx] += 1
            start_idx_rc = ranks[start_pt_rc]
            end_idx_rc = ranks[end_pt_rc]
            seq_pts.extend(points[start_idx_rc:end_idx_rc])
            copy_num[start_idx_rc] += 1
            copy_num[end_idx_rc] += 1
        seqs_break[k] = seq_pts
        seqs_break_start.setdefault(seq_pts[0], []).append(k)
        seqs_break_end.setdefault(seq_pts[-1], []).append(k)

    return copy_num, seqs_break, seqs_break_start, seqs_break_end

def build_graph(points, copy_num, seqs_break, seqs_break_start, seqs_break_end):
    blocks = []
    blocks_cn = []
    pt_block_map = {}
    block_ind = 0
    no_cov_ind = np.nonzero(copy_num == 0),
    for i in range(1, no_cov_ind.size):
        a = no_cov_ind[i-1] + 1
        b = no_cov_ind[i]
        if a == b:
            continue

        block = points[a:b]
        block_cn = copy_num[a:b]
        blocks.append(block)
        blocks_cn.append(block_cn)

        for b in block:
            pt_block_map[b] = block_ind

        block_ind += 1

    block_break_fwd = [{} for _ in range(len(blocks))]
    break_block_bwd = {}
    for k, v in seqs_break_start.items():
        block = pt_block_map[k]
        block_break_fwd[block][k] = v
        break_block_bwd[k] = block

    block_break_bwd = [{} for _ in range(len(blocks))]
    break_block_fwd = {}
    for k, v in seqs_break_end.items():
        block = pt_block_map[k]
        block_break_bwd[block][k] = v
        break_block_fwd[k] = block
    
    node_data = (blocks, blocks_cn, pt_block_map)
    edge_data = (block_break_fwd, block_break_bwd, break_block_fwd, break_block_bwd)
    return node_data, edge_data

def prune_graph(node_data, edge_data):
    blocks, blocks_cn, pt_block_map = node_data
    blocks_p = dict(enumerate(blocks))
    blocks_cn_p = dict(enumerate(blocks_cn))

    block_break_fwd, block_break_bwd, break_block_fwd, break_block_bwd = edge_data
    block_break_fwd_p = dict(enumerate(block_break_fwd))
    block_break_bwd_p = dict(enumerate(block_break_bwd))
    break_block_fwd_p = copy(break_block_fwd)
    break_block_bwd_p = copy(break_block_bwd)

    consider_set = set(range(len(blocks)))
    while consider_set:
        consider_set_next = set()
        delete_node_set = set()
        delete_edge_set = set()
        for i in consider_set:
            inbounds = block_break_bwd_p[i]
            outbounds = block_break_fwd_p[i]
            if inbounds and outbounds:
                min_inbound = min(inbounds.keys())
                max_outbound = max(outbounds.keys())
                if max_outbound >= min_inbound:
                    continue

            delete_node_set.add(i)
            in_deletes = set(inbounds.keys())
            out_deletes = set(outbounds.keys())
            delete_edge_set |= (in_deletes | out_deletes)

            in_consider = set(break_block_bwd_p[i] for i in in_deletes)
            out_consider = set(break_block_fwd_p[i] for i in out_deletes)
            consider_set_next |= (in_consider | out_consider)

        consider_set = consider_set_next - delete_node_set
        for n in delete_node_set:
            blocks_p.pop(n)
            blocks_cn_p.pop(n)
            block_break_fwd_p.pop(n)
            block_break_bwd_p.pop(n)
        for e in delete_edge_set:
            break_block_fwd_p.pop(e)
            break_block_bwd_p.pop(e)

    node_data = (blocks_p, blocks_cn_p, pt_block_map)
    edge_data = (block_break_fwd_p, block_break_bwd_p, break_block_fwd_p, break_block_bwd_p)
    return node_data, edge_data

def create_coarse_grained(node_data, edge_data):
    blocks, blocks_cn, pt_block_map = node_data
    block_break_fwd, block_break_bwd, break_block_fwd, break_block_bwd = edge_data

    coarse_edges = {}
    for k in blocks.keys():
        nxt = break_block_fwd[block_break_fwd[k]]
        edge = (k, nxt)
        coarse_edges.setdefault(edge)
        coarse_edges[edge] += 1

    return coarse_edges
        
def reconstruct_amplicons(in_path, out_path, resolution):
    with open(in_path, "rb") as in_file:
        seqs, lens = pickle.load(in_file)

    points, ranks = parse_points(seqs, resolution)
    copy_num, seqs_break, seqs_break_start, seqs_break_end = merge_seqs(seqs, points, ranks, resolution)
    print(copy_num[:10]) ####
    node_data, edge_data = build_graph(points, copy_num, seqs_break, seqs_break_start, seqs_break_end)
    print(len(node_data)) ####
    node_data_p, edge_data_p = prune_graph(node_data, edge_data)
    print(len(node_data_p)) ####
    coarse_edges = create_coarse_grained(node_data, edge_data)
    print(len(coarse_edges)) ####

    res = (node_data_p, edge_data_p, coarse_edges)
    with open(out_path, "wb") as out_file:
        pickle.dump(res, out_file)

if __name__ == "__main__":
    resolution = 1e3
    data_dir = "/oak/stanford/groups/wjg/atwang/ecdna/data"

    # in_path = os.path.join(data_dir, "COLO320DM_gDNA_nanopore_guppy_4.4_all.pickle")
    # out_path = os.path.join(data_dir, "COLO320DM_gDNA_nanopore_guppy_4.4_graph.pickle")

    # breakpoints_df(in_path, out_path, resolution)

    in_path = os.path.join(data_dir, "PC3_gDNA_combined_all.pickle")
    out_path = os.path.join(data_dir, "PC3_gDNA_combined_graph.pickle")

    reconstruct_amplicons(in_path, out_path, resolution)