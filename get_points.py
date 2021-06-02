import pickle
import os
import numpy as np

def parse_seq(seq, resolution):
    strand = i[3]
    chrom = i[4]
    start = i[5]
    end = i[6]
    strand_sgn = 1 if strand == "+" else -1
    start_bucket = int(start // resolution) * resolution * strand_sgn
    end_bucket = int(end // resolution) * resolution * strand_sgn
    start_pt = (chrom, strand, start_bucket)
    end_pt = (chrom, strand, end_bucket)
    return start_pt, end_pt

def parse_points(seqs, resolution):
    points_set = set()
    # adj = set()
    for k, v in seqs.items():
        for i in v:
            start_pt, end_pt = parse_seq(i, resolution)
            points_set.add(start_pt)
            points_set.add(end_pt)
            # adj.add((start_pt, end_pt),)

    points = sorted(points_set)
    ranks = {val: ind for ind, val in points}

    return points, ranks

def merge_seqs(seqs, points, ranks, resolution):
    copy_num = np.zeros(len(points), dtype=int)
    seqs_break = {}
    seqs_break_start = {}
    seqs_break_end = {}
    for k, v in seqs.items():
        if len(v) == 1:
            i = v,
            start_pt, end_pt = parse_seq(i, resolution)
            start_idx = ranks[start_pt]
            end_idx = ranks[end_pt]
            copy_num[start_idx:end_idx] += 1
            continue

        seq_pts = []
        for i in v:
            start_pt, end_pt = parse_seq(i, resolution)
            start_idx = ranks[start_pt]
            end_idx = ranks[end_pt]
            seq_pts.extend(points[start_idx:end_idx])
            copy_num[start_idx] += 1
            copy_num[end_idx] += 1
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
    
    return block_break_fwd, block_break_bwd, break_block_fwd, break_block_bwd

def prune_graph():



def separate_repeats(seqs, overlap_thresh):
    reads_repeat = {}
    reads_foldback = {}
    reads_distal = {}

    for k, v in seqs.items():
        has_repeat = False
        has_foldback = False
        has_distal = False

        for i in range(1, len(v)):
            if v[i][3] != v[i-1][3]:
                has_distal = True
                continue
            overlap = (min(v[i-1][5], v[i][5]) - max(v[i-1][4], v[i][4])) / min(v[i][5] - v[i][4], v[i-1][5] - v[i-1][4])
            if overlap > overlap_thresh:
                if v[i][2] != v[i-1][2]:
                    has_foldback = True
                else:
                    has_repeat = True
            else:
                has_distal = True

        if has_repeat:
            reads_repeat[k] = v
        if has_foldback:
            reads_foldback[k] = v
        if has_distal:
            reads_distal[k] = v

    return reads_repeat, reads_foldback, reads_distal

def filter_breakpoints(in_path, out_path, overlap_thresh):
    with open(in_path, "rb") as in_file:
        seqs, lens = pickle.load(in_file)

    reads_repeat, reads_foldback, reads_distal = separate_repeats(seqs, overlap_thresh)

    for k, v in reads_distal.items(): ####
        print(k, v) ####
    print(len(reads_repeat), len(reads_foldback), len(reads_distal))

    res = reads_repeat, reads_foldback, reads_distal, lens
    with open(out_path, "wb") as out_file:
        pickle.dump(res, out_file)

if __name__ == "__main__":
    data_dir = "/oak/stanford/groups/wjg/atwang/ecdna/data"
    in_path = os.path.join(data_dir, "COLO320DM_gDNA_nanopore_guppy_4.4_splits.pickle")
    out_path = os.path.join(data_dir, "COLO320DM_gDNA_nanopore_guppy_4.4_cats.pickle")

    overlap_thresh = 0.5
    filter_breakpoints(in_path, out_path, overlap_thresh)