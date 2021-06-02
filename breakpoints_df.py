import pickle
import os
import numpy as np
import pandas as pd

def parse_seq(seq, resolution):
    print(seq) ####
    strand = seq[3]
    chrom = seq[4]
    start = seq[5]
    end = seq[6]
    strand_sgn = 1 if strand == "+" else -1
    start_bucket = int(start // resolution) * resolution * strand_sgn
    end_bucket = int(end // resolution) * resolution * strand_sgn
    start_pt = (chrom, strand, start_bucket)
    end_pt = (chrom, strand, end_bucket)
    return start_pt, end_pt

def get_breaks(seqs, resolution):
    points_set = set()
    breaks = {}
    for k, v in seqs.items():
        pts = []
        for i in v:
            p = parse_seq(i, resolution)
            pts.append(p)
        for x in range(len(pts) - 1):
            a = pts[i][1]
            b = pts[i+1][0]
            br = (a, b)
            breaks.setdefault(br, 0)
            breaks[br] += 1

    return breaks

def get_breaks_df(breaks):
    breaks_lst = []
    for k, v in breaks.items():
        a, b = k
        chrom_from, strand_from, pos_from = a
        chrom_to, strand_to, pos_to = b
        freq = b
        breaks_lst.append([chrom_from, strand_from, pos_from, chrom_to, strand_to, pos_to, freq])

    cols = ["chrom_from", "strand_from", "pos_from", "chrom_to", "strand_to", "pos_to", "frequency"]
    breaks_df = pd.DataFrame.from_records(breaks_lst, cols=cols)
    return breaks_df

def breakpoints_df(in_path, out_path, resolution):
    with open(in_path, "rb") as in_file:
        seqs, lens = pickle.load(in_file)

    breaks = get_breaks(seqs, resolution)
    breaks_df = get_breaks_df(breaks)

    print(breaks_df) ####
    print(breaks_df[breaks_df["frequency"] >= 2]) ####

    res = breaks_hc, breaks
    with open(out_path, "wb") as out_file:
        pickle.dump(res, out_file)

if __name__ == "__main__":
    data_dir = "/oak/stanford/groups/wjg/atwang/ecdna/data"
    in_path = os.path.join(data_dir, "COLO320DM_gDNA_nanopore_guppy_4.4_splits.pickle")
    out_path = os.path.join(data_dir, "COLO320DM_gDNA_nanopore_guppy_4.4_breaks_df.pickle")

    resolution = 1e4
    breakpoints_df(in_path, out_path, resolution)