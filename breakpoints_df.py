import pickle
import os
import numpy as np
import pandas as pd

def parse_seq(seq, resolution):
    strand = seq[2]
    chrom = seq[3]
    start = seq[4]
    end = seq[5]
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
            a = pts[x][1]
            b = pts[x+1][0]
            br = (a, b)
            breaks.setdefault(br, 0)
            breaks[br] += 1

    return breaks

def get_breaks_df(breaks):
    breaks_lst = []
    for k, v in breaks.items():
        a, b = k
        chrom_from, strand_from, pos_from = a
        seq_from = (chrom_from, strand_from)
        chrom_to, strand_to, pos_to = b
        seq_to = (chrom_to, strand_to)
        freq = v
        breaks_lst.append([seq_from, chrom_from, strand_from, pos_from, seq_to, chrom_to, strand_to, pos_to, freq])

    cols = ["seq_from", "chrom_from", "strand_from", "seq_to", "pos_from", "chrom_to", "strand_to", "pos_to", "frequency"]
    breaks_df = pd.DataFrame.from_records(breaks_lst, columns=cols)
    return breaks_df

def breakpoints_df(in_path, out_path, resolution):
    with open(in_path, "rb") as in_file:
        seqs, lens = pickle.load(in_file)

    breaks = get_breaks(seqs, resolution)
    breaks_df = get_breaks_df(breaks)

    print(breaks_df) ####
    print(breaks_df[breaks_df["frequency"] >= 2]) ####

    res = breaks_df, breaks
    with open(out_path, "wb") as out_file:
        pickle.dump(res, out_file)

if __name__ == "__main__":
    data_dir = "/oak/stanford/groups/wjg/atwang/ecdna/data"
    in_path = os.path.join(data_dir, "COLO320DM_gDNA_nanopore_guppy_4.4_splits.pickle")
    out_path = os.path.join(data_dir, "COLO320DM_gDNA_nanopore_guppy_4.4_breaks_proc.pickle")

    resolution = 1e3
    breakpoints_df(in_path, out_path, resolution)