import pickle
import os

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