import pickle
import os

def parse_paf(paf_path, q_cutoff):
    seqs = {}
    lens = {}
    with open(paf_path, "r") as paf_file:
        for line in paf_file:
            cols = line.strip().split()
            qual = int(cols[11])
            if qual < q_cutoff or qual == 255:
                continue

            qname = cols[0]
            qlen, qstart, qend = int(cols[1]), int(cols[2]), int(cols[3])
            strand = cols[4]
            rname = cols[5]
            rstart, rend = int(cols[7]), int(cols[8])

            seq_data = (qstart, qend, strand, rname, rstart, rend)
            seqs.setdefault(qname, []).append(seq_data)

            lens[qname] = qlen

    return seqs, lens

def get_breakpoints(paf_path, out_path, out_path_all, q_cutoff):
    seqs, lens = parse_paf(paf_path, q_cutoff)
    seqs_filtered = {}
    lens_filtered = {}
    for k, v in seqs.items():
        if len(v) >= 2:
            s = sorted(v)
            print(s) ####
            seqs_filtered[k] = s
            lens_filtered[k] = lens[k]

    print(len(seqs), len(seqs_filtered))
    res = (seqs_filtered, lens_filtered)
    with open(out_path, "wb") as out_file:
        pickle.dump(res, out_file)
    res_all = (seqs, lens)
    with open(out_path_all, "wb") as out_file:
        pickle.dump(res_all, out_file)

if __name__ == "__main__":
    data_dir = "/oak/stanford/groups/wjg/atwang/ecdna/data"
    q_cutoff = 50

    paf_path = os.path.join(data_dir, "COLO320DM_gDNA_nanopore_guppy_4.4_map.paf")
    out_path = os.path.join(data_dir, "COLO320DM_gDNA_nanopore_guppy_4.4_splits.pickle")
    out_path_all = os.path.join(data_dir, "COLO320DM_gDNA_nanopore_guppy_4.4_all.pickle")

    get_breakpoints(paf_path, out_path, out_path_all, q_cutoff)

    # paf_path = os.path.join(data_dir, "PC3_gDNA_combined_map.paf")
    # out_path = os.path.join(data_dir, "PC3_gDNA_combined_splits.pickle")
    # out_path_all = os.path.join(data_dir, "PC3_gDNA_combined_all.pickle")

    # get_breakpoints(paf_path, out_path, out_path_all, q_cutoff)