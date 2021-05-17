import pickle
import os

def parse_paf(paf_path, q_cutoff):
    seqs = {}
    lens = {}
    with open(paf_path, "r") as paf_file:
        for line in paf_file:
            cols = line.strip.split()
            qual = int(line[11])
            if qual < q_cutoff or qual == 255:
                continue

            qname = line[0]
            qlen, qstart, qend = int(line[1]), int(line[2]), int(line[3])
            strand = line[4]
            rname = line[5]
            rstart, rend = int(line[6]), int(line[7])

            seq_data = (qstart, qend, strand, rname, rstart, rend)
            seqs.setdefault(qname, []).append(seq_data)

            lens[qname] = qlen

    return seqs, lens

def get_breakpoints(paf_path, out_path, q_cutoff):
    seqs, lens = parse_paf(paf_path, q_cutoff)
    seqs_filtered = {}
    lens_filtered = {}
    for k, v in seqs.items():
        if len(v) >= 2:
            s = sorted(v)
            print(s) ####
            seqs_filtered[k] = s
            lens_filtered[k] = lens[k]

    res = (seqs_filtered, lens_filtered)
    with open(out_path, "wb") as out_file:
        pickle.dump(res, out_file)

if __name__ == "__main__":
    data_dir = "/oak/stanford/groups/wjg/atwang/ecdna/data"
    paf_path = os.path.join(data_dir, "COLO320DM_gDNA_nanopore_guppy_4.4_map.paf")
    out_path = os.path.join(data_dir, "COLO320DM_gDNA_nanopore_guppy_4.4_splits.pickle")

    q_cutoff = 50
    get_breakpoints(paf_path, out_path, q_cutoff)