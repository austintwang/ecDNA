import pickle
import os


def gather_distal_breaks(reads_distal, gene_contig, gene_strand, gene_pos, radius, overlap_thresh):
    min_pos = gene_pos - radius
    max_pos = gene_pos + radius

    in_breaks = []
    out_breaks = []

    for k, v in reads_distal.items():
        for i in range(len(v)):
            if (v[i][3] != gene_contig) or (v[i][4] < min_pos) or (v[i][5] > max_pos):
                continue
            
            if i > 0:
                overlap = (min(v[i-1][5], v[i][5]) - max(v[i-1][4], v[i][4])) / min(v[i][5] - v[i][4], v[i-1][5] - v[i-1][4])
                if overlap <= overlap_thresh:
                    a_contig = v[i][3]
                    a_strand = v[i][2]
                    a_pos = v[i][4] if a_strand == "+" else v[i][5]
                    b_contig = v[i-1][3]
                    b_strand = v[i-1][2]
                    b_pos = v[i-1][5] if b_strand == "+" else v[i-1][4]

                    if gene_strand == "+" and a_strand == "+":
                        in_breaks.append((b_contig, b_strand, b_pos, a_contig, a_strand, a_pos, k),)
                    elif gene_strand == "-" and b_strand == "-":
                        out_breaks.append((a_contig, a_strand, a_pos, b_contig, b_strand, b_pos, k),)

            if i < len(v) - 1:
                overlap = (min(v[i][5], v[i+1][5]) - max(v[i][4], v[i+1][4])) / min(v[i+1][5] - v[i+1][4], v[i][5] - v[i][4])
                if overlap <= overlap_thresh:
                    a_contig = v[i][3]
                    a_strand = v[i][2]
                    a_pos = v[i][4] if a_strand == "-" else v[i][5]
                    b_contig = v[i+1][3]
                    b_strand = v[i+1][2]
                    b_pos = v[i+1][5] if b_strand == "-" else v[i+1][4]

                    if gene_strand == "+" and a_strand == "+":
                        out_breaks.append((a_contig, a_strand, a_pos, b_contig, b_strand, b_pos, k),)
                    elif gene_strand == "-" and b_strand == "-":
                        in_breaks.append((b_contig, b_strand, b_pos, a_contig, a_strand, a_pos, k),)

    return in_breaks, out_breaks

def group_breaks(in_breaks, out_breaks, bucket_size):
    in_grouped = {}
    for v in in_breaks:
        contig, strand, pos = v[0], v[1], v[2]
        bucket = (contig, strand, pos // bucket_size * bucket_size)
        in_grouped.setdefault(bucket, 0)
        in_grouped[bucket] += 1

    out_grouped = {}
    for v in out_breaks:
        contig, strand, pos = v[3], v[4], v[5]
        bucket = (contig, strand, pos // bucket_size * bucket_size)
        out_grouped.setdefault(bucket, 0)
        out_grouped[bucket] += 1

    return in_grouped, out_grouped


def analyze_loci(reads_path, out_dir, genes, radius, overlap_thresh, bucket_size):
    with open(reads_path, "rb") as reads_file:
        reads_repeat, reads_foldback, reads_distal, _ = pickle.load(reads_file)

    for name, contig, strand, pos in genes:
        in_breaks, out_breaks = gather_distal_breaks(reads_distal, contig, strand, pos, radius, overlap_thresh)
        in_grouped, out_grouped = group_breaks(in_breaks, out_breaks, bucket_size)
        for k in sorted(in_grouped): ####
            v = in_grouped[k]
            print(k, v) ####
        for k in sorted(out_grouped): ####
            v = out_grouped[k]
            print(k, v) ####

        res = (in_breaks, out_breaks)
        out_path = os.path.join(out_dir, f"{name}_breaks_distal.pickle")
        with open(out_path, "wb") as out_file:
            pickle.dump(res, out_file)

if __name__ == "__main__":
    data_dir = "/oak/stanford/groups/wjg/atwang/ecdna/data"
    reads_path = os.path.join(data_dir, "COLO320DM_gDNA_nanopore_guppy_4.4_cats.pickle")
    out_dir = os.path.join(data_dir, "genes")
    os.makedirs(out_dir, exist_ok=True)

    genes = [
        ("ENSG00000136997", "NC_000008.11", "+", 127735434)
    ]

    radius = 1e5
    overlap_thresh = 0.
    bucket_size = 1e4
    analyze_loci(reads_path, out_dir, genes, radius, overlap_thresh, bucket_size)