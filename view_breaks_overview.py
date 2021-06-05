import pickle
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def plot_facet_breaks(df, freq_col, seq_order, result_path):
    # print(df) ####
    df_plot = pd.pivot_table(df, values=freq_col, index="seq_from", columns="seq_to", aggfunc=np.sum, fill_value=0)
    sns.set(style="whitegrid", font="Roboto")
    plt.figure(figsize=(50,50))
    sns.heatmap(df_plot, annot=True, robust=True)
    plt.savefig(result_path, bbox_inches='tight')
    plt.clf()
    plt.close()

def view_breaks(in_path, min_supports, results_dir):
    with open(in_path, "rb") as in_file:
        breaks_df, breaks, freqs_from, freqs_to = pickle.load(in_file)

    seq_order = sorted(list(set(breaks_df["seq_from"]) | set(breaks_df["seq_to"])))
    print(seq_order) ####
    facet_dir = os.path.join(results_dir, "genome_overview")
    os.makedirs(facet_dir, exist_ok=True)

    for min_support in min_supports:
        df_from = breaks_df[breaks_df["freq_from"] >= min_support]
        df_to = breaks_df[breaks_df["freq_to"] >= min_support]
        df_pair = breaks_df[breaks_df["freq_pair"] >= min_support]

        result_path_from = os.path.join(facet_dir, f"ov_from_s_{min_support}.svg")
        result_path_to = os.path.join(facet_dir, f"ov_to_s_{min_support}.svg")
        result_path_pair = os.path.join(facet_dir, f"ov_pair_s_{min_support}.svg") 

        plot_facet_breaks(df_from, "freq_from", seq_order, result_path_from)
        plot_facet_breaks(df_to, "freq_to", seq_order, result_path_to)
        plot_facet_breaks(df_pair, "freq_pair", seq_order, result_path_pair)

if __name__ == "__main__":
    data_dir = "/oak/stanford/groups/wjg/atwang/ecdna/data"
    in_path = os.path.join(data_dir, "COLO320DM_gDNA_nanopore_guppy_4.4_breaks_proc.pickle")

    results_dir = "/oak/stanford/groups/wjg/atwang/ecdna/results"
    os.makedirs(results_dir, exist_ok=True)

    min_supports = [1, 2]
    view_breaks(in_path, min_supports, results_dir)