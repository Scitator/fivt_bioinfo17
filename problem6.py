import numpy as np
import pandas as pd
import click
import collections
import itertools


dnas = "ACGT"

dna2id = {ch:i for i, ch in enumerate(dnas)}
id2dna = {i:ch for i, ch in enumerate(dnas)}


def one_hot_encoding(labels, n_label, muptiplixator=0.0):
    n_batch = len(labels)
    one_hot_label = np.ones((n_batch, n_label), dtype=np.float32) * muptiplixator
    for idx, lbl in enumerate(labels):
        one_hot_label[idx, int(lbl)] = 1.0
    return one_hot_label


def make_profile(motifs):
    motifs = np.array([x for x in motifs]).reshape((len(motifs), -1))
    profiles = {}
    for i in range(len(dnas)):
        dna_i = motifs == i
        profiles[i] = dna_i.sum(axis=0) / float(motifs.shape[0])
    profiles = np.vstack([value for key, value in profiles.items()])
    return profiles


def score(motifs):
    profile = make_profile(motifs).T
    consensus = profile.argmax(axis=1)
    result = np.sum((motifs - consensus) != 0)
    return result


def most_prob_row_kmer(row, profile, k):
    row_motif_probs = []
    for i in range(len(row) - k + 1):
        motif = one_hot_encoding(row[i:i+k], 4)
        row_motif_probs.append(np.product(np.sum(motif * profile, axis=1)))
    most_prob_motif_index = np.argmax(row_motif_probs)
    most_prob_motif = row[most_prob_motif_index:most_prob_motif_index+k]
    return most_prob_motif


def greedy_motif_search(k, t, dna):
    """
    - great time to start write docs!
    - No.
    """
    dna = np.array([dna2id[x] for x in dna]).reshape((t, -1))
    best_motifs = dna[:, :k]
    row_len = dna.shape[1]
    for i in range(row_len - k + 1):
        motifs = [dna[0, i:i+k]]
        for row in dna[1:]:
            profile = make_profile(motifs).T
            most_prob_motif = most_prob_row_kmer(row, profile, k)
            motifs.append(most_prob_motif)
        motifs = np.array(motifs)
        if score(motifs) < score(best_motifs):
            best_motifs = motifs
    result = ["".join(id2dna[x] for x in motif) for motif in best_motifs]
    return "\n".join(result)


@click.command()
@click.option(
    "--fin",
    type=str,
    default="problem6_input.tsv")
def main(fin):
    df = pd.read_csv(fin, sep="\t")
    assert all(x in df.columns.values.tolist() for x in ["k", "t", "dna"])
    for i, row in df.iterrows():
        print(greedy_motif_search(row["k"], row["t"], row["dna"]))


if __name__ == '__main__':
    main()
