import numpy as np
import pandas as pd
import click
import collections
import itertools
from joblib import Parallel, delayed


dnas = "ACGT"

dna2id = {ch:i for i, ch in enumerate(dnas)}
id2dna = {i:ch for i, ch in enumerate(dnas)}


def one_hot_encoding(labels, n_label):
    n_batch = len(labels)
    one_hot_label = np.zeros((n_batch, n_label), dtype=np.float32)
    for idx, lbl in enumerate(labels):
        one_hot_label[idx, int(lbl)] = 1.0
    return one_hot_label


def make_profile(motifs, rseudocounts=True):
    motifs = np.array([x for x in motifs]).reshape((len(motifs), -1))
    profiles = {}
    for i in range(len(dnas)):
        dna_i = motifs == i
        profiles[i] = dna_i.sum(axis=0) / float(motifs.shape[0])
    profiles = np.vstack([value for key, value in profiles.items()])
    if rseudocounts:
        profiles += np.ones_like(profiles)
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


def random_kmer(row, k):
    index = np.random.randint(len(row) - k + 1)
    return row[index:index+k]


def make_motifs(profile, dna, k):
    motifs = []
    for row in dna:
        most_prob_motif = most_prob_row_kmer(row, profile, k)
        motifs.append(most_prob_motif)
    motifs = np.array(motifs)
    return motifs


def randomized_motif_search(k, t, dna):
    """
    - great time to start write docs!
    - No.
    """
    motifs = np.vstack(list(map(lambda x: random_kmer(x, k), dna)))
    best_motifs = motifs.copy()
    while True:
        profile = make_profile(motifs).T
        motifs = make_motifs(profile, dna, k)
        curr_score = score(best_motifs)
        new_score = score(motifs)
        if new_score < curr_score:
            best_motifs = motifs
        else:
            return best_motifs, curr_score


def randomized_motif_search_runner(k, t, dna, t_max=1000, seed=42, n_jobs=-1):
    np.random.seed(seed)
    dna = np.array([dna2id[x] for x in dna]).reshape((t, -1))
    mofits = Parallel(n_jobs)(t_max * [delayed(randomized_motif_search)(k, t, dna)])
    mofits, scores = map(np.array, zip(*mofits))
    best_motifs = mofits[scores.argmin()] 
    result = ["".join(id2dna[x] for x in motif) for motif in best_motifs]
    return "\n".join(result)


@click.command()
@click.option(
    "--fin",
    type=str,
    default="problem8_input.tsv")
def main(fin):
    df = pd.read_csv(fin, sep="\t")
    assert all(x in df.columns.values.tolist() for x in ["k", "t", "dna"])
    for i, row in df.iterrows():
        print(randomized_motif_search_runner(row["k"], row["t"], row["dna"]))


if __name__ == '__main__':
    main()
