import numpy as np
import pandas as pd
import click
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


def make_profile(motifs, pseudocounts=True):
    motifs = np.array([x for x in motifs]).reshape((len(motifs), -1))
    profiles = []
    for i in range(len(dnas)):
        dna_i = motifs == i
        profiles.append(dna_i.sum(axis=0) / float(motifs.shape[0]))
    profiles = np.vstack(profiles)
    if pseudocounts:
        profiles += np.ones_like(profiles)
    return profiles


def score(motifs):
    profile = make_profile(motifs).T
    consensus = profile.argmax(axis=1)
    result = np.sum((motifs - consensus) != 0)
    return result


def most_prob_row_kmer(row, profile, k, use_probs=True):
    row_motif_probs = []
    for i in range(len(row) - k + 1):
        motif = one_hot_encoding(row[i:i+k], 4)
        row_motif_probs.append(np.product(np.sum(motif * profile, axis=1)))
    if not use_probs:
        most_prob_motif_index = np.argmax(row_motif_probs)
    else:
        row_motif_probs = np.array(row_motif_probs, dtype=np.float32) / float(sum(row_motif_probs))
        most_prob_motif_index = np.random.choice(len(row) - k + 1, p=row_motif_probs)
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


def gibbs_sampler(k, t, N, dna):
    motifs = np.vstack(list(map(lambda x: random_kmer(x, k), dna)))
    best_motifs = motifs.copy()
    for i in range(N):
        index = np.random.randint(t)
        profile = make_profile(np.delete(motifs, index, axis=0)).T / float(t)
        motif_i = most_prob_row_kmer(dna[index], profile, k)
        motifs[index] = motif_i
        if score(motifs) < score(best_motifs):
            best_motifs = motifs.copy()
    return best_motifs, score(best_motifs)


def gibbs_sampler_runner(k, t, N, dna, n_jobs=-1, t_max=40, seed=42):
    # np.random.seed(seed)
    dna = np.array([dna2id[x] for x in dna]).reshape((t, -1))
    mofits = Parallel(n_jobs)(t_max * [delayed(gibbs_sampler)(k, t, N, dna)])
    # mofits = [gibbs_sampler(k, t, N, dna)]
    mofits, scores = map(np.array, zip(*mofits))
    best_motifs = mofits[scores.argmin()]
    result = ["".join(id2dna[x] for x in motif) for motif in best_motifs]
    return "\n".join(result)


@click.command()
@click.option(
    "--fin",
    type=str,
    default="problem9_input.tsv")
def main(fin):
    df = pd.read_csv(fin, sep="\t")
    assert all(x in df.columns.values.tolist() for x in ["k", "t", "N", "dna"])
    for i, row in df.iterrows():
        print(gibbs_sampler_runner(row["k"], row["t"], row["N"], row["dna"]))


if __name__ == '__main__':
    main()
