import numpy as np
import pandas as pd
import click
import collections
import itertools
from joblib import Parallel, delayed


def pattern_vs_strings_distance(pattern, dna):
    k = len(pattern)
    dna = dna.split()
    distance = 0.0
    for row in dna:
        hamming_distance = np.inf
        for i in range(len(row) - k + 1):
            kmer = row[i:i+k]
            curr_dist = np.sum([x != y for x, y in zip(kmer, pattern)])
            if hamming_distance > curr_dist:
                hamming_distance = curr_dist
        distance += hamming_distance
    result = [str(int(distance))]
    return "\n".join(result)


@click.command()
@click.option(
    "--fin",
    type=str,
    default="problem10_input.tsv")
def main(fin):
    df = pd.read_csv(fin, sep="\t")
    assert all(x in df.columns.values.tolist() for x in ["pattern", "dna"])
    for i, row in df.iterrows():
        print(pattern_vs_strings_distance(row["pattern"], row["dna"]))


if __name__ == '__main__':
    main()
