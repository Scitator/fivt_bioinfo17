"""
Given integers L and t, a string Pattern forms an (L, t)-clump inside a (larger) string Genome
    if there is an interval of Genome of length L in which Pattern appears at least t times.
    For example, TGCATGCA forms a (25,3)-clump in the following Genome:
    gatcagcataagggtcccTGCAATGCATGACAAGCCTGCAgttgttttacgatcagcataagggtcccTGCAATGCATGACAAGCCTGCAgttgttttac.

Clump Finding Problem

Find patterns forming clumps in a string.

Given: A string Genome, and integers k, L, and t.

Return: All distinct k-mers forming (L, t)-clumps in Genome.

Sample Dataset

CGGACTCGACAGATGTGAAGAAATGTGAAGACTGAGTGAAGAGAAGAGGAAACACGACACGACATTGCGACATAATGTACGAATGTAATGTGCCTATGGC
5 75 4
Sample Output

CGACA GAAGA AATGT
"""


import pandas as pd
import click


def clump_finding(genome, k, L, t):
    result = set()
    for i in range(len(genome) - L):
        count = {}
        for j in range(L - k):
            kmer = genome[i+j:i+j+k]
            count[kmer] = count.get(kmer, 0) + 1
        for kmer in count:
            if count[kmer] >= t:
                result.add(kmer)
    result = list(result)
    return "\t".join(result)


@click.command()
@click.option(
    "--fin",
    type=str,
    default="problem1_input.tsv")
def main(fin):
    df = pd.read_csv(fin, sep="\t")
    assert all(x in df.columns.values.tolist() for x in ["genome", "k", "L", "t"])
    for i, row in df.iterrows():
        print(clump_finding(row["genome"],  row["k"], row["L"], row["t"]))


if __name__ == '__main__':
    main()
