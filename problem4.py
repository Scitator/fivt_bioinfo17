"""
We defined a mismatch in “Compute the Hamming Distance Between Two Strings”.
    We now generalize “Find the Most Frequent Words in a String” to incorporate mismatches as well.

Given strings Text and Pattern as well as an integer d, we define Countd(Text, Pattern) as
    the total number of occurrences of Pattern in Text with at most d mismatches.
    For example, Count1(AACAAGCTGATAAACATTTAAAGAG, AAAAA) = 4
    because AAAAA appears four times in this string with at most one mismatch:
    AACAA, ATAAA, AAACA, and AAAGA. Note that two of these occurrences overlap.

A most frequent k-mer with up to d mismatches in
    Text is simply a string Pattern maximizing Countd(Text, Pattern) among all k-mers.
    Note that Pattern does not need to actually appear as a substring of Text;
    for example, AAAAA is the most frequent 5-mer with 1 mismatch in AACAAGCTGATAAACATTTAAAGAG,
    even though AAAAA does not appear exactly in this string.
    Keep this in mind while solving the following problem.

Frequent Words with Mismatches Problem

Find the most frequent k-mers with mismatches in a string.

Given: A string Text as well as integers k and d.

Return: All most frequent k-mers with up to d mismatches in Text.

Sample Dataset

ACGTTGCATGTCGCATGATGCATGAGAGCT
4 1
Sample Output

GATG ATGC ATGT
"""

import pandas as pd
import click
import collections
import itertools

mis_table = {
    "A": "CGT",
    "C": "AGT",
    "G": "ACT",
    "T": "ACG"
}


def kmer_variations(kmer, d):
    variations = [kmer]

    for delta in range(1, d + 1):
        for delta_idx in itertools.combinations(range(len(kmer)), delta):
            for variants in itertools.product(*[mis_table[kmer[i]] for i in delta_idx]):
                result = list(kmer)
                for idx, val in zip(delta_idx, variants):
                    result[idx] = val
                variations.append("".join(result))

    return variations


def approximate_pattern_matching(text, k, d):
    words = collections.defaultdict(int)
    for i in range(len(text) - k + 1):
        words[text[i:i + k]] += 1

    # now all mismatch
    variations = collections.defaultdict(int)
    for key, value in words.items():
        for variation in kmer_variations(key, d):
            variations[variation] += value

    maximum = max(variations.values())

    result = [key for key, value in variations.items() if value == maximum]
    return "\t".join(result)


@click.command()
@click.option(
    "--fin",
    type=str,
    default="problem4_input.tsv")
def main(fin):
    df = pd.read_csv(fin, sep="\t")
    assert all(x in df.columns.values.tolist() for x in ["text", "k", "d"])
    for i, row in df.iterrows():
        print(approximate_pattern_matching(row["text"], row["k"], row["d"]))


if __name__ == '__main__':
    main()
