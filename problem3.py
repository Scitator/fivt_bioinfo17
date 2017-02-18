"""
We say that a k-mer Pattern appears as a substring of Text with at most d mismatches
    if there is some k-mer substring Pattern' of Text having d or fewer mismatches with Pattern,
    i.e., HammingDistance(Pattern, Pattern') ≤ d.
    Our observation that a DnaA box may appear with slight variations leads to
        the following generalization of the Pattern Matching Problem.

Approximate Pattern Matching Problem

Find all approximate occurrences of a pattern in a string.

Given: Strings Pattern and Text along with an integer d.

Return: All starting positions where Pattern appears as a substring of Text with at most d mismatches.

Sample Dataset

ATTCTGGA
CGCCCGAATCCAGAACGCATTCCCATATTTCGGGACCACTGGCCTCCACGGTACGGACGTCAATCAAATGCCTAGCGGCTTGTGGTTTCTCCTACGCTCC
3
Sample Output

6 7 26 27 78
"""


import pandas as pd
import click


def approximate_pattern_atching(pattern, text, d):
    result = []
    n = len(pattern)
    for i in range(len(text) - n):
        delta = sum([x != y for x, y in zip(pattern, text[i:i+n])])
        if delta <= d:
            result.append(str(i))
    return "\t".join(result)


@click.command()
@click.option(
    "--fin",
    type=str,
    default="problem3_input.tsv")
def main(fin):
    df = pd.read_csv(fin, sep="\t")
    assert all(x in df.columns.values.tolist() for x in ["pattern", "text", "d"])
    for i, row in df.iterrows():
        print(approximate_pattern_atching(row["pattern"], row["text"], row["d"]))


if __name__ == '__main__':
    main()
