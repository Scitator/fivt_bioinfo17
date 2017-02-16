"""
Define the skew of a DNA string Genome, denoted Skew(Genome),
    as the difference between the total number of occurrences of 'G' and 'C' in Genome.
    Let Prefixi (Genome) denote the prefix (i.e., initial substring) of Genome of length i.
    For example, the values of Skew(Prefixi ("CATGGGCATCGGCCATACGCCCATGGGCATCGGCCATACGCC")) are:
        0 -1 -1 -1 0 1 2 1 1 1 0 1 2 1 0 0 0 0 -1 0 -1 -2

Minimum Skew Problem

Find a position in a genome minimizing the skew.

Given: A DNA string Genome.

Return: All integer(s) i minimizing Skew(Prefixi (Text)) over all values of i (from 0 to |Genome|).

Sample Dataset

CCTATCGGTGGATTAGCATGTCCCTGTACGTTTCGCCGCGAACTAGTTCACACGGCTTGATGGCAAATGGTTTTTCCGGCGACCGTAATCGTCCACCGAG
Sample Output

53 97
"""


import pandas as pd
import click

genome_dict = {
    "C": -1,
    "G": 1
}


def minimum_skew(genome):
    result_sum = 0

    def complex_update_func(x):
        nonlocal result_sum
        result_sum += genome_dict.get(x, 0)
        return result_sum

    result = list(map(lambda x: complex_update_func(x), genome))
    minimum = min(result)
    result = [str(i+1) for i, v in enumerate(result) if v == minimum]
    return "\t".join(result)


@click.command()
@click.option(
    "--fin",
    type=str,
    default="problem2_input.tsv")
def main(fin):
    df = pd.read_csv(fin, sep="\t")
    assert all(x in df.columns.values.tolist() for x in ["genome"])
    for i, row in df.iterrows():
        print(minimum_skew(row["genome"]))


if __name__ == '__main__':
    main()
