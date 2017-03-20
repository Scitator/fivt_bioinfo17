import click
import collections
import itertools


def dna_reconstruction(k, d, kmers):
    path1 = [x[0] for x in kmers]
    path2 = [x[1] for x in kmers]
    result1 = "".join([x[0] for x in path1] + [path1[-1][1:]])
    result2 = "".join([x[0] for x in path2] + [path2[-1][1:]])
    result = result1[:k+d] + result2
    return result


@click.command()
@click.option(
    "--fin",
    type=str,
    default="problem15_input.txt")
def main(fin):
    with open(fin) as fin:
        data = fin.readlines()
    k, d = tuple(map(int, data[0].split()))
    kmers = list(map(lambda x: x.replace("\n", "").split("|"), data[1:]))
    print(dna_reconstruction(k, d, kmers))


if __name__ == '__main__':
    main()
