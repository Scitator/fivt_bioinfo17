import click
import numpy as np


def breakpoints_number(permutation):
    return sum(map(
        lambda x, y: x - y != 1,
        permutation + [len(permutation) + 1],
        [0] + permutation))


@click.command()
@click.option(
    "--fin",
    type=str,
    default="problem27_input.txt")
def main(fin):
    with open(fin) as fin:
        data = fin.readlines()[0]
    permutation = list(map(int, data[1:-1].replace("\n", "").replace("+", "").split()))
    result = breakpoints_number(permutation)
    print(result)


if __name__ == '__main__':
    main()
