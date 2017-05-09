import click
import numpy as np


def greedy_sorting(permutation):
    permutations = []
    perm_index_fn = lambda perm, var: list(map(abs, perm)).index(var)
    perm_sort_fn = lambda perm, i, j: \
        perm[:i] + list(map(lambda x: -x, perm[i:j + 1][::-1])) + perm[j + 1:]

    i = 0
    while i < len(permutation):
        if permutation[i] == i + 1:
            i += 1
        else:
            permutation = perm_sort_fn(
                permutation,
                i, perm_index_fn(permutation, i + 1))
            permutations.append(permutation)

    return permutations


@click.command()
@click.option(
    "--fin",
    type=str,
    default="problem26_input.txt")
def main(fin):
    with open(fin) as fin:
        data = fin.readlines()[0]
    permutation = list(map(int, data[1:-1].replace("\n", "").split()))
    permutations = greedy_sorting(permutation)

    pretty_value_fn = lambda value: ['-', '+'][value > 0] + str(abs(value))

    result = "\n".join(
        ['(' + ' '.join([pretty_value_fn(value) for value in perm]) + ')'
         for perm in permutations])

    print(result)


if __name__ == '__main__':
    main()
