import click
import numpy as np


def manhattan_tourist(n, m, down, right):
    path = np.zeros((n + 1, m + 1), dtype=int)

    for i in range(1, n + 1):
        path[i][0] = path[i - 1][0] + down[i - 1][0]
    for j in range(1, m + 1):
        path[0][j] = path[0][j - 1] + right[0][j - 1]

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            path[i][j] = max(path[i - 1][j] + down[i - 1][j], path[i][j - 1] + right[i][j - 1])

    return path[n][m]


@click.command()
@click.option(
    "--fin",
    type=str,
    default="problem22_input.txt")
def main(fin):
    with open(fin) as fin:
        data = fin.readlines()
    n, m = tuple(map(int, data[0].split()))
    down = [list(map(int, data[i].strip().split())) for i in range(1, 1 + n)]
    right = [list(map(int, data[i].strip().split())) for i in range(2 + n, 3 + n + n)]
    result = manhattan_tourist(n, m, down, right)
    print(result)


if __name__ == '__main__':
    main()
