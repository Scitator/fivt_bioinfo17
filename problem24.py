import click
import numpy as np


def load_matrix():
    with open("./PAM250.txt") as fin:
        data = [line.strip().split() for line in fin.readlines()]
        score_matrix = {(item[0], item[1]): int(item[2]) for item in data}
        return score_matrix


def alignment(s1, s2, scoring_matrix, sigma=5):
    (height, width) = (len(s1) + 1, len(s2) + 1)

    forward = np.zeros((height, width), dtype=np.int32)
    backward = np.zeros((height, width), dtype=np.int32)
    # forward[:, 0] = np.linspace(0, height-1, height) * (-sigma)
    # forward[0, :] = np.linspace(0, width-1, width) * (-sigma)

    for i in range(1, len(s1)+1):
        for j in range(1, len(s2)+1):
            scores = [
                forward[i-1, j] - sigma,
                forward[i, j-1] - sigma,
                forward[i-1, j-1] + scoring_matrix[s1[i - 1], s2[j - 1]],
                0]
            forward[i, j] = max(scores)
            backward[i, j] = scores.index(forward[i, j])

    max_score = str(forward.max())
    l1, l2 = np.unravel_index(forward.argmax(), forward.shape)
    s1, s2 = s1[:l1], s2[:l2]

    add_indel = lambda word, i: word[:i] + '-' + word[i:]
    while l1 * l2 != 0 and backward[l1, l2] != 3:
        if backward[l1, l2] == 0:
            l1 -= 1
            s2 = add_indel(s2, l2)
        elif backward[l1, l2] == 1:
            l2 -= 1
            s1 = add_indel(s1, l1)
        elif backward[l1, l2] == 2:
            l1 -= 1
            l2 -= 1

    # for _ in range(l1):
    #     s2 = add_indel(s2, 0)
    # for _ in range(l2):
    #     s1 = add_indel(s1, 0)
    s1, s2 = s1[l1:], s2[l2:]

    return "\n".join([max_score, s1, s2])


@click.command()
@click.option(
    "--fin",
    type=str,
    default="problem24_input.txt")
def main(fin):
    with open(fin) as fin:
        data = fin.readlines()
    s1 = data[0].strip()
    s2 = data[1].strip()
    result = alignment(s1, s2, load_matrix())
    print(result)


if __name__ == '__main__':
    main()
