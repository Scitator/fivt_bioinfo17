import click
import numpy as np
from collections import defaultdict

mass = [57, 71, 87, 97, 99, 101, 103, 113, 114, 115, 128, 129, 131, 137, 147, 156, 163, 186]


def expand_list(peptides, masses):
    if len(peptides) == 0:
        return [[m] for m in masses]

    return [p + [m] for p in peptides for m in masses]


def count_dict(peptide):
    peptide_dct = defaultdict(int)
    for p in peptide:
        peptide_dct[p] += 1
    return peptide_dct


def score(peptide_dct, spectrum_dct):
    return sum([min(peptide_dct[k], spectrum_dct[k]) for k in peptide_dct.keys()])


def cyclopeptide(peptide):
    def subpeptide(peptide, pos, length):
        if pos + length <= len(peptide):
            return peptide[pos: pos+length]
        else:
            return peptide[pos:] + peptide[:length + pos - len(peptide)]
    return [subpeptide(peptide, p, l)
            for p in range(len(peptide))
            for l in range(1, len(peptide) + 1)]


def cyclospectrum(peptide):
    return sorted([sum(p) for p in cyclopeptide(peptide)] + [0])


def cut(leaderboard, spectrum_dct, n):
    scores = list(map(
        lambda x: (x, score(count_dict(cyclospectrum(x)), spectrum_dct)),
        leaderboard))
    scores = sorted(scores, key=lambda x: x[1], reverse=True)
    threshold = scores[n][1] if n < len(scores) else -np.inf
    leaders = [s[0] for s in scores if s[1] >= threshold]
    return leaders


def LeaderboardCyclopeptideSequencing(spectrum, n, mass_replacer=None):
    used_mass = (mass_replacer or mass)
    spectrum_dct = count_dict(spectrum)
    leaderboard = []
    leader_peptide = (0, 0)
    done = False

    while not done:
        leaderboard = expand_list(leaderboard, used_mass)
        # print('len(expanded): '.format(len(leaderboard)))
        leaderboard = cut(leaderboard, spectrum_dct, n)

        leaderboard = [p for p in leaderboard if sum(p) <= max(spectrum)]
        # print('len(cutted): '.format(len(leaderboard)))

        if len(leaderboard) == 0:
            done = True
            break

        best_score = score(count_dict(cyclospectrum(leaderboard[0])), spectrum_dct)
        if best_score > leader_peptide[1]:
            leader_peptide = (leaderboard[0], best_score)
        # print('leader_peptide: {}: {}'.format(
        #     leader_peptide[1], '-'.join(str(i) for i in leader_peptide[0])))





    return leader_peptide


@click.command()
@click.option(
    "--fin",
    type=str,
    default="problem19_input.txt")
def main(fin):
    with open(fin) as fin:
        data = fin.readlines()
    n = tuple(map(int, data[0].split()))[0]
    spectrum = tuple(map(int, data[1].split()))
    result = LeaderboardCyclopeptideSequencing(spectrum, n)
    print("-".join(str(i) for i in result[0]))


if __name__ == '__main__':
    main()
