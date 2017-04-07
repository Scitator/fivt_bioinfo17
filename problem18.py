import click

char2mass = dict(G=57, A=71, S=87, P=97, V=99, T=101, C=103, I=113, L=113, N=114, D=115, K=128,
                 Q=128, E=129, M=131, H=137, F=147, R=156, Y=163, W=186)

mass = [57, 71, 87, 97, 99, 101, 103, 113, 114, 115, 128, 129, 131, 137, 147, 156, 163, 186]


def expand_list(peptides, masses):
    if len(peptides) == 0:
        return [([m], [0, m]) for m in masses]

    def combine_spectrum(peptide, mass):
        def extend_spectrum(masses, m):
            return masses + [m] + [(sum(peptide[0][i:]) + m) for i in range(len(peptide[0]))]
        return peptide[0] + [mass], extend_spectrum(peptide[1], mass)
    return [combine_spectrum(p, m) for p in peptides for m in masses]


def consistent(peptide, spectrum):
    for e in peptide:
        if peptide.count(e) > spectrum.count(e):
            return False
    return True


def CyclopeptideSequencing(spectrum):
    result = []
    curr_list = []
    done = False

    while not done:
        result = [c[0] for c in curr_list]
        curr_list = [cand for cand in expand_list(curr_list, mass) if consistent(cand[1], spectrum)]
        if len(curr_list) == 0:
            done = True

    return result


@click.command()
@click.option(
    "--fin",
    type=str,
    default="problem18_input.txt")
def main(fin):
    with open(fin) as fin:
        data = fin.readlines()
    spectrum = tuple(map(int, data[0].split()))
    result = CyclopeptideSequencing(spectrum)
    # and finally, I'm too lazy
    result = " ".join(["-".join(map(str, line)) for line in reversed(result)])
    print(result)


if __name__ == '__main__':
    main()
