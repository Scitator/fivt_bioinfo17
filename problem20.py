import click
import numpy as np

from problem19 import LeaderboardCyclopeptideSequencing


def spectrum_mass(spectrum, m):
    spectrum = sorted(spectrum)

    spectral_convolution = [
        spectrum[i] - spectrum[j]
        for i in range(len(spectrum))
        for j in range(i)]
    spectral_convolution.extend(spectrum)

    spectral_convolution = list(filter(lambda x: 57 <= x <= 200, spectral_convolution))
    spectral_convolution_counted = list(map(
        lambda x: (x, spectral_convolution.count(x)),
        set(spectral_convolution)))
    spectral_convolution_counted = sorted(
        spectral_convolution_counted, key=lambda x: x[1], reverse=True)

    threshold = spectral_convolution_counted[m][1] \
        if m < len(spectral_convolution_counted) \
        else -np.inf
    spectral_mass = [s[0] for s in spectral_convolution_counted if s[1] >= threshold]
    return spectral_mass


def ConvolutionCyclopeptideSequencing(m, n, spectrum):
    new_mass = spectrum_mass(spectrum, m)
    result = LeaderboardCyclopeptideSequencing(spectrum, n, new_mass)
    return result


@click.command()
@click.option(
    "--fin",
    type=str,
    default="problem20_input.txt")
def main(fin):
    with open(fin) as fin:
        data = fin.readlines()
    m = tuple(map(int, data[0].split()))[0]
    n = tuple(map(int, data[1].split()))[0]
    spectrum = list(map(int, data[2].split()))
    result = ConvolutionCyclopeptideSequencing(m, n, spectrum)
    print("-".join(str(i) for i in result[0]))


if __name__ == '__main__':
    main()
