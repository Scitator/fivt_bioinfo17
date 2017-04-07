import click

char2mass = dict(G=57, A=71, S=87, P=97, V=99, T=101, C=103, I=113, L=113, N=114, D=115, K=128,
                 Q=128, E=129, M=131, H=137, F=147, R=156, Y=163, W=186)

correct_codon_dict = {
    "AAA": "K",
    "AAC": "N",
    "AAG": "K",
    "AAT": "N",
    "ACA": "T",
    "ACC": "T",
    "ACG": "T",
    "ACT": "T",
    "AGA": "R",
    "AGC": "S",
    "AGG": "R",
    "AGT": "S",
    "ATA": "I",
    "ATC": "I",
    "ATG": "M",
    "ATT": "I",
    "CAA": "Q",
    "CAC": "H",
    "CAG": "Q",
    "CAT": "H",
    "CCA": "P",
    "CCC": "P",
    "CCG": "P",
    "CCT": "P",
    "CGA": "R",
    "CGC": "R",
    "CGG": "R",
    "CGT": "R",
    "CTA": "L",
    "CTC": "L",
    "CTG": "L",
    "CTT": "L",
    "GAA": "E",
    "GAC": "D",
    "GAG": "E",
    "GAT": "D",
    "GCA": "A",
    "GCC": "A",
    "GCG": "A",
    "GCT": "A",
    "GGA": "G",
    "GGC": "G",
    "GGG": "G",
    "GGT": "G",
    "GTA": "V",
    "GTC": "V",
    "GTG": "V",
    "GTT": "V",
    "TAA": "X",
    "TAC": "Y",
    "TAG": "X",
    "TAT": "Y",
    "TCA": "S",
    "TCC": "S",
    "TCG": "S",
    "TCT": "S",
    "TGA": "X",
    "TGC": "C",
    "TGG": "W",
    "TGT": "C",
    "TTA": "L",
    "TTC": "F",
    "TTG": "L",
    "TTT": "F",
}

dna_dict = {
    "A": "T",
    "C": "G",
    "G": "C",
    "T": "A"}


def reverse_complement(dna):
    return "".join(reversed([dna_dict[x] for x in dna]))


def encoding(dna, peptide, k=3):
    results = []
    for f in range(k):
        translated = "".join(
            [correct_codon_dict[dna[i:i + k]]
             for i in range(f, int(len(dna)-f), k) if len(dna[i:i + k]) == k])
        results.extend(
            [dna[i * k + f:i * k + f + k * len(peptide)]
             for i in range(len(translated) - len(peptide) + 1)
             if translated[i:i + len(peptide)] == peptide])
    return results


def PeptideEncodingProblem(dna, peptide):
    results = []
    results.extend(encoding(dna, peptide))
    results.extend([reverse_complement(x) for x in encoding(reverse_complement(dna), peptide)])
    return "\n".join(results)


@click.command()
@click.option(
    "--fin",
    type=str,
    default="problem16_input.txt")
def main(fin):
    with open(fin) as fin:
        data = fin.readlines()
    dna = data[0].strip()
    peptide = data[1].strip()
    print(PeptideEncodingProblem(dna, peptide))


if __name__ == "__main__":
    main()
