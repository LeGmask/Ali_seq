from tabulate import tabulate

from src.sequence import Sequence
from .fasta import readFastaMul


def getMatrix(seq1: Sequence, seq2: Sequence):
    return [[" " if i != j else "*" for j in seq2.sequence] for i in seq1.sequence]


def dotPlot(file):
    sequences = readFastaMul(file)
    matrix = getMatrix(sequences[0], sequences[1])
    print(f" {sequences[1].sequence}")
    for i, val in enumerate(matrix):
        print(sequences[0].sequence[i] + "".join(val))


def drawMatrix(matrix):
    print(tabulate(matrix, tablefmt="fancy_grid"))


def saveMatrix(matrix, filename):
    with open(filename, "w") as f:
        f.write(tabulate(matrix, tablefmt="fancy_grid"))
