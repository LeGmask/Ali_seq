from src.sequence import Sequence
from .fasta import readFastaMul

# def dotPlot(file):
#     seq = readFastaMul(file)
#     P1 = seq[0].sequence
#     P2 = seq[1].sequence

#     c = " "
#     for i in P1:
#         c += i
#     print(c)
#     for i in P2:
#         l = f'{i}'
#         for j in P1:
#             l += "*" if i == j else " "
#         print(l)


def getMatrix(seq1: Sequence, seq2: Sequence):
    return [[" " if i != j else "*" for j in seq2.sequence] for i in seq1.sequence]


def dotPlot(file):
    sequences = readFastaMul(file)
    matrix = getMatrix(sequences[0], sequences[1])
    print(f" {sequences[1].sequence}")
    for i, val in enumerate(matrix):
        print(sequences[0].sequence[i] + "".join(val))
