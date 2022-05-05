from typing import List

from src.sequence import Sequence
from .fasta import readFastaMul

PlotMatrix = List[List[str]]


def getMatrix(seq1: Sequence, seq2: Sequence) -> PlotMatrix:
    """
    The getMatrix function takes two sequences and returns a matrix of the
    alignment scores between each position in the two sequences.  The first sequence
    is vertical, and the second is horizontal.  If it's a match

    :param seq1:Sequence: Specify that the first parameter is a sequence
    :param seq2:Sequence: Specify the second sequence to be compared with seq2
    :return: A matrix of size len(seq2) + 1 by len(seq2) + 1
    """
    return [[" " if i != j else "*" for j in seq2.sequence] for i in seq1.sequence]


def dotPlot(file):
    """
    The dotPlot function takes two sequences and returns a matrix of the dot plot between them.
    The second sequence is printed along the top, and the first sequence is printed along
    the side. The matrix contains all of the matches between each pair of letters in both
    sequences.

    :param file: Open the file
    :return: A matrix of the dotplot
    """
    sequences = readFastaMul(file)
    matrix = getMatrix(sequences[0], sequences[1])
    print(f" {sequences[1].sequence}")
    for i, val in enumerate(matrix):
        print(sequences[0].sequence[i] + "".join(val))


def drawMatrix(matrix):
    """
    The drawMatrix function draws a matrix of numbers to the screen.
    It takes one argument, a list of lists (a 2D array).
    The function will print each row on its own line.

    :param matrix: Specify the matrix to be drawn
    :return: A string representation of the matrix
    """
    return "".join("".join(i) + "\n" for i in matrix)


def saveMatrix(matrix, filename):
    """
    The saveMatrix function takes a matrix and saves it to the specified file.

    :param matrix: Save the matrix to a file
    :param filename: Specify the name of the file to which you want to save your matrix
    :return: The number of values written to the file
    """
    with open(filename, "w") as f:
        f.write(drawMatrix(matrix))


def saveToFile(string, filename) -> None:
    """
    The saveToFile function saves a string to the specified file.

    :param string: Store the string that is to be written into the file
    :param filename: Specify the name of the file to which data is saved
    :return: Nothing
    """
    with open(filename, "w") as f:
        f.write(string)
