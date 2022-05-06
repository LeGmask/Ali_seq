from ctypes import alignment
from itertools import combinations
from src.alignement import Alignment
from src.msa import Algorithm, Msa
from src.fasta import readFasta, readFastaMul
from src.utils import dotPlot, drawMatrix, saveToFile


def part_1_exercise_3():
    print(
        "-----------------------PART 1 — Exercise 1-------------------------------",
        end="\n\n",
    )
    # Here, function return at the same time the sequence as a plain string and a Sequence object.
    # Sequence object is a class from src/sequence.py used to store the sequence and its metadata.
    sequenceAsString, sequence = readFasta("data/HBA_HUMAN.fasta")
    print("1)", sequenceAsString)
    print()
    print("Optional, the sequence object:", sequence)

    print()
    sequences = readFastaMul("data/HBA_HUMAN-LGB2_LUPLU.fasta")
    print("2)", sequences)

    print()
    print("3)")
    dotPlot("data/dotPlot.fasta")
    print(
        "--------------------------------------------------------------------------",
        end="\n\n",
    )


def part_1_exercise_4():
    print(
        "-----------------------PART 1 — Exercise 4-------------------------------",
        end="\n\n",
    )
    alignment = Alignment("LAFLALMEE", "LAME", match=0, mismatch=-1, gap=-1)
    alignment.NWSIterFill()
    alignment.NWSBacktrack()

    print("For 'LAFLALMEE' and 'LAME', the algo give the following alignement :")
    print(alignment.aliSeqs)

    print("with the following matScores :")
    print(drawMatrix(alignment.matScores))
    print("and the following matDir :")
    print(drawMatrix(alignment.matDir))
    print(
        "--------------------------------------------------------------------------",
        end="\n\n",
    )


def part_2_exercise_5():
    print(
        "-----------------------PART 2 — Exercise 5-------------------------------",
        end="\n\n",
    )
    for i in combinations(["VSNS", "SNA", "AS"], 2):

        alignment = Alignment(i[0], i[1], match=0, mismatch=-1, gap=-1)
        alignment.NWSIterFill()
        alignment.NWSBacktrack()
        print("Pairwise alignment between", i[0], "and", i[1], "is :")
        print(alignment.aliSeqs)
    print(
        "--------------------------------------------------------------------------",
        end="\n\n",
    )


def part_3_exercise_6():
    print(
        "-----------------------PART 3 — Exercise 6-------------------------------",
        end="\n\n",
    )
    # Here we didn't rewrite the whole NWS funtion, instead we use an argument to use or not the blosum matrix.
    sequences = readFastaMul("data/HBA_HUMAN-LGB2_LUPLU.fasta")
    alignment = Alignment(sequences[0], sequences[1], gap=-7)
    alignment.NWSIterFill(useBlosum=True)
    alignment.NWSBacktrack()

    withBlosum = str(alignment.aliSeqs)
    print(
        "Pairwise alignment using BLOSUM62 matrix between",
        sequences[0],
        "and",
        sequences[1],
        "is :",
    )
    print(alignment)

    alignment.NWSIterFill(useBlosum=False)
    alignment.NWSBacktrack()
    withoutBlosum = str(alignment.aliSeqs)
    print(
        "Pairwise alignment without BLOSUM62 matrix between",
        sequences[0],
        "and",
        sequences[1],
        "is :",
    )
    print(alignment)
    print(
        "For a more readable output, the output has been save to a file, see 'results/part_3_exercise_6.result'"
    )
    saveToFile(
        "With BLOSUM62:\n" + withBlosum + "\nWithout BLOSUM62:\n" + withoutBlosum,
        "results/part_3_exercise_6.result",
    )
    print(
        "--------------------------------------------------------------------------",
        end="\n\n",
    )


def part_3_exercise_7():
    print(
        "-----------------------PART 3 — Exercise 7-------------------------------",
        end="\n\n",
    )
    # SWBacktrack function take as parameter full, which is a boolean to indicate if we want to use the full alignment or not.
    # If full is True, the function will return the full alignement, if not, the function will return the best local portion.

    sequences = readFastaMul("data/HBA_HUMAN-LGB2_LUPLU.fasta")
    alignment = Alignment(sequences[0], sequences[1], gap=-7)
    alignment.SWIter(useBlosum=True)
    alignment.SWBacktrack(full=True)
    SWresult = str(alignment.aliSeqs)
    print(
        "Pairwise alignment using SWIter and SWBacktrack between",
        sequences[0],
        "and",
        sequences[1],
        "is :",
    )
    print(alignment)

    alignment.NWSIterFill(useBlosum=True)
    alignment.NWSBacktrack()
    print(
        "Pairwise alignment using NWSIterFill and NWSBacktrack between",
        sequences[0],
        "and",
        sequences[1],
        "is :",
    )
    print(alignment)
    NWSresult = str(alignment.aliSeqs)

    print(
        "For a more readable output, the output has been save to a file, see 'results/part_3_exercise_7.result'"
    )
    saveToFile(
        "With SWIter:\n" + SWresult + "\nWith NWSIterFill:\n" + NWSresult,
        "results/part_3_exercise_7.result",
    )
    print(
        "--------------------------------------------------------------------------",
        end="\n\n",
    )


def part_3_exercise_8():
    print(
        "-----------------------PART 3 — Exercise 8-------------------------------",
        end="\n\n",
    )
    # here everything is bundled in the Msa class, which is a class from src/msa.py.
    # The Msa class is used to store the multiple sequence alignment and its metadata.
    # The Msa class has a function to align the sequences, which is the align function.
    # The Class take as parameter the algorithm to use, the different score and if we use the blosum matrix.
    sequences = readFastaMul("data/sequences_synthetases.fasta")
    msa = Msa(sequences, algorithm=Algorithm.NeedlemanWunsch, blosum=True)
    msa.align()
    print("MSA of the sequences is :")
    print(msa)

    print(
        "For a more readable output, the output has been save to a file, see 'results/part_3_exercise_8.result'"
    )
    saveToFile(str(msa), "results/part_3_exercise_8.result")
    print(
        "--------------------------------------------------------------------------",
        end="\n\n",
    )


if __name__ == "__main__":
    part_1_exercise_3()
    part_1_exercise_4()
    part_2_exercise_5()
    part_3_exercise_6()
    part_3_exercise_7()
    part_3_exercise_8()
