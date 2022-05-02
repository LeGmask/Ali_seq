from operator import methodcaller
from typing import List
from tabulate import tabulate

from src.sequence import GroupSequences, Sequence
from src.alignement import Alignment
from src.fasta import readFasta, readFastaMul
from src.utils import drawMatrix


def multAlign(seqs):
    multAlign = []
    scores = pairScores(seqs)
    m = max(scores, key=lambda x: scores.get(x))
    ali = Alignment(m[0].sequence, m[1].sequence, match=0, mismatch=-1, gap=-7)
    ali.NWSIterFill(useBlosum=True)
    ali.NWSBacktrack()
    multAlign.append(ali.aliSeqs[0])
    multAlign.append(ali.aliSeqs[1])

    print(multAlign)


def getScore(sesq1, sesq2):
    ali = Alignment(sesq1, sesq2, match=0, mismatch=-1, gap=-7)
    ali.NWSIterFill(useBlosum=True)
    ali.NWSBacktrack()
    return ali.bestScore[0]


def scoreMatrix(seqs: List[Sequence]):
    mat = [[None for _ in range(len(seqs))] for _ in range(len(seqs))]
    for i in range(len(seqs)):
        for j in range(i + 1, len(seqs)):
            ali = Alignment(seqs[i].sequence, seqs[j].sequence, gap=-7)
            ali.NWSIterFill(useBlosum=True)
            a = ali.bestScore[0]
            mat[i][j] = a
            mat[j][i] = a
    return mat
    # return {i: getScore(i[0].sequence, i[1].sequence) for i in combinations(seqs, 2)}


def main():
    # _, hba_Human = readFasta("data/HBA_HUMAN.fasta")
    # _, lgb2_luplu = readFasta("data/LGB2_LUPLU.fasta")

    ali = Alignment("VSNS", "SNA", match=1, mismatch=-1, gap=-7)
    ali.SWIter(useBlosum=False)
    ali.SWBacktrack(full=True)
    # ali.NWSIterFill(useBlosum=False)
    # ali.NWSBacktrack()
    drawMatrix(ali.matScores)
    print(ali.aliSeqs)

    sequences = readFastaMul("data/toaster.fasta")

    # gpe = GroupSequences(sequences)
    # gpe.resetForAlignment()
    # gpe.insertFromBacktrace(["S", "-", "S"])

    # print(gpe)

    # print(pairScores(sequences))
    # print(tabulate(scoreMatrix(sequences)))


if __name__ == "__main__":
    main()
