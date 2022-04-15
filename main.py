from src.alignement import Alignment
from src.fasta import readFasta, readFastaMul
from src.utils import dotPlot, drawMatrix, saveMatrix
from itertools import combinations


def pairScores(seqs):
    scores = {}
    for i in combinations(seqs, 2):
        ali = Alignment(i[0].sequence, i[1].sequence, match=0, mismatch=-1, gap=-7)
        ali.NWSIterFill(useBlosum=True)
        scores[i] = ali.bestScore[0]
    return scores


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


def main():
    _, hba_Human = readFasta("data/HBA_HUMAN.fasta")
    _, lgb2_luplu = readFasta("data/LGB2_LUPLU.fasta")

    ali = Alignment("LALMEE", "LAME", match=1, mismatch=-1, gap=-2)
    ali.SWIter(useBlosum=True)
    ali.SWBacktrack()
    # ali.NWSIterFill(useBlosum=True)
    # ali.NWSBacktrack()
    print(ali.aliSeqs)

    # sequences = readFastaMul("data/toaster.fasta")

    # # print(pairScores(sequences))
    # multAlign(sequences)


if __name__ == "__main__":
    main()
