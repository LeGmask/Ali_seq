from ctypes import alignment
from itertools import combinations
from operator import methodcaller
from typing import List, Tuple
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
    multAlign.append(ali.aliSeqs    [1])

    print(multAlign)


def getScore(sesq1, sesq2):
    ali = Alignment(sesq1, sesq2, match=0, mismatch=-1, gap=-7)
    ali.NWSIterFill(useBlosum=True)
    ali.NWSBacktrack()
    
    # number of matches / length of the sequences
    return sum(ali.aliSeqs.getSequence(0).getIndex(k) == ali.aliSeqs.getSequence(1).getIndex(k) for k in range(ali.aliSeqs.getLength())) / ali.aliSeqs.getLength()


def scoreMatrix(seqs: List[Sequence]):
    # mat = [[None for _ in range(len(seqs))] for _ in range(len(seqs))]
    # for i in range(len(seqs)):
    #     for j in range(i + 1, len(seqs)):
    #         # on aligne les deux séquence,

    #         ali = Alignment(seqs[i].sequence, seqs[j].sequence, gap=-7)
    #         ali.SWIter(useBlosum=True)
    #         ali.SWBacktrack(full=True)
    #         print(ali.aliSeqs)
    #         # number of matches / length of the sequences
    #         score = sum(ali.aliSeqs.getSequence(0).getIndex(k) == ali.aliSeqs.getSequence(1).getIndex(k) for k in range(ali.aliSeqs.getLength())) / ali.aliSeqs.getLength()


    #         mat[i][j] = score
    #         mat[j][i] = score
    # return mat
    
    return {i: getScore(i[0].sequence, i[1].sequence) for i in combinations(seqs, 2)}

def saveToFile(string, filename):
    with open(filename, "w") as f:
        f.write(string)

def main():
    # _, hba_Human = readFasta("data/HBA_HUMAN.fasta")
    # _, lgb2_luplu = readFasta("data/LGB2_LUPLU.fasta")

    # ali = Alignment(Sequence("toaster","VSN-S"), "SNA-", match=0, mismatch=-1, gap=-7)
    # ali.NWSIterFill(useBlosum=True)
    # drawMatrix(ali.matScores)
    # ali.NWSBacktrack()
    # ali.NWSIterFill(useBlosum=False)
    # ali.NWSBacktrack()
    # print(ali.aliSeqs)

    sequences = readFastaMul("data/sequences_synthetases.fasta")


    dictionary = scoreMatrix(sequences)
    dictionary = list(sorted(dictionary.items(), key=methodcaller('__getitem__', 1), reverse=True))
    print(tabulate(dictionary, headers=["Sequences", "Score"]))
    # for i in dictionary:


    seqManagement = {seq : seq for seq in sequences}
    scoreDict = {}
    for i in dictionary:
        for j in i[0]:
            scoreDict[j] = scoreDict.get(j, 0) + i[1]
    alignment = None

    for i in dictionary:
        if seqManagement[i[0][0]] != seqManagement[i[0][1]]:
            ali = Alignment(seqManagement[i[0][0]], seqManagement[i[0][1]], match=1, mismatch=-1, gap=-7)
            ali.NWSIterFill(useBlosum=True)
            ali.NWSBacktrack()

            alignment = ali.aliSeqs

            currectScoreDict = {i : scoreDict[i] for i in ali.aliSeqs.getOriginalSequences()}
            ali.aliSeqs.setBestSequence(max(currectScoreDict, key=currectScoreDict.get))


            for key in ali.aliSeqs.getOriginalSequences():
                seqManagement[key] = ali.aliSeqs


    saveToFile(str(alignment), "data/test.fasta.result")

    # ali = Alignment(dictionary[0][0][0], dictionary[0][0][1], match=1, mismatch=-1, gap=-7)
    # ali.SWIter(useBlosum=False)
    # ali.SWBacktrack(full=True)
    # group = ali.aliSeqs
    # print(group)
    # ali = Alignment(group, dictionary[1][0][1], match=1, mismatch=-1, gap=-7)
    # ali.SWIter(useBlosum=False)
    # ali.SWBacktrack(full=True)
    # print(ali.aliSeqs)


    # dictionary[0]
    



    # gpe = GroupSequences(sequences)
    # gpe.resetForAlignment()
    # gpe.insertFromBacktrace(["S", "-", "S"])

    # print(gpe)

    # print(pairScores(sequences))
    # print(tabulate(scoreMatrix(sequences)))


if __name__ == "__main__":
    main()
