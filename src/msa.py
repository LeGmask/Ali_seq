from enum import Enum
from itertools import combinations
from operator import methodcaller
from typing import List

from .alignement import Alignment
from .sequence import Sequence


class Algorithm(Enum):
    NeedlemanWunsch = 1
    SmithWaterman = 2


class Msa:
    def __init__(
        self,
        sequences: List[Sequence],
        algorithm: Algorithm,
        match: int = 1,
        mismatch: int = -1,
        gap: int = -7,
        blosum: bool = True,
    ) -> None:
        self.sequences = sequences
        self.scoreMatrix = None
        self.scoreDict = None
        # sequence tracking, for each sequence the current state, ie group it belongs to or sequence it belongs to
        self.seqManagement = {seq: seq for seq in sequences}

        self.algorithm = algorithm
        self.mismatch = mismatch
        self.blosum = blosum
        self.match = match
        self.gap = gap

        self.alignments = None

    def generateScoreMatrix(self) -> None:
        """
        The generateScoreMatrix function generates a score matrix for the sequences in the alignment.
        The score matrix is a dictionary of tuples, where each tuple contains two sequences and their pairwise alignment score.
        The function sorts this list by the scores (in descending order), so that it can be used to find an optimal solution for
        the sequence alignments.

        :param self: Access the attributes and methods of the class in python
        :return: None
        """
        self.scoreMatrix = {
            i: self.getPairwiseScore(i[0].sequence, i[1].sequence)
            for i in combinations(self.sequences, 2)
        }
        # Sort the score matrix by the score
        self.scoreMatrix = list(
            sorted(
                self.scoreMatrix.items(),
                key=methodcaller("__getitem__", 1),
                reverse=True,
            )
        )

    def generateScoreDict(self) -> None:
        """
        Generate the score dict in order to find the best representing sequence,
        the score of each key, aka each sequence, is the sum of pairwise score.

        :param self: Refer to the object itself
        :return: None
        """
        self.scoreDict = {}
        for i in self.scoreMatrix:
            for j in i[0]:
                self.scoreDict[j] = self.scoreDict.get(j, 0) + i[1]

    def getPairwiseScore(self, sesq1: Sequence, sesq2: Sequence) -> int:
        """
        The getPairwiseScore function takes two sequences and returns the number of matches between them.
        The function returns a score which is equal to the number of matches divided by length of sequence.

        :param self: Refer to the object of the class
        :param sesq1: First sequence
        :param sesq2: Second sequence
        :return: Score which is equal to the number of matches divided by length of sequence.
        """
        alignment = Alignment(
            sesq1, sesq2, match=self.match, mismatch=self.mismatch, gap=self.gap
        )

        if self.algorithm == Algorithm.NeedlemanWunsch:
            alignment.NWSIterFill(useBlosum=self.blosum)
            alignment.NWSBacktrack()
        elif self.algorithm == Algorithm.SmithWaterman:
            alignment.SWIter(useBlosum=self.blosum)
            alignment.SWBacktrack(full=True)
        else:
            raise ValueError("Algorithm not supported")

        # number of matches / length of the sequences
        return (
            sum(
                alignment.aliSeqs.getSequence(0).getIndex(k)
                == alignment.aliSeqs.getSequence(1).getIndex(k)
                for k in range(alignment.aliSeqs.getLength())
            )
            / alignment.aliSeqs.getLength()
        )

    def align(self) -> None:
        """
        Align the sequences by following the score matrix.

        :param self: Access the attributes and methods of the class
        :return: Nothing
        """
        # First, we need to generate the score matrix
        self.generateScoreMatrix()
        # Then, we need to generate the score dict in order to find the best representing sequence
        self.generateScoreDict()

        # We parse paires of sequences from best to worst
        for pair in self.scoreMatrix:
            # We check if the first both sequence are not already in the same group
            if self.seqManagement[pair[0][0]] != self.seqManagement[pair[0][1]]:
                alignment = Alignment(
                    self.seqManagement[pair[0][0]],
                    self.seqManagement[pair[0][1]],
                    match=self.match,
                    mismatch=self.mismatch,
                    gap=self.gap,
                )
                if self.algorithm == Algorithm.NeedlemanWunsch:
                    alignment.NWSIterFill(useBlosum=self.blosum)
                    alignment.NWSBacktrack()
                elif self.algorithm == Algorithm.SmithWaterman:
                    alignment.SWIter(useBlosum=self.blosum)
                    alignment.SWBacktrack(full=True)
                else:
                    raise ValueError("Algorithm not supported")

                self.alignment = alignment.aliSeqs

                # We set the best sequence of group in fucntion of sum of pairwise score (we take the best score)
                currectScoreDict = {
                    i: self.scoreDict[i]
                    for i in alignment.aliSeqs.getOriginalSequences()
                }
                alignment.aliSeqs.setBestSequence(
                    max(currectScoreDict, key=currectScoreDict.get)
                )

                # We update the seqManagement dict to keep track of current position of each sequence
                for key in alignment.aliSeqs.getOriginalSequences():
                    self.seqManagement[key] = alignment.aliSeqs

    def __repr__(self) -> str:
        return str(self.alignment)
