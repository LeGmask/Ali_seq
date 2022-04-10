import itertools
from enum import Enum
from typing import List, Tuple


class Direction(Enum):
    DOWN = "down"
    RIGHT = "right"
    DIAG = "diag"

    def __repr__(self):
        return self.value


# Type of alignment
Coord = Tuple[int, int]

# @TODO: add docstring and tests
class Alignment:
    def __init__(self, seq1, seq2, match=2, mismatch=-1, gap=-1):
        self.seqs = (seq1, seq2)
        self.match = match
        self.mismatch = mismatch
        self.gap = gap
        self.matScores = [
            [0 for _ in range(len(self.seqs[0]) + 1)]
            for _ in range(len(self.seqs[1]) + 1)
        ]
        self.matDir = [
            [[] for _ in range(len(self.seqs[0]) + 1)]
            for _ in range(len(self.seqs[1]) + 1)
        ]
        self.aliSeqs = ["", ""]

    def bestAction(self, coord: Coord) -> List[Tuple[int, Direction]]:
        """
        Return the bestaction with coord, maximal score and direction
        """
        previousScore = self.matScores[coord[0] - 1][coord[1] - 1]
        score = [
            (previousScore + self.match, Direction.DIAG)
            if self.seqs[1][coord[0] - 1] == self.seqs[0][coord[1] - 1]
            else (previousScore + self.mismatch, Direction.DIAG)
        ]

        previousScore = self.matScores[coord[0] - 1][coord[1]]
        score.append((previousScore + self.gap, Direction.DOWN))

        previousScore = self.matScores[coord[0]][coord[1] - 1]
        score.append((previousScore + self.gap, Direction.RIGHT))

        maxScore = max(score, key=lambda x: x[0])[0]
        return [i for i in score if i[0] == maxScore]

    def NWSIterFill(self):
        for j, i in itertools.product(
            range(len(self.seqs[0]) + 1), range(len(self.seqs[1]) + 1)
        ):
            if j == 0 or i == 0:
                self.matScores[i][j] = self.gap * max(i, j)
                if j == 0:
                    self.matDir[i][j].append(Direction.DOWN)
                else:
                    self.matDir[i][j].append(Direction.RIGHT)
            else:
                for score, direction in self.bestAction((i, j)):
                    self.matScores[i][j] = score
                    self.matDir[i][j].append(direction)

    def NWSBacktrack(self):
        # we start from the end of the matrix
        i, j = map(len, self.seqs)
        while i > 0 or j > 0:
            # Here we use pattern matching to get the direction (python 3.10)
            match self.matDir[j][i][-1]:
                case Direction.DOWN:
                    self.aliSeqs[0] = "-" + self.aliSeqs[0]
                    self.aliSeqs[1] = self.seqs[1][i - 1] + self.aliSeqs[1]
                    j, i = j - 1, i
                case Direction.RIGHT:
                    self.aliSeqs[0] = self.seqs[0][i - 1] + self.aliSeqs[0]
                    self.aliSeqs[1] = "-" + self.aliSeqs[1]
                    j, i = j, i - 1
                case Direction.DIAG:
                    self.aliSeqs[0] = self.seqs[0][i - 1] + self.aliSeqs[0]
                    self.aliSeqs[1] = self.seqs[1][j - 1] + self.aliSeqs[1]
                    j, i = j - 1, i - 1

    def __repr__(self):
        return str(self.aliSeqs[0] + "\n" + self.aliSeqs[1])
