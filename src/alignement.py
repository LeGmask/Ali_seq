import itertools
from enum import Enum
from typing import List, Tuple


class Direction(Enum):
    DOWN = "down"
    RIGHT = "right"
    DIAG = "diag"

    def __repr__(self):
        return self.value


# Alias de type:
Coord = Tuple[int, int]


class Alignement:
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

    def bestAction(self, coord: Coord) -> List[Tuple[Coord, int, Direction]]:
        """ "
        Return the bestaction with coord, maximal score and direction
        """
        score = []
        for direction in Direction:
            match direction:
                case Direction.DIAG:
                    oldScore = self.matScores[coord[0] - 1][coord[1] - 1]
                    score.append(
                        (coord, oldScore + self.match, Direction.DIAG)
                        if self.seqs[1][coord[0] - 1] == self.seqs[0][coord[1] - 1]
                        else (coord, oldScore + self.mismatch, Direction.DIAG)
                    )
                case Direction.DOWN:
                    oldScore = self.matScores[coord[0] - 1][coord[1]]
                    score.append((coord, oldScore + self.gap, Direction.DOWN))
                case Direction.RIGHT:  # right
                    oldScore = self.matScores[coord[0]][coord[1] - 1]
                    score.append((coord, oldScore + self.gap, Direction.RIGHT))
        maxScore = max(score, key=lambda x: x[1])[1]
        return [i for i in score if i[1] == maxScore]

    def NWSIterFill(self):
        # Initialisation le gap est proportionel a l'index de la sequence
        for i in range(len(self.matScores)):
            self.matScores[i][0] = self.gap * i
            self.matDir[i][0].append(Direction.DOWN)
        for j in range(len(self.matScores[0])):
            self.matScores[0][j] = self.gap * j
            self.matDir[0][j].append(Direction.RIGHT)
        for j, i in itertools.product(
            range(1, len(self.seqs[0]) + 1), range(1, len(self.seqs[1]) + 1)
        ):
            for _, score, direction in self.bestAction((i, j)):
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
