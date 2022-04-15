from dataclasses import dataclass
import itertools
from enum import Enum
from typing import List, Tuple
from .BLOSUM import BLOSUM


class Direction(Enum):
    UP = "up"
    LEFT = "left"
    DIAG = "diag"

    def __repr__(self):
        return self.value

@dataclass
class Coord(object):
    """Class coordinates for position in a matrix"""
    x: int
    y: int

    def split(self):
        """
        The split function takes a string and returns a list of strings,
        split on the specified delimiter. If no delimiter is specified,
        whitespace is used.
        
        :param self: Reference the class itself
        :return: A list of strings
        """
        return self.x, self.y

    def __repr__(self):
        return f"({self.x}, {self.y})"


# # Type of alignment
# Coord = Tuple[int, int]

# @TODO: add docstring and tests
class Alignment:
    """
    Alignment class to perform global or local alignment of two sequences
    """
    def __init__(self, seq1, seq2, match=2, mismatch=-1, gap=-1):
        """
        The __init__ function initializes the class with two sequences, a match score,
        a mismatch score, and a gap penalty. It also initializes matScores to be an empty
        list of lists (the same size as seqs) and matDir to be an empty list of lists (also
        the same size as seqs). The for loops in __init__ fill these two variables with 0's.
        
        :param self: Reference the object itself
        :param seq1: Store the first sequence
        :param seq2: Store the second sequence
        :param match=2: Set the match score
        :param mismatch=-1: Set the score for mismatches
        :param gap=-1: Set the gap penalty
        :return: The list of lists matscores, which is the matrix that contains the scores for each pair of letters in seqs
        """
        self.seqs = (seq1, seq2)
        self.bestScore = None
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

    def __resetMatrix(self)-> None:
        """
        The __resetMatrix function is a helper function that resets the matScores and matDir
        matrices to all zeros. This is done so that when we call alignement function, 
        we can reset these two matrices to their original state before we start our alignment.
        
        :param self: Represent the instance of the object itself
        :return: A matrix of lists with the same dimensions as the input sequences, but instead of being filled with numbers, it is filled with empty lists
        """
        
        self.matScores = [
            [0 for _ in range(len(self.seqs[0]) + 1)]
            for _ in range(len(self.seqs[1]) + 1)
        ]
        self.matDir = [
            [[] for _ in range(len(self.seqs[0]) + 1)]
            for _ in range(len(self.seqs[1]) + 1)
        ]

    def __calculateScore(self, coord: Coord, previousScore: int, match: bool, useBlosum: bool) -> Tuple[int, Direction]:
        """
        The __calculateScore function calculates the score of a given coordinate in the matrix.
        The function takes as input:
            coord - The coordinate to be calculated, which is an object of type Coord.
            previousScore - The score at the previous position in the matrix, which is an integer. 
            match - A boolean value indicating whether or not there was a match between two letters at this position in alignment (i.e., did they both come from self.seqs[0][coord] and self.seqs[0][coord]).  This value will be True if there was a match and False otherwise (i..e, if they were mismatched).  
            useBlosum: A boolean value indicating whether or not we should use BLOSUM62 scoring instead of simple scoring.
        
        :param self: Reference the class object itself
        :param coord:Coord: Store the x and y coordinates of the current cell
        :param previousScore: Store the score of the cell that is above and to the left of current cell
        :param match: Determine whether the current cell is a match or not
        :param useBlosum:bool: Determine whether the score should be calculated using the blosum matrix or not
        :return: The score of the current cell and the direction to get there
        """
        if useBlosum:
            return (
                previousScore
                + BLOSUM[self.seqs[0][coord.y - 1]][self.seqs[1][coord.x - 1]],
                Direction.DIAG,
            )
        else:
            return (
                previousScore + (self.match if match else self.mismatch),
                Direction.DIAG,
            )

    def __bestAction(
        self, coord: Coord, useBlosum: bool, use0: bool = False
    ) -> List[Tuple[int, Direction]]:
        """
        The __bestAction function returns the best action with coord, maximal score and direction.
        The function takes as input:
        -coord: a coordinate in the matrix of scores (x,y)
        -useBlosum: boolean to use blosum62 or not. If True it will use it otherwise not. 
        -use0 : boolean to add 0 or not in the calculation of maxScore
        
        :param self: Reference the class itself
        :param coord:Coord: Specify the current position in the matrix
        :param useBlosum:bool: Determine whether the blosum matrix should be used or not
        :param use0:bool=False: Determine whether or not to use the 0 score
        :return: The best action with coord, maximal score and direction
        """

        previousScore = self.matScores[coord.x - 1][coord.y - 1]
        score = [
            self.__calculateScore(
                coord,
                previousScore,
                self.seqs[1][coord.x - 1] == self.seqs[0][coord.y - 1],
                useBlosum,
            ),
        ]

        previousScore = self.matScores[coord.x - 1][coord.y]
        score.append((previousScore + self.gap, Direction.UP))

        previousScore = self.matScores[coord.x][coord.y - 1]
        score.append((previousScore + self.gap, Direction.LEFT))

        if use0:
            previousScore = self.matScores[coord.x - 1][coord.y - 1]
            score.append((0, None))

        maxScore = max(score, key=lambda x: x[0])[0]
        return [i for i in score if i[0] == maxScore]

    def NWSIterFill(self, useBlosum=False):
        """
        The NWSIterFill function fills the matScores and matDir arrays with the
        alignment scores for each cell in the matrix. The function iterates through
        each cell of the matrix, starting at (0, 0) and ending at (len(seq2), len(seq2)).
        The function fills each cell by calculating all possible alignment scores for that 
        cell based on its neighbors' values. The score is calculated using a linear gap penalty 
        and either a blosum62 substitution matrix or an affine gap penalty. If useBlosum is True, then 
        the blosum62 substitution matrix will be used to calculate alignment scores; otherwise it will use an affine gap penalty.
        
        :param self: Access the class variables, such as self
        :param useBlosum=False: Indicate whether the blosum matrix should be used or not
        :return: The best score for the given coordinate
        """
        self.__resetMatrix()
        for i, j in itertools.product(
            range(len(self.seqs[1]) + 1), range(len(self.seqs[0]) + 1)
        ):
            if j == 0 or i == 0:
                self.matScores[i][j] = self.gap * max(i, j)
                if j == 0:
                    self.matDir[i][j].append(Direction.UP)
                else:
                    self.matDir[i][j].append(Direction.LEFT)
            else:
                for score, direction in self.__bestAction(Coord(i, j), useBlosum):
                    self.matScores[i][j] = score
                    self.matDir[i][j].append(direction)

        self.bestScore = (self.matScores[-1][-1], Coord(len(self.seqs[1]), len(self.seqs[0])))

    def NWSBacktrack(self):
        """
        The NWSBacktrack function takes the matrix of scores and directions as input,
        and returns a list of aligned sequences. The function starts from the end of the
        matrix (the bottom right corner), and works its way back to the origin. At each step,
        it adds an '-' to one sequence if it is moving up in the matrix (i.e., going from 
        columns to rows), and adds a letter from another sequence if it is moving left in 
        the matrix (going from rows to columns). This process continues until we reach either 
        the top or left edge.
        
        :param self: Access the class attributes
        :return: The alignment of the two sequences in a list
        """
        # we start from the end of the matrix
        j, i = map(len, self.seqs)
        while i > 0 or j > 0:
            match self.matDir[i][j][-1]:
                case Direction.UP:
                    self.aliSeqs[0] = "-" + self.aliSeqs[0]
                    self.aliSeqs[1] = self.seqs[1][i - 1] + self.aliSeqs[1]
                    i -= 1
                case Direction.LEFT:
                    self.aliSeqs[0] = self.seqs[0][j - 1] + self.aliSeqs[0]
                    self.aliSeqs[1] = "-" + self.aliSeqs[1]
                    j -= 1
                case Direction.DIAG:
                    self.aliSeqs[0] = self.seqs[0][j - 1] + self.aliSeqs[0]
                    self.aliSeqs[1] = self.seqs[1][i - 1] + self.aliSeqs[1]
                    i -= 1
                    j -= 1

    
    def SWIter(self, useBlosum=False):
        """
        The SWIter function is a generator that yields the best alignment of two sequences.
        It does so by iteratively applying the Smith-Waterman algorithm to align two sequences.
        The function takes one argument, useBlosum, which is a boolean value indicating whether or not to use BLOSUM62 as the substitution matrix.
        
        :param self: Access the class's attributes and methods
        :param useBlosum=False: Indicate whether the blosum62 matrix should be used or not
        :return: The best score for the current cell
        """
        self.__resetMatrix()
        self.bestScore = None # we need to reset the best score
        for i, j in itertools.product(
            range(len(self.seqs[1]) + 1), range(len(self.seqs[0]) + 1)
        ):
            if j == 0 or i == 0:
                self.matScores[i][j] = max(0, self.gap) * max(i, j) # initialisation we use 0 to reset score
                if j == 0:
                    self.matDir[i][j].append(Direction.UP if self.gap > 0 else None)
                else:
                    self.matDir[i][j].append(Direction.LEFT if self.gap > 0 else None)
            else:
                for score, direction in self.__bestAction(Coord(i, j), useBlosum=useBlosum, use0=True):
                    self.matScores[i][j] = score
                    self.matDir[i][j].append(direction)
                    if not self.bestScore or score > self.bestScore[0]:
                        self.bestScore = (score, Coord(i, j))

    def SWBacktrack(self):
        """
        The SWBacktrack function takes a list of sequences and the scoring matrix as input.
        It then uses the Smith-Waterman algorithm to find the best local alignment between all pairs of sequences in seqs.
        The function returns a list of aligned sequence strings.
        
        :param self: Access the attributes and methods of the class in python
        :return: The alignment of the two sequences
        """
        # we start from the end of the matrix and go to the best score with only gap
        self.__goToCoord(Coord(len(self.seqs[0]), len(self.seqs[1])), self.bestScore[1])

        i, j = self.bestScore[1].x, self.bestScore[1].y

        while i > 0 or j > 0:
            match self.matDir[i][j][-1]:
                case Direction.UP:
                    self.aliSeqs[0] = "-" + self.aliSeqs[0]
                    self.aliSeqs[1] = self.seqs[1][i - 1] + self.aliSeqs[1]
                    i -= 1
                case Direction.LEFT:
                    self.aliSeqs[0] = self.seqs[0][j - 1] + self.aliSeqs[0]
                    self.aliSeqs[1] = "-" + self.aliSeqs[1]
                    j -= 1
                case Direction.DIAG:
                    self.aliSeqs[0] = self.seqs[0][j - 1] + self.aliSeqs[0]
                    self.aliSeqs[1] = self.seqs[1][i - 1] + self.aliSeqs[1]
                    i -= 1
                    j -= 1
                case None:
                    break
        
        # now we have the best local alignment we go to start of the matrix
        self.__goToCoord(Coord(j, i), Coord(0,0))

    def __goToCoord(self, coord: Coord, target: Coord):
        """
        The __goToCoord function takes in a coordinate and a target coordinate. It then moves the sequence that is currently being aligned to the target coordinate, by moving up or down rows until it reaches the target.
        
        
        :param self: Access the attributes and methods of the class in python
        :param coord:Coord: Store the coordinates of the current position in the sequence
        :param target:Coord: Specify the target coordinate
        :return: The alignment of the two sequences from the current coordinate to target coordinate
        """
        i, j = coord.x, coord.y
        for k in range(i, target.y, -1):
            self.aliSeqs[0] = self.seqs[0][k - 1] + self.aliSeqs[0]
            self.aliSeqs[1] = "-" + self.aliSeqs[1]
        for k in range(j, target.x, -1):
            self.aliSeqs[0] = "-" + self.aliSeqs[0]
            self.aliSeqs[1] = self.seqs[1][k - 1] + self.aliSeqs[1]
        
        return target


    def __repr__(self):
        return str(self.aliSeqs[0] + "\n" + self.aliSeqs[1])
