from typing import List
from tabulate import tabulate


class Sequence:
    def __init__(self, identifier: str, sequence: str) -> None:
        self.id = identifier
        self.sequence = sequence

    def getIndex(self, index: int) -> str:
        return self.sequence[index]

    def getLength(self) -> int:
        return len(self.sequence)

    def getEmptySequence(self):
        return Sequence(self.id, "")

    def insertFromBacktrace(self, value: str) -> None:
        self.sequence = "".join([value, self.sequence])

    def __repr__(self) -> str:
        return f"{self.id} : {self.sequence}"


class GroupSequences:
    def __init__(self, sequences: List[Sequence]) -> None:
        """
        :param sequences: number of sequences in the group or list of sequences
        """
        self.sequences = sequences

    def getIndex(self, index: int) -> str:
        return [sequence.getIndex(index) for sequence in self.sequences]

    def resetForAlignment(self) -> None:
        """
        The resetForAlignment function reset the sequences to insert data obtained
        from alignement, and move the previous sequence in the beforeAlignment
        attribute to be able to check wich sequence is the best one for future
        alignement.

        :param self: Access variables that belongs to the class
        :return: None
        """
        self.beforeAlignment = self.sequences
        self.sequences = [sequence.getEmptySequence() for sequence in self.sequences]

    def getLength(self) -> int:
        return len(self.sequences)

    def insertFromBacktrace(self, value: str) -> None:
        """
        The insertFromBacktrace function insert a value in the sequences starting
        from the end.

        :param self: Access variables that belongs to the class
        :param value: value to insert
        :return: None
        """
        for index, sequence in enumerate(self.sequences):
            sequence.insertFromBacktrace(value[index])

    def __repr__(self) -> str:
        return f"{self.sequences}"
