from textwrap import shorten
from typing import List
from tabulate import tabulate


class Sequence:
    def __init__(self, identifier: str, sequence: str) -> None:
        self.id = identifier
        self.sequence = sequence
        self.numberOfSequences = 1

    def getIndex(self, index: int) -> str:
        return self.sequence[index]

    def getLength(self) -> int:
        return len(self.sequence)

    def getEmptySequence(self):
        return Sequence(self.id, "")

    def insertFromBacktrace(self, value: str) -> None:
        self.sequence = "".join([value, self.sequence])

    def __repr__(self) -> str:
        return shorten(f"{self.id} : {self.sequence}", width=80, placeholder="...")


class GroupSequences:
    def __init__(self, sequences: List[Sequence]) -> None:
        """
        :param sequences: number of sequences in the group or list of sequences
        """
        self.originalSequences: List[Sequence] = []
        self.sequences: List[Sequence] = []

        for i in sequences:
            if isinstance(i, GroupSequences):
                self.originalSequences = self.originalSequences + i.originalSequences
                self.sequences = self.sequences + i.sequences
            else:
                self.originalSequences.append(i)
                self.sequences.append(i)

        self.numberOfSequences = len(self.sequences)
        self.indexOfBestSequence = -1
        
        self.resetForAlignment()

    def setBestSequence(self, index: Sequence) -> None:
        """
        The setBestSequence function set the index of the best sequence.
        
        :param index: index of the best sequence
        :return: None
        """
        self.indexOfBestSequence = self.originalSequences.index(index)


    def getOriginalSequences(self) -> List[Sequence]:
        return self.originalSequences            

    def getIndex(self, index: int, sequence: int = None) -> Sequence:
        return self.beforeAlignment[sequence or self.indexOfBestSequence].getIndex(index)

    def getSequence(self, index: int) -> Sequence:
        return self.beforeAlignment[index]

    def getSequencesAtIndex(self, index: int) -> List[str]:
        return [sequence.getIndex(index) for sequence in self.beforeAlignment]

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

    def setAsAligned(self) -> None:
        """
        The setAsAligned function set beforeAlignment the to be the same as the
        sequences attribute.

        :param self: Access variables that belongs to the class
        :return: None
        """
        self.beforeAlignment = self.sequences

    def getGroupLength(self) -> int:
        return len(self.sequences)
    
    def getLength(self, sequence = None) -> int:
        return self.beforeAlignment[sequence or self.indexOfBestSequence].getLength()

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
        names, seqs = [], []
        for sequence in self.sequences:
            names.append(sequence.id)
            seqs.append(sequence.sequence)

        # get max lenght in names
        max_len = max(len(name) for name in names)
        # add spaces to names
        names = [name + " " * (max_len - len(name)) for name in names]
        return ''.join((f"{names[i]}   {seqs[i]}" + "\n") for i in range(len(names)))
