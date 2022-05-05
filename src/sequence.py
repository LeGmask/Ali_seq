from textwrap import shorten
from typing import List


class Sequence:
    def __init__(self, identifier: str, sequence: str) -> None:
        self.id = identifier
        self.sequence = sequence
        self.numberOfSequences = 1

    def getIndex(self, index: int) -> str:
        """
        The getIndex function returns the index of a given value in the list.

        :param self: Access the attributes and methods of the class
        :param index:int: Specify the index of the element in the list to be returned
        :return: The value of the indexth element in the list
        """
        return self.sequence[index]

    def getLength(self) -> int:
        return len(self.sequence)

    def getEmptySequence(self):
        """
        The getEmptySequence function returns the sequence with is sequence fields empty, keep only name.


        :param self: Access the attributes and methods of the class in python
        :return: The empty sequence
        """
        return Sequence(self.id, "")

    def insertFromBacktrace(self, value: str) -> None:
        """
        The insertFromBacktrace function inserts a value into the start of the sequence.

        :param self: Reference the object itself
        :param value: Value that is being inserted
        :return: None
        """
        self.sequence = "".join([value, self.sequence])

    def __repr__(self) -> str:
        """
        The __repr__ function is a special method that returns a string representation of the object.
        It is called by the repr() function, and is also used implicitly in other contexts where Python expects a string, such as tracebacks or inside of print statements.

        :param self: Refer to the instance of the class
        :return: A string that represents the object
        """
        return shorten(f"{self.id} : {self.sequence}", width=80, placeholder="...")


class GroupSequences:
    def __init__(self, sequences: List[Sequence]) -> None:
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
        """
        The getOriginalSequences function returns a list of the original sequences that was originally in the file.

        :param self: Refer to the object itself
        :return: A list of the original sequences
        """
        return self.originalSequences

    def getIndex(self, index: int, sequence: int = None) -> Sequence:
        """
        The getIndex function returns the index of a sequence in the alignment.
        The getIndex function takes two arguments:
            1) The index of the residu to be returned, and
            2) The index of the sequences were we get the residu (defaults to None).
        If no argument is provided for sequence, it will default to using
        the value stored as self.indexOfBestSequence. And if self.indexOfBestSequence is not set,
        it will return the last sequence (defaults for self.indexOfBestSequence).

        :param self: Refer to the object of the class
        :param index: Get the index of a specific sequence in the alignment
        :param sequence: Specify the sequence number
        :return: The sequence at the index
        """
        return self.beforeAlignment[sequence or self.indexOfBestSequence].getIndex(
            index
        )

    def getSequence(self, index: int) -> Sequence:
        """
        The getSequence function returns the sequence at the given index.

        :param self: Access the attributes and methods of the class
        :param index: Specify which sequence to return
        :return: The sequence at the specified index
        """
        return self.beforeAlignment[index]

    def getSequencesAtIndex(self, index: int) -> List[str]:
        """
        The getSequencesAtIndex function returns a list of residu that are at the given index for all sequences.



        :param self: Access the attributes and methods of the class in which it is used
        :param index:int: Specify the index of the residu to be returned
        :return: A list of residu at a given index for all the sequences
        """
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
        """
        The getGroupLength function returns the length of the group.

        :param self: Access the attributes and methods of the class in python
        :return: The length of the group
        """

        return len(self.sequences)

    def getLength(self, sequence: int = None) -> int:
        """
        The getLength function returns the length of a sequence.
        If no sequence is provided, it returns the length of the representing sequence.
        and if the representing sequence is not set, it returns the last sequence.

        :param self: Access the attributes and methods of the class in python
        :param sequence=None: Determine if the sequence is being passed in as a parameter or not
        :return: The length of the sequence
        """
        return self.beforeAlignment[sequence or self.indexOfBestSequence].getLength()

    def insertFromBacktrace(self, value: List[str]) -> None:
        """
        The insertFromBacktrace function insert a in all the sequences starting
        from the end.

        :param self: Access variables that belongs to the class
        :param value: List of values to be inserted
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
        return "".join((f"{names[i]}   {seqs[i]}" + "\n") for i in range(len(names)))
