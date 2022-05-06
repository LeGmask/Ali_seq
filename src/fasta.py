import re
from typing import Tuple

from .sequence import Sequence


def readFasta(file: str) -> Tuple[str, Sequence]:
    """
    The readFasta function reads a fasta file and returns the sequence in it.
    It takes one argument, which is the name of the fasta file to be read.
    The function reads this file line by line, starting with the first line
    (which should have a '>' at its start), and ending with an empty string.

    :param file:str: Specify the file to read from
    :return: A tuple of two values, the first is the sequence as a string, the second is a sequence object
    """
    with open(file, "r") as f:
        content = "".join(f.readlines())
        m = re.search(">(.*)\n([\s\S]+?(?=>|$))", content)  # find the first fasta input
        return m[2].replace("\n", ""), Sequence(m[1], m[2].replace("\n", ""))


def readFastaMul(file: str) -> Tuple[Sequence]:
    """
    The readFastaMul function reads a FASTA file and returns a tuple of Sequence objects.
    The function takes one argument, the name of the FASTA file to be read.

    :param file: Specify the file that is to be read
    :return: A tuple of sequence objects
    """
    with open(file, "r") as f:
        content = "".join(f.readlines())
        return tuple(
            Sequence(m.group(1), m.group(2).replace("\n", ""))
            for m in re.finditer(">(.*)\n([\s\S]+?(?=>|$))", content)
        )
