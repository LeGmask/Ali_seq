from src.msa import Algorithm, Msa
from src.fasta import readFastaMul
from src.utils import saveToFile


def main():

    sequences = readFastaMul("data/sequences_synthetases.fasta")
    msa = Msa(sequences, algorithm=Algorithm.NeedlemanWunsch, blosum=True)
    msa.align()

    saveToFile(str(msa), "data/test.fasta.result")


if __name__ == "__main__":
    main()
