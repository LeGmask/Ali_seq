from src.alignement import Alignment
from src.fasta import readFasta, readFastaMul
from src.utils import dotPlot, drawMatrix, saveMatrix


def main():
    _, hba_Human = readFasta("data/HBA_HUMAN.fasta")
    _, lgb2_luplu = readFasta("data/LGB2_LUPLU.fasta")
    dotPlot("./data/toaster.fasta")

    ali = Alignment(hba_Human.sequence, lgb2_luplu.sequence,  match=0, mismatch=-1, gap=-7)
    ali.NWSIterFill(useBlosum=True)
    ali.NWSBacktrack()

    saveMatrix(ali.matScores, "./data/ali.fasta.matDir")
    print(ali)


if __name__ == "__main__":
    main()
