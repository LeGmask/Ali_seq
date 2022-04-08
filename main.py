from src.alignement import Alignement
from src.fasta import readFastaMul
from src.utils import dotPlot, drawMatrix


def main():
    # print(readFastaMul("./data/sequences_synthetases.fasta"))
    dotPlot("./data/toaster.fasta")

    ali = Alignement("LAFLALMEE", "LAME", match=0, mismatch=-1, gap=-2)
    ali.NWSIterFill()

    drawMatrix(ali.matScores)
    drawMatrix(ali.matDir)
    ali.NWSBacktrack()
    print(ali)


if __name__ == "__main__":
    main()
