from src.fasta import readFastaMul
from src.utils import dotPlot


def main():
    # print(readFastaMul("./data/sequences_synthetases.fasta"))
    dotPlot("./data/toaster.fasta")


if __name__ == "__main__":
    main()
