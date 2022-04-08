class Sequence:
    def __init__(self, identifier: str, sequence: str) -> None:
        self.id = identifier
        self.sequence = sequence

    def __repr__(self) -> str:
        return f"{self.id} : {self.sequence}"
