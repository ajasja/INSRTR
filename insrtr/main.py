"""Main module."""


def insert_sequence(seq, one_based_position, insert):
    """Inserts an AA sequence into another AA sequence. Uses 1 based indexing"""
    pos = one_based_position - 1
    return  seq[:pos] + insert + seq[pos:]
