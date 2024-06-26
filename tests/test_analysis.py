import insrtr


def test_get_loops_from_annotation():

    # basic funct
    loops = insrtr.get_loops_from_annotation("HHHLLLLHHHH")
    assert loops == [[4, 5, 6, 7]]

    # multiple loops
    loops = insrtr.get_loops_from_annotation("HHHLLLHHHLLLLHHHHLLLLHHHH")
    assert loops == [[4, 5, 6], [10, 11, 12, 13], [18, 19, 20, 21]]

    # filtering
    loops = insrtr.get_loops_from_annotation(
        "HHHLLHHHLLLLHHHHLLHHHH", min_length=3)
    assert loops == [[9, 10, 11, 12]]

    # filtering
    loops = insrtr.get_loops_from_annotation(
        "LLLHHHLLHHHLLLLHHHHLLHHHHLLL", min_length=3)
    assert loops == [[12, 13, 14, 15]]

    # filtering at ends
    loops = insrtr.get_loops_from_annotation(
        "LLLHHHLLHHHLLLLHHHHLLHHHHLLL", min_length=3, skip_ends=True)
    assert loops == [[12, 13, 14, 15]]

    loops = insrtr.get_loops_from_annotation(
        "LLLHHHLLHHHLLLLHHHHLLHHHHLLL", min_length=0, skip_ends=True)
    assert loops == [[7, 8], [12, 13, 14, 15], [20, 21]]


def test_loops_to_0_based():
    loops = [[4, 5, 6], [10, 11, 12, 13], [18, 19, 20, 21]]
    loops0 = [[3, 4, 5], [9, 10, 11, 12], [17, 18, 19, 20]]
    assert insrtr.loops_to_0_based(loops) == loops0

    # TOOD test assertion if index goes below 0


if __name__ == "__main__":
    test_get_loops_from_annotation()
