"""
Module to analyze the loops
"""


def get_loops_from_annotation(annot: str, min_length=2, skip_ends=True, loop_char="L"):
    """
    Returns a list of loop indices (looks for LLL in the string)

    Parameters
    ----------
    annot : str
        Annotation containg L
    min_length : int, optional
        _description_, by default 2
    skip_ends : bool, optional
        _description_, by default True
    """

    in_loop = False
    loops = []

    for n, c in enumerate(annot, start=1):
        if c == loop_char:
            if in_loop:
                single_loop.append(n)
            else:  # starting new loop
                single_loop = []
                single_loop.append(n)
                in_loop = True
        else:
            if in_loop:  # end loop and push loop to loops
                loops.append(single_loop)
            in_loop = False
    else:  # end of loop don't forget to add last loop
        if in_loop:
            loops.append(single_loop)

    # filter "loops" at the begging and end
    if skip_ends:
        # If string ends with loop char, it's just some ....LLLL at the end and not the loop
        if annot[-1] == loop_char:
            del loops[-1]

        # same for beggining
        if annot[0] == loop_char:
            del loops[0]

    # filter on min size
    if min_length > 0:
        loops = [loop for loop in loops if len(loop) > min_length]

    return loops
