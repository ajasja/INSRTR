"""
Module to analyze the loops
"""


def get_loops_from_annotation(annot: str, min_length=2, skip_ends=True, loop_char="L"):
    """
    Returns a list of list of loop indices (looks for `loop_char` in the string)

    Parameters
    ----------
    annot : str
        String containing annotation of secondary structure
    min_length : int, optional
       Minimum length to pick up as loop. At least `min_length` residues must be present.
    skip_ends : bool, optional
        Do not include "loops" at the end, by default True
    loop_char : str, optional
        What does the loop char look like, by default "L"

    Returns
    -------
    list of list of list indices
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
        loops = [loop for loop in loops if len(loop) >= min_length]

    return loops


def loops_to_0_based(loops):
    """Decrements a list of lists. Returns a deep copy"""
    import copy

    # must make copy so results are nt modified in place
    res = copy.deepcopy(loops)

    for x in range(len(res)):
        for i in range(len(res[x])):
            res[x][i] = res[x][i] - 1
            assert res[x][i] > 0, "0 based indexing can not be less than 0"
    return res

import mdtraj as md
class LoopAnalyzer:
    def __init__(self, struct_file_path):
        self.struct_file_path = struct_file_path

        self.traj = md.load(struct_file_path)
        self.topology = self.traj.topology
        self.dssp = md.compute_dssp(self.traj, simplified=True)[0]
        self.loops = get_loops_from_annotation(self.dssp, loop_char='C', skip_ends=True)
        self.loops0 = loops_to_0_based(self.loops)
        self.full_sasa =  md.shrake_rupley(self.traj)[0]

    _loop_features = []
    _resi_features = []



    def analyze_structure(self):
        """Analyze the structure"""

        self.get_loop_features()

        self.get_resi_features()

    def get_loop_features(self):
        self._loop_features = []
        for li, loop  in enumerate(self.loops0):
            self._loop_features.append({}) # make new dict
            for loop_analyzer in self._loop_analyzers:
                res = loop_analyzer(self, li)
                self._loop_features[li].update(res)



    def get_resi_features(self):
        pass

    def get_loop_sasa(self, loop_index0):
        print(f'Hi {loop_index0}')
        return(dict(test=1))


    _loop_analyzers = [get_loop_sasa]
