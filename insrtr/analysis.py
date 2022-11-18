"""
Module to analyze the loops
"""
from .utils import *

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
        self.loops = get_loops_from_annotation(self.dssp, loop_char="C", skip_ends=True)
        self.loops0 = loops_to_0_based(self.loops)
        self.sasa_A = md.shrake_rupley(self.traj)[0] * 100  # make in in angstrom
        self.total_sasa_A = sum(self.sasa_A)

    _loop_features = []
    _resi_features = []

    def analyze_structure(self):
        """Analyze the structure"""

        self.get_loop_features()

        self.get_resi_features()

    def get_loop_features(self):
        self._loop_features = []
        for li, loop in enumerate(self.loops0):
            self._loop_features.append(dict(loop_index0=li, loop_length_AA=len(loop)))  # make new dict
            for loop_analyzer in self._loop_analyzers:
                res = loop_analyzer(self, li)
                self._loop_features[li].update(res)

    def get_resi_features(self):
        pass

    def get_loop_geometry(self, loop_index0):
        loop_residues = self.loops0[loop_index0]
        loop_ids = self.topology.select(f"resid {loop_residues[0]} to {loop_residues[-1]}")
        loop_isolation_traj = self.traj.atom_slice(loop_ids)

        first_CA = self.topology.select(f"resid {loop_residues[0]} and name CA")[0]
        last_CA =  self.topology.select(f"resid {loop_residues[-1]} and name CA")[0]
        # returns a set of frames, but we only have one frame, so [0] is needed 
        loop_start_end_distance_A =  md.compute_distances(self.traj, [[first_CA, last_CA]] )[0][0]*10  

        loop_radius_gyration_A = md.compute_rg(loop_isolation_traj)[0]*10
        
        #TODO: calculate distance to active site

        return dict(
            loop_start_end_distance_A=loop_start_end_distance_A,
            loop_radius_gyration_A=loop_radius_gyration_A
        )

    def get_loop_sequence_features(self, loop_index0):
        """Returns sequence features, such as percent """
        loop_residues = self.loops0[loop_index0]
        #loop_ids = self.topology.select(f"resid {loop_residues[0]} to {loop_residues[-1]} and name CA")
        
        #get three letter names
        seq = resname_3to1([self.topology.residue(lid).name for lid in loop_residues])
        seq = "".join(seq)

        ll = len(seq)
        loop_G_percent = seq.count('G')/ll*100
        loop_P_percent = seq.count('P')/ll*100
        loop_S_percent = seq.count('S')/ll*100
        loop_T_percent = seq.count('T')/ll*100
        return dict(
            loop_seq=seq,
            loop_G_percent=loop_G_percent,
            loop_P_percent=loop_P_percent,
            loop_S_percent=loop_S_percent,
            loop_T_percent=loop_T_percent,
        )        


    def get_loop_sasa(self, loop_index0):
        """Returns loop sasa , loop sasa in isolation and relative loop sasa"""
        loop_residues = self.loops0[loop_index0]

        loop_ids = self.topology.select(f"resid {loop_residues[0]} to {loop_residues[-1]}")

        loop_sasa_A = sum(self.sasa_A[loop_ids])  # get sasa just for loop
        loop_sasa_A_per_res = loop_sasa_A / len(loop_residues)

        # get SASA if the loop was on it's own, without the rest of the protein
        loop_isolation_traj = self.traj.atom_slice(loop_ids)
        loop_isolation_SASA_A = sum(md.shrake_rupley(loop_isolation_traj)[0] * 100)  # make in in angstrom
        loop_burial_percent = loop_sasa_A/loop_isolation_SASA_A*100
        loop_percent_of_total_surface = loop_sasa_A/self.total_sasa_A*100
       
        return dict(
            loop_sasa_A=loop_sasa_A,
            loop_sasa_A_per_res=loop_sasa_A_per_res,
            loop_isolation_SASA_A=loop_isolation_SASA_A,
            loop_burial_percent=loop_burial_percent,
            loop_percent_of_total_surface=loop_percent_of_total_surface
        )

    _loop_analyzers = [get_loop_geometry, get_loop_sasa, get_loop_sequence_features]
