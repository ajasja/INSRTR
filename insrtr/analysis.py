"""
Module to analyze the loops
"""
from .utils import *
import numpy as np
import pathlib
from collections import Counter
import pandas as pd
import itertools


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
import pandas as pd


class LoopAnalyzer:
    def __init__(
        self,
        struct_file_path,
        struct_name=None,
        active_res_index1=None,
        always_include_sites1=None,
        include_dssp="LHE",
        skip_ends=True,
    ):
        self.struct_file_path = struct_file_path

        if struct_name is None:  # If no name given take it from the struct file
            struct_name = pathlib.Path(struct_file_path).stem

        if always_include_sites1 is None:
            always_include_sites1 = []

        self.struct_name = str(struct_name)
        self.traj = md.load(struct_file_path)
        self.topology = self.traj.topology
        self.dssp = md.compute_dssp(self.traj, simplified=True)[0]
        self.dssp = np.char.replace(self.dssp, "C", "L")
        self.seq = "".join(resname_3to1([res.name for res in self.topology.residues]))
        # self.loops = get_loops_from_annotation(self.dssp, loop_char="L", skip_ends=True) + always_include_sites1
        loop_array = [get_loops_from_annotation(self.dssp, loop_char=dssp_char, skip_ends=skip_ends) for dssp_char in include_dssp]
        self.loops = list(itertools.chain(*loop_array))
        # print(self.loops)
        self.loops0 = loops_to_0_based(self.loops)
        self.sasa_atoms_A = md.shrake_rupley(self.traj)[0] * 100  # make in in angstrom
        self.total_sasa_A = sum(self.sasa_atoms_A)
        self.residue_features_table = None
        if active_res_index1:
            self.active_res_index0 = [int(resid) - 1 for resid in active_res_index1]
        else:
            self.active_res_index0 = []

    _loop_features = []
    loop_feature_descriptions = {}
    _resi_features = []
    resi_feature_descriptions = {}
    feature_descriptions_table = None

    def analyze_structure(self):
        """Analyze the structure"""

        self.get_loop_features()

        self.get_resi_features()

        resi_table = pd.DataFrame(self._resi_features)
        loop_table = pd.DataFrame(self._loop_features)
        self.residue_features_table = resi_table.merge(loop_table, on="loop_index0", how="left")

        # Merge the descriptions as well
        self.feature_descriptions_table = self.resi_feature_descriptions.copy()
        self.feature_descriptions_table.update(self.loop_feature_descriptions)
        self.feature_descriptions_table = pd.DataFrame(self.feature_descriptions_table, index=["description"])
        self.feature_descriptions_table = self.feature_descriptions_table.T

        return self.residue_features_table

    def get_loop_features(self):
        self._loop_features = []
        for li, loop in enumerate(self.loops0):
            self._loop_features.append(dict(loop_index0=li, loop_length_AA=len(loop)))  # make new dict
            if li == 0:  # description need to be added on first pass only
                self.loop_feature_descriptions["loop_index0"] = "The zero based index of the loop."
                self.loop_feature_descriptions["loop_length_AA"] = "The length of the loop."
            for loop_analyzer in self._loop_analyzers:
                res = loop_analyzer(self, li)
                for f in res.keys():
                    # add the values and descriptions to different lists
                    self._loop_features[li][f] = res[f][0]
                    if li == 0:  # add the description if this is the first loop pass
                        self.loop_feature_descriptions[f] = res[f][1]

    def get_loop_geometry(self, loop_index0):
        loop_residues = self.loops0[loop_index0]
        loop_ids = self.topology.select(f"resid {loop_residues[0]} to {loop_residues[-1]}")
        loop_isolation_traj = self.traj.atom_slice(loop_ids)

        first_CA = self.topology.select(f"resid {loop_residues[0]} and name CA")[0]
        last_CA = self.topology.select(f"resid {loop_residues[-1]} and name CA")[0]
        # returns a set of frames, but we only have one frame, so [0] is needed
        loop_start_end_distance_A = md.compute_distances(self.traj, [[first_CA, last_CA]])[0][0] * 10

        loop_radius_gyration_A = md.compute_rg(loop_isolation_traj)[0] * 10

        # TODO: calculate distance to active site

        return dict(
            loop_start_end_distance_A=(
                loop_start_end_distance_A,
                "Distance between the CA atom of the first loop residue and CA atom of last loop residue.",
            ),
            loop_radius_gyration_A=(loop_radius_gyration_A, "Radius of Gyration of the loop residues."),
        )

    def get_loop_sequence_features(self, loop_index0):
        """Returns sequence features, such as percent"""
        loop_residues = self.loops0[loop_index0]
        # loop_ids = self.topology.select(f"resid {loop_residues[0]} to {loop_residues[-1]} and name CA")

        # get three letter names
        seq = resname_3to1([self.topology.residue(lid).name for lid in loop_residues])
        seq = "".join(seq)

        ll = len(seq)
        loop_G_percent = seq.count("G") / ll * 100
        loop_P_percent = seq.count("P") / ll * 100
        loop_S_percent = seq.count("S") / ll * 100
        loop_T_percent = seq.count("T") / ll * 100
        return dict(
            loop_seq=(seq, "Aminoacid sequence of the loop in one letter code."),
            loop_G_percent=(loop_G_percent, "Percent og Glycine residues in loop."),
            loop_P_percent=(loop_P_percent, "Percent of Proline residues in loop."),
            loop_S_percent=(loop_S_percent, "Percent of Serine residues in loop."),
            loop_T_percent=(loop_T_percent, "Percent of Threonine residues in loop."),
        )

    def get_loop_sasa(self, loop_index0):
        """Returns loop sasa , loop sasa in isolation and relative loop sasa"""
        loop_residues = self.loops0[loop_index0]

        loop_ids = self.topology.select(f"resid {loop_residues[0]} to {loop_residues[-1]}")

        loop_sasa_A = sum(self.sasa_atoms_A[loop_ids])  # get sasa just for loop
        loop_sasa_A_per_res = loop_sasa_A / len(loop_residues)

        # get SASA if the loop was on it's own, without the rest of the protein
        loop_isolation_traj = self.traj.atom_slice(loop_ids)
        loop_isolation_SASA_A = sum(md.shrake_rupley(loop_isolation_traj)[0] * 100)  # make in in angstrom
        loop_burial_percent = (1 - loop_sasa_A / loop_isolation_SASA_A) * 100
        loop_percent_of_total_surface = loop_sasa_A / self.total_sasa_A * 100

        return dict(
            loop_sasa_A=(loop_sasa_A, "Surface accessible area of loop in A**2"),
            loop_sasa_A_per_res=(
                loop_sasa_A_per_res,
                "Surface accessible area divided by number of residues of loop in A**2",
            ),
            loop_isolation_SASA_A=(
                loop_isolation_SASA_A,
                "Surface accessible area of loop without the rest of the protein in A**2",
            ),
            loop_burial_percent=(
                loop_burial_percent,
                "1-SASA/SASA_isolation, i.e. the percent of the loop surface covered by the rest of the protein",
            ),
            loop_percent_of_total_surface=(
                loop_percent_of_total_surface,
                "SASA/SASA_of_whole_protein, i.e. how big is the loop relative to the rest of the protein",
            ),
        )

    _loop_analyzers = [get_loop_geometry, get_loop_sasa, get_loop_sequence_features]

    def get_resi_features(self):
        self._resi_features = []
        for li, loop in enumerate(self.loops0):
            for ri, resi in enumerate(loop):
                self._resi_features.append(
                    dict(struct_name=self.struct_name, resi_loop_index0=ri, loop_index0=li, resi_index0=resi)
                )  # make new dict
                if li == 0 and ri == 0:  # description need to be added on first pass only
                    self.resi_feature_descriptions["struct_name"] = "Name of the structure"
                    self.resi_feature_descriptions["resi_index0"] = "The zero based index of the residue."
                    self.resi_feature_descriptions[
                        "resi_loop_index0"
                    ] = "The zero based index of the residue inside the loop."
                    self.resi_feature_descriptions["loop_index0"] = "The zero based index of the loop."

                for resi_analyzer in self._resi_analyzers:
                    res = resi_analyzer(self, li, ri, loop)
                    for f in res.keys():
                        # Add the values and descriptions to different lists
                        # Append to the last
                        self._resi_features[-1][f] = res[f][0]
                        if ri == 0:  # add the description if this is the first loop pass
                            self.resi_feature_descriptions[f] = res[f][1]

    def get_resi_geometry(self, loop_index0, resi_loop_index0, loop_residues):
        first_CA = self.topology.select(f"resid {loop_residues[0]} and name CA")[0]
        last_CA = self.topology.select(f"resid {loop_residues[-1]} and name CA")[0]
        residue_CA = self.topology.select(f"resid {loop_residues[resi_loop_index0]} and name CA")[0]
        # returns a set of frames, but we only have one frame, so [0] is needed
        resi_distance_to_N_term_A = md.compute_distances(self.traj, [[first_CA, residue_CA]])[0][0] * 10
        resi_distance_to_C_term_A = md.compute_distances(self.traj, [[last_CA, residue_CA]])[0][0] * 10

        return dict(
            resi_distance_to_N_term_A=(
                resi_distance_to_N_term_A,
                "distance of residue CA atom to start of loop (N_term first residue of loop) in Angstrom",
            ),
            resi_distance_to_C_term_A=(
                resi_distance_to_C_term_A,
                "distance of residue CA atom to end of loop (C_term last residue of loop) in Angstrom",
            ),
        )

    def get_resi_seq_features(self, loop_index0, resi_loop_index0, loop_residues):
        resi_index0 = loop_residues[resi_loop_index0]

        prev_resi_index0 = resi_index0 - 1
        next_resi_index0 = resi_index0 + 1

        resi_type = self.seq[resi_index0]
        resi_dssp = self.dssp[resi_index0]

        if prev_resi_index0 >= 0:
            prev_resi_type = self.seq[prev_resi_index0]
            prev_resi_dssp = self.dssp[prev_resi_index0]

        if next_resi_index0 < self.topology.n_residues:
            next_resi_type = self.seq[next_resi_index0]
            next_resi_dssp = self.dssp[next_resi_index0]

        return dict(
            resi_type=(resi_type, "Residue type"),
            resi_dssp=(resi_dssp, "Residue secondary structure"),
            prev_resi_type=(prev_resi_type, "Type of previous residue"),
            prev_resi_dssp=(prev_resi_dssp, "Secondary structure of previous residue"),
            next_resi_type=(next_resi_type, "Type of next residue"),
            next_resi_dssp=(next_resi_dssp, "Secondary structure of next residue"),
        )

    def get_resi_sasa(self, loop_index0, resi_loop_index0, loop_residues):
        resi_index0 = loop_residues[resi_loop_index0]

        resi_atoms = self.topology.select(f"resid {resi_index0}")

        resi_sasa_A = sum(self.sasa_atoms_A[resi_atoms])  # get sasa just for residue

        # get SASA if residue  was on it's own, without the rest of the protein
        resi_isolation_traj = self.traj.atom_slice(resi_atoms)
        resi_isolation_SASA_A = sum(md.shrake_rupley(resi_isolation_traj)[0] * 100)  # make in in angstrom
        resi_burial_percent = (1 - resi_sasa_A / resi_isolation_SASA_A) * 100
        resi_percent_of_total_surface = resi_sasa_A / self.total_sasa_A * 100

        return dict(
            resi_sasa_A=(
                resi_sasa_A,
                "Surface accessible area of residue in A**2 in the context of the protein",
            ),
            resi_isolation_SASA_A=(
                resi_isolation_SASA_A,
                "Surface accessible area of residue in A**2 in in isolation",
            ),
            resi_burial_percent=(
                resi_burial_percent,
                "1-SASA_residue/SASA_residue_isolation, i.e. the percent of the residue surface covered by the rest of the protein",
            ),
            resi_percent_of_total_surface=(
                resi_percent_of_total_surface,
                "SAS_resi/SASA_of_whole_protein, i.e. how big is the residue is relative to the rest of the protein",
            ),
        )

    def get_distances(self, resid0, target_resids0):
        """Returns an array of CA distances to resid0"""
        residue_CA = self.topology.select(f"resid {resid0} and name CA")[0]
        targets = [
            atom.index
            for atom in self.topology.atoms
            if ((atom.residue.index in target_resids0) and (atom.name == "CA"))
        ]

        pairs = [(residue_CA, target) for target in targets]
        dists_A = md.compute_distances(self.traj, pairs)[0] * 10
        return dists_A

    def get_active_geometry_res_info(self, loop_index0, resi_loop_index0, loop_residues):
        if not self.active_res_index0:
            return dict()
        resi_index0 = loop_residues[resi_loop_index0]
        dists = self.get_distances(resi_index0, self.active_res_index0)
        resi_active_site_dist_min_A = np.min(dists)
        resi_active_site_dist_avg_A = np.average(dists)
        resi_active_site_dist_max_A = np.min(dists)

        return dict(
            resi_active_site_dist_min_A=(
                resi_active_site_dist_min_A,
                "Minimum distance to one of the active site residues in A",
            ),
            resi_active_site_dist_avg_A=(
                resi_active_site_dist_avg_A,
                "Average distance to the active site residues in A",
            ),
            resi_active_site_dist_max_A=(
                resi_active_site_dist_max_A,
                "Maximum distance to one of the active site residues in A",
            ),
        )

    def get_active_seq_res_info(self, loop_index0, resi_loop_index0, loop_residues):
        if not self.active_res_index0:
            return dict()
        resi_index0 = loop_residues[resi_loop_index0]
        seq_dist = []
        for target in self.active_res_index0:
            seq_dist.append(abs(resi_index0 - target))

        resi_active_site_seq_dist_min = np.min(seq_dist)
        resi_active_site_seq_dist_avg = np.average(seq_dist)
        resi_active_site_seq_dist_max = np.min(seq_dist)

        return dict(
            resi_active_site_seq_dist_min=(
                resi_active_site_seq_dist_min,
                "Minimum distance (in number of residues) to one of the active site residues.",
            ),
            resi_active_site_seq_dist_avg=(
                resi_active_site_seq_dist_avg,
                "Average distance (in number of residues) to active site residues.",
            ),
            resi_active_site_seq_dist_max=(
                resi_active_site_seq_dist_max,
                "Maximum distance (in number of residues) to one of the active site residues.",
            ),
        )

    def get_resi_active_site_dssp_info(self, loop_index0, resi_loop_index0, loop_residues):
        if not self.active_res_index0:
            return dict()
        resi_index0 = loop_residues[resi_loop_index0]
        counts = []
        # dssps = []
        for target in self.active_res_index0:
            start, end = sorted((target, resi_index0))
            dssp_subset = self.dssp[start + 1 : end]
            # dssps.append()
            counter = Counter(dssp_subset, H=0, L=0, E=0)  # init counts to 0
            counts.append(counter)

        cdf = pd.DataFrame(counts)

        return dict(
            resi_active_site_num_H_min=(
                cdf.H.min(),
                "Minimum number of helical residues between this residue and one of the active site residues,",
            ),
            resi_active_site_num_H_avg=(
                cdf.H.mean(),
                "Average number of helical residues between this residue and one of the active site residues,",
            ),
            resi_active_site_num_H_max=(
                cdf.H.max(),
                "Maximum number of helical residues between this residue and one of the active site residues,",
            ),
            resi_active_site_num_E_min=(
                cdf.E.min(),
                "Minimum number of extended (beta) residues between this residue and one of the active site residues,",
            ),
            resi_active_site_num_E_avg=(
                cdf.E.mean(),
                "Average number of extended (beta) residues between this residue and one of the active site residues,",
            ),
            resi_active_site_num_E_max=(
                cdf.E.max(),
                "Maximum number of extended (beta) residues between this residue and one of the active site residues,",
            ),
            resi_active_site_num_L_min=(
                cdf.L.min(),
                "Minimum number of loop residues between this residue and one of the active site residues,",
            ),
            resi_active_site_num_L_avg=(
                cdf.L.mean(),
                "Average number of loop residues between this residue and one of the active site residues,",
            ),
            resi_active_site_num_L_max=(
                cdf.L.max(),
                "Maximum number of loop residues between this residue and one of the active site residues,",
            ),
        )

    _resi_analyzers = [
        get_resi_geometry,
        get_resi_seq_features,
        get_resi_sasa,
        get_active_geometry_res_info,
        get_active_seq_res_info,
        get_resi_active_site_dssp_info,
    ]
