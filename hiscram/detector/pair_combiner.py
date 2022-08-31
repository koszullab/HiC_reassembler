from math import exp
from operator import index
import os
import pypairix
import numpy as np
import sys
from hiscram.reassemble.svs import SVs
from tqdm import tqdm
from typing import Dict, List

class PairCombiner(object):
    """
    Provides a class to combine the detected breakpoints in order to form structural variations.
    The combining is based on a score matrix. Each side of each breakpoint is associated with
    another side of a breakpoint.
    Following this association, the SV are determined.
    For example if 3 and 4 are 2 different breakpoints.
    3+/4+ and 3-/4- -> 3-4 is INVersion.
    """
    def __init__(self, chrom, detected_bps, pairfile, binsize, readlength, tmpdir:str="./tmpdir"):
        self.bp_to_associate = detected_bps
        self.score_matrix = np.zeros((2*len(detected_bps),2*len(detected_bps)))
        self.associated_sv = []
        
        self.sv_name = []
        self.sv_type = []
        self.coordsBP1 = []
        self.coordsBP2 = []
        self.coordsBP3 = []
        self.sv_size = []
        self.sv_index = 0

        self.chrom = chrom
        self.pairfile = pairfile
        self.tmpdir = tmpdir
        self.binsize = binsize
        self.readlength = readlength
        self.load_pairs()
        
    def load_pairs(self):
        # pairix only works with compressed pair files. So we use bgzip to compress them
        try:
            os.mkdir(self.tmpdir)
        except:
            pass
        os.system(f"pairtools sort {self.pairfile} > {self.tmpdir}/{self.pairfile.parts[-1]}.sorted")
        os.system(f"bgzip {self.tmpdir}/{self.pairfile.parts[-1]}.sorted -c > {self.tmpdir}/{self.pairfile.parts[-1]}.sorted.gz")
        os.system(f"pairix {self.tmpdir}/{self.pairfile.parts[-1]}.sorted.gz")
        self.pairs_data = pypairix.open(f"{self.tmpdir}/{self.pairfile.parts[-1]}.sorted.gz")
    
    def create_coverage(self, it_for, it_rev) -> Dict : 
        """
        Computes the positions of the associated pair reads as a coverage.
        Returns a Dictionary. Each position is a key and the associated word is the number
        of pair reads mapping at this position.
        """
        coverage = {}
        for read in it_for:
            pair_loc = int(read[4])
            sign = read[6]
            for i in range(self.readlength):
                if sign == '+':
                    index = pair_loc + i
                else:
                    index = pair_loc - i
                if coverage.get(index) is None:
                    coverage[index] = 0
                coverage[index] += 1

        for read in it_rev:
            pair_loc = int(read[2])
            sign = read[5]
            for i in range(self.readlength):
                if sign == '+':
                    index = pair_loc + i
                else:
                    index = pair_loc - i
                if coverage.get(index) is None:
                    coverage[index] = 0
                coverage[index] += 1
        return coverage

    def find_closest_bp(self, position):
        distances = [(position - bp) for bp in self.bp_to_associate]
        min_index = np.argmin(np.absolute(distances))
        return min_index, distances[min_index]

    # Gives more weight to pair reads close to a breakpoint
    def score(self,n,d):
        # Dis is based on exp (seems to work pretty well)
        return exp(-abs(d)/n**2)
        

        """
        # Dis is the basic invert function
        if d == 0 : return 2*n
        else: return abs(n/d)
        """

    def associate_one_sv(self,bp_index):
        """
        Compute the 2 lines of the score matrix corresponding to the 2 sides of the given breakpoint.
        """
        shift_window = 0
        window_size = 1000

        bp = self.bp_to_associate[bp_index]

        start_query_left = bp - shift_window - window_size
        end_query_left = bp - shift_window
        query_left = f"{self.chrom}:{start_query_left}-{end_query_left}"
        start_query_right = bp + shift_window
        end_query_right = bp + shift_window + window_size
        query_right = f"{self.chrom}:{start_query_right}-{end_query_right}"

        it_for_left = self.pairs_data.querys2D(f"{query_left}|*")
        it_rev_left = self.pairs_data.querys2D(f"*|{query_left}")

        it_for_right = self.pairs_data.querys2D(f"{query_right}|*")
        it_rev_right = self.pairs_data.querys2D(f"*|{query_right}")

        # for each side + - of the bp, create the pair coverage
        coverage_left = self.create_coverage(it_for_left, it_rev_left)
        coverage_right = self.create_coverage(it_for_right, it_rev_right)

        # Score each bp as a function of their proximity to the pairs
        for position in coverage_left:
            closest_index, distance_to_closest = self.find_closest_bp(position)
            if(distance_to_closest < 0):
                self.score_matrix[2*bp_index][2*closest_index] += self.score(coverage_left[position],distance_to_closest)
            if(distance_to_closest > 0):
                self.score_matrix[2*bp_index][2*closest_index + 1] += self.score(coverage_left[position],distance_to_closest)
                    
        for position in coverage_right:
            closest_index, distance_to_closest = self.find_closest_bp(position)
            if(distance_to_closest < 0):
                self.score_matrix[2*bp_index + 1][2*closest_index] += self.score(coverage_right[position],distance_to_closest)
                # self.score_matrix[2*closest_index][2*bp_index + 1] += round(abs(coverage_right[position]/distance_to_closest),2)
            if(distance_to_closest > 0):
                self.score_matrix[2*bp_index + 1][2*closest_index + 1] += self.score(coverage_right[position],distance_to_closest)
                # self.score_matrix[2*closest_index + 1][2*bp_index + 1] += round(abs(coverage_right[position]/distance_to_closest),2)

    def create_score_matrix(self):
        """
        Compute the whole score matrix, for all the detected breapoints.
        """
        for i in tqdm(range(len(self.bp_to_associate))):
            self.associate_one_sv(i)
        for i in range(2*len(self.bp_to_associate)):
            self.score_matrix[i][i] = 0
    
    def solve_matrix(self,debug=False):
        """
        Associate the breapoints following the score matrix.
        """
        threshold = 0.5
        self.mat = np.array(self.score_matrix)
        print("score matrix =\n", np.array(self.mat,dtype=np.int32))
        while True:
            if(debug):
                print("score matrix :\n", np.array(self.mat, dtype=np.int32))
            # Find the highest association score
            maxindexes = np.argmax(self.mat, axis=1)
            maxes = np.max(self.mat, axis = 1)
            line = np.argmax(maxes)
            col = maxindexes[line]
            score = self.mat[line, col]

            if(debug):
                print("maxindexes = ", maxindexes)
                print("maxes = ", maxes)
                print("line = ", line)
                print("col = ", col)

            # Condition to stop the solving: the scores left are too low
            print("INFO : score = ", score)
            if(score < threshold):
                break

            # Faut regarder le score max de l'autre orientation de chaque bp
            if(line%2 == 0): oo_line = line + 1
            else: oo_line = line - 1
            if(col%2 == 0): oo_col = col + 1
            else: oo_col = col -1

            # Max of other orientation
            ass_line = maxindexes[oo_line]
            ass_line_score = maxes[oo_line]
            ass_col = maxindexes[oo_col]
            ass_col_score = maxes[oo_col]

            if(debug):
                print("ass_line = ", ass_line)
                print("ass_col = ", ass_col)
                print("ass_line_score = ", ass_line_score)
                print("ass_col_score = ", ass_col_score)
            
            # Interprete max score
            bp1_index = oo_line // 2 # (also = line // 2)
            bp2_index = oo_col // 2 # (also = col // 2)
            bp3_index = ass_line // 2 # (also = ass_col // 2)
            if(bp1_index == bp2_index):
                self.associate_DEL(bp1_index)
                if(debug): print(f"Added DEL at index {bp1_index} : {self.bp_to_associate[bp1_index]}")
            elif(ass_line == oo_col):
                self.associate_INV(bp1_index, bp2_index)
                if(debug): print(f"Added INV at indexes {bp1_index}-{bp2_index} : {self.bp_to_associate[bp1_index]}-{self.bp_to_associate[bp2_index]}")
            elif(abs(ass_line-ass_col) == 1):
                # In that case, we have a TRAnslocation
                self.associate_TRA(bp1_index,bp2_index, bp3_index)
                if(debug): print(f"Added TRA at indexes {bp1_index}-{bp2_index}-{bp3_index} : {self.bp_to_associate[bp1_index]}-{self.bp_to_associate[bp2_index]}-{self.bp_to_associate[bp3_index]}")
            else:
                self.associate_DEL(bp1_index)
                if(debug): print(f"Added DEL at index {bp1_index} : {self.bp_to_associate[bp1_index]}")


    def create_SVs_info(self) -> SVs:
        self.info_sv = SVs(
            np.array(self.sv_name),
            np.array(self.sv_type),
            np.array(self.coordsBP1),
            np.array(self.coordsBP2),
            np.array(self.coordsBP3),
            np.array(self.sv_size),
        ) # No sign because not implemented yet


    def print_sv(self):
        for sv in self.associated_sv:
            print(sv)
    
    
    def associate_DEL(self, index_bp1):
        # Size of a DEL is unkown yet. We set it a 0 for now.
        # There won't be any reassembly of this SV.
        self.associated_sv.append([self.bp_to_associate[index_bp1], "DEL"])
        # Ignore DEL for the moment
        """
        self.add_SV(
            self.sv_index,
            "INS",
            self.bp_to_associate[index_bp1],
            -1,
            -1,
            0,
        )
        """
        self.delete_bp(index_bp1)
    
    def associate_INV(self, index_bp1, index_bp2):
        self.associated_sv.append([self.bp_to_associate[index_bp1], self.bp_to_associate[index_bp2], "INV"])
        self.add_SV(
            self.sv_index,
            "INV",
            self.bp_to_associate[index_bp1],
            self.bp_to_associate[index_bp2],
            -1,
            abs(self.bp_to_associate[index_bp1]-self.bp_to_associate[index_bp2])
        )
        self.delete_bp(index_bp1)
        self.delete_bp(index_bp2)
    
    def associate_TRA(self, index_bp1, index_bp2, index_bp3):
        # Sort the breakpoints by position
        bp_sorted = np.sort([
            self.bp_to_associate[index_bp1], 
            self.bp_to_associate[index_bp2], 
            self.bp_to_associate[index_bp3]
        ])
        bp1 = bp_sorted[0]
        bp2 = bp_sorted[1]
        bp3 = bp_sorted[2]

        # Record the associated bp
        self.associated_sv.append(
            [
                bp1,
                bp2,
                bp3,
                "TRA"
            ]
        )
        self.add_SV(
            self.sv_index,
            "TRA_back",
            bp1,
            bp2,
            bp3,
            abs(bp1-bp2)
        )

        self.delete_bp(index_bp1)
        self.delete_bp(index_bp2)
        self.delete_bp(index_bp3)

    def delete_bp(self, i):
        N = len(self.mat)
        for k in range(N):
            self.mat[2*i][k] = 0
            self.mat[2*i+1][k] = 0
            self.mat[k][2*i] = 0
            self.mat[k][2*i+1] = 0

    def add_SV(self, name, type, coord1, coord2, coord3, size):
        self.sv_name.append(name)
        self.sv_type.append(type)
        self.coordsBP1.append(coord1)
        self.coordsBP2.append(coord2)
        self.coordsBP3.append(coord3)
        self.sv_size.append(size)
        self.sv_index += 1
        
    def combine(self):
        self.create_score_matrix()
        self.solve_matrix()
        self.create_SVs_info()
        self.print_sv()
        self.save_svs()
        return self.info_sv

    def save_svs(self):
        with open("data/output/assembled_sv.txt", 'w') as sv_file:
            for sv in self.associated_sv:
                for i in sv:
                    sv_file.write(f"{i} ")
                sv_file.write("\n")
