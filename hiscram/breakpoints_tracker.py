import numpy as np

class BPTracker(object):
    """
    Helps to track the breakpoints during the scrambling of a genome.
    Is suppose to yield a correct list of the breakpoints at base pair in the SCRAMBLED genome,
    not in the reference genome.
    """
    def __init__(self):
        self.bps = []
        self.pairs = [] ## To be implemented, maybe
    
    def INV(self, start:int, end:int):
        self.bps.append(start)
        self.bps.append(end)
        for i in range(len(self.bps)):
            if self.bps[i] > start and self.bps[i] < end:
                self.bps[i] = end - (self.bps[i] - start)

    def DEL(self, start:int, end:int):
        self.bps.append(start)
        for i in range(len(self.bps)):
            if self.bps[i] > start and self.bps[i] < end:
                self.bps[i] = start
            if self.bps[i] > end:
                self.bps[i] = self.bps[i] - (end - start)
    

    def TRA(self, start:int, end:int, target:int, invert:bool):
        # target is the position of the insertion IN THE DELETED GENOME
        
        # First, treat Deletion
        self.bps.append(start)
        bp_to_translocate = []
        for i in range(len(self.bps)):
            if self.bps[i] > start and self.bps[i] < end:
                bp_to_translocate.append(self.bps[i])
                self.bps[i] = start
            if self.bps[i] > end:
                    self.bps[i] = self.bps[i] - (end - start)
        
        # Then treat insertion
        ## The bps after the insertion
        self.bps.append(target)
        for i in range(len(self.bps)):
            if self.bps[i] > target:
                self.bps[i] += end - start

        ## The bps in the translocation
        for bp in bp_to_translocate:
            if invert:
                self.bps.append(target + end - (bp - start))
            else:
                self.bps.append(target + (bp - start))
        self.bps.append(target + end - start)
    
    def get_bps(self):
        return np.unique(self.bps)
        

