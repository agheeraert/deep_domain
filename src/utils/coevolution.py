import numpy as np
from Bio.PDB.Polypeptide import aa1
from scipy.stats import entropy
import torch

class Coevolution():
    def __init__(self, filename):
        aa = list(aa1)
        aa.append('-')
        aa_to_idx = {}
        for i, a in enumerate(aa):
            aa_to_idx[a] = i
        self.q = len(aa)
        _msa = []
        with open(filename) as msa_file:
            for line in msa_file:
                _L = np.asarray(list(line)[:-1])
                _ = []
                for i, elt in enumerate(_L):
                    if elt in aa_to_idx:
                        _.append(aa_to_idx[elt])
                    else:
                        _.append(aa_to_idx['-'])
                _msa.append(_)    
        self.B = len(_msa)
        self.L = len(_msa[0])
        msa = np.stack(_msa)
        msa = msa.astype(np.int)
        self.msa = msa
                
    def frequency(self, i, k):
        line = self.msa[i,:]
        compare = k*np.ones(self.B)
        return 1./self.B*np.sum(np.eq(line, compare))
    
    def cofrequency(self, i, j, k, l):
        line1, line2 = self.msa[i:,], self.msa[j:,]
        compare1, compare2 = k*np.ones(self.B), l*np.ones(self.B)
        return 1./self.B*np.sum(np.eq(line1, compare1)*line2(compare2))
    
    def matrix_frequency(self):
        freq_mat = np.zeros([self.L, self.q])
        for k in range(self.q):
            compare = k*np.ones([self.B, self.L])
            freq_mat[:,k] = 1./self.B*np.sum(np.equal(self.msa, compare), axis=0)
        return freq_mat
    
    def create_comatrix(self):
        _ = []
        always_add = self.q*self.msa
        for i in range(1, self.L):
            tom = np.roll(self.msa, -i, axis=1)
            _.append(tom+always_add)
        return np.concatenate(_,axis=1)
    
    def matrix_cofrequency(self):
        """
        cofreq(x,y) = cofreq(i(L-1)+j, kq+l) = fijkl
        works but sadly useless
        """
        comatrix = self.create_comatrix()
        cofreq_mat = np.zeros([self.L*(self.L-1), self.q**2])
        for k in range(self.q):
            for l in range(self.q):
                compare = (k*self.q+l)*np.ones([self.B, self.L*(self.L-1)])
                cofreq_mat[:,k*self.q+l] = 1/self.B*np.sum(np.equal(comatrix, compare), axis=0)
        return cofreq_mat
    
    def entropy_line(self):
        return torch.from_numpy(entropy(np.swapaxes(self.matrix_frequency(), 0, 1),base=2))



            
class Coev():
    def __init__(self, aln):
        self.msa = aln.numpy()
        self.q = len(aa1)+1
        self.B, self.L = self.msa.shape
                
    def frequency(self, i, k):
        line = self.msa[i,:]
        compare = k*np.ones(self.B)
        return 1./self.B*np.sum(np.eq(line, compare))
    
    def cofrequency(self, i, j, k, l):
        line1, line2 = self.msa[i:,], self.msa[j:,]
        compare1, compare2 = k*np.ones(self.B), l*np.ones(self.B)
        return 1./self.B*np.sum(np.eq(line1, compare1)*line2(compare2))
    
    def matrix_frequency(self):
        freq_mat = np.zeros([self.L, self.q])
        for k in range(self.q):
            compare = k*np.ones([self.B, self.L])
            freq_mat[:,k] = 1./self.B*np.sum(np.equal(self.msa, compare), axis=0)
        return freq_mat
    
    def create_comatrix(self):
        _ = []
        always_add = self.q*self.msa
        for i in range(1, self.L):
            tom = np.roll(self.msa, -i, axis=1)
            _.append(tom+always_add)
        return np.concatenate(_,axis=1)
    
    def matrix_cofrequency(self):
        """
        cofreq(x,y) = cofreq(i(L-1)+j, kq+l) = fijkl
        works but sadly useless
        """
        comatrix = self.create_comatrix()
        cofreq_mat = np.zeros([self.L*(self.L-1), self.q**2])
        for k in range(self.q):
            for l in range(self.q):
                compare = (k*self.q+l)*np.ones([self.B, self.L*(self.L-1)])
                cofreq_mat[:,k*self.q+l] = 1/self.B*np.sum(np.equal(comatrix, compare), axis=0)
        return cofreq_mat
    
    def entropy_line(self):
        return torch.from_numpy(entropy(np.swapaxes(self.matrix_frequency(), 0, 1),base=2))
        
        
        
#if __name__ == '__main__':             
#    filename = '/home/aria/Stage/deep_domain/dataset/1BDO_A.aln'
#    msa = Coevolution(filename).msa
#    H = Coevolution(filename).entropy_line()
