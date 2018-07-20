from os.path import join
from os import listdir
import numpy as np
import torch
from glob import glob
from Bio.PDB.Polypeptide import aa1
from utils.coevolution import Coev
from setup import folder
from tqdm import tqdm



DATA_DIR = join(folder, 'dataset')
ALN_DIR = join(DATA_DIR, 'aln/')
FAS_DIR = join(DATA_DIR, 'fas/')
MAT_DIR = join(DATA_DIR, 'mat/')

aa = list(aa1)
aa.append('-')
q = len(aa)
aa_to_idx = {}
for i, a in enumerate(aa):
    aa_to_idx[a] = i

for file in tqdm(listdir(FAS_DIR)):
    with open(join(FAS_DIR, file), "r") as seq_file:
        first = True
        _seq = []
        for line in seq_file:
            if not first:
                _seq += list(line)[:-1]
            else:
                first = False
        seq = []
        for elt in _seq:
            seq.append(aa_to_idx[elt])
        seq = np.asarray(seq, dtype='int')
        _ = []
        for k in range(q):
            _.append(np.equal(seq, k*np.ones_like(seq))*1)
        seq = np.stack(_)
        seq = torch.from_numpy(seq).float()
        torch.save(seq, join(DATA_DIR, file))
        
for file in tqdm(listdir(ALN_DIR)):
    with open(join(ALN_DIR, file), "r") as aln_file:
        _msa = []
        for line in aln_file:
            _msa.append(list(line)[:-1])
        B = len(_msa)
        L = len(_msa[0])
        msa = np.zeros([B, L], dtype='int')
        for i in range(L):
            for b in range(B):
                if _msa[b][i] in aa_to_idx:
                    msa[b, i] = aa_to_idx[_msa[b][i]]
                else:
                    msa[b, i] = q
        msa = torch.from_numpy(msa).float()
        _seq = torch.load(join(DATA_DIR, file[:-4]+'.fas'))
        _column = Coev(msa).entropy_line().float()
        seq = torch.cat([_seq, _column.unsqueeze(dim=0)],dim=0).float()
        torch.save(msa, join(DATA_DIR, file))
        torch.save(seq, join(DATA_DIR, file[:-4]+'.fas'))

for file in tqdm(listdir(MAT_DIR)):
    with open(join(MAT_DIR, file), "r") as mat_file:
        _mat = []
        for line in mat_file:
            _mat.append(line.split())
        mat = torch.from_numpy(np.asarray(_mat, dtype=float)).float()
        torch.save(mat, join(DATA_DIR, file))
        
