import os
from os.path import join
from setup import folder
import torch
import lzma
import numpy as np
from tqdm import tqdm
from scipy.spatial.distance import cdist

DATA_DIR = join(folder, 'dataset')
PDB_XZ_DIR = join(DATA_DIR, 'pdbxz')
LOG_DIR = join(DATA_DIR, 'log')

if __name__ == '__main__':
    for filename in tqdm(os.listdir(PDB_XZ_DIR)):
        with lzma.open(os.path.join(PDB_XZ_DIR, filename), 'rt') as file:
            x, y, z = [], [], []
            for line in file:
                if line.split()[0] == 'ATOM' and line.split()[2]=='CA':
                    _x, _y, _z = line.split()[5:8]
                    x.append(_x)
                    y.append(_y)
                    z.append(_z)
            for elt in [x, y, z]:
                elt = np.asarray(elt)
            coords = np.stack([x,y,z])
            try:
                dist_map = cdist(coords.transpose(), coords.transpose())
                cmap = dist_map < 12.0
                cmap = torch.from_numpy(1*cmap)
                torch.save(cmap, os.path.join(DATA_DIR, filename[:-7]+'.cmap'))
            except ValueError:
                with open(join(LOG_DIR, 'badfiles.log'), 'w') as badfiles:
                    badfiles.write(filename+'\n')
