import torch
import numpy as np
import os
import pandas as pd
import glob
import sys
from torch.utils.data import Dataset
from tqdm import tqdm
sys.path.append(os.path.join(os.path.dirname(__file__), ".."))
from utils.coevolution import Coevolution
import random
random.seed(42)


from Bio.PDB.Polypeptide import aa1

class Loader(Dataset):
    """
    The dataset that loads everything
    """
    def __init__(self, directory, csv_file=None):
        """
        A csv file that describes the database is used to load the database in the directory
        The algorithm is more convenient if the csv file is in the same folder as the database
        and that its name is index.csv
        """
        self.directory = directory
        if csv_file == None:
            if os.path.isfile(os.path.join(directory, 'index.csv')):
                csv_file = os.path.join(directory, 'index.csv')
        self.fac = pd.read_csv(csv_file)        
        os.chdir(self.directory)            
        self.length = len(glob.glob('*.fas'))

            
    def __getitem__(self, idx):
        """ Each kind of data needed for the training is loaded in a dictionnary thanks to the 
        csv_file that describes our database.
        The data needed are:
            fas: the binary torch tensor that represents the sequence
            mat: the binary torch tensor that represents CCMpred input
            cmap: the binary torch tensor that represents the experimental contact map
        """
        os.chdir(self.directory)
        fas, mat, cmap = self.fac.iloc[idx, :]
        fas, mat, cmap = torch.load(fas), torch.load(mat), torch.load(cmap)
        sample = {'fas': fas, 'mat': mat, 'cmap': cmap}
        return sample
        
        
    def __len__(self):
        """
        Returns length of the dataset
        """
        return self.length
    
 
    
#Former useless bits of code that computed everything on the fly
        
#        self.data = {'input':[], 'label':[]}
#        os.chdir(self.directory)            
#        filelist = glob.glob('*.fas')
#        for fasfile in filelist:
#            msafile = fasfile[:-4]+'.aln'
#            with open(fasfile, "r") as seq_file:
#                first = True
#                seq = []
#                for line in seq_file:
#                    if not first:
#                        seq += list(line)[:-1]
#                    else:
#                        first = False
#                L = len(seq)
#            rep_seq = np.zeros([L, self.q])
#            for i, aa in enumerate(seq):
#                aa = self.aa_to_idx[aa]
#                rep_seq[i, aa] = 1
#            _seq = torch.from_numpy(rep_seq)
#            _column = Coevolution(msafile).entropy_line().unsqueeze(dim=1)
#            self.data['input'].append(torch.cat([_seq, _column], dim=1))
#            self.data['label'].append(torch.load(fasfile[:-4]+'.cmap'))
#        return self.data
    