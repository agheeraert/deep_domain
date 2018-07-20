import torch
from os.path import join
import argparse
from setup import folder
from Training import DeepDomainTrainer
from Dataset.Loader import Loader
from tqdm import tqdm
import numpy as np
import matplotlib.pyplot as plt
from torch.utils.data import DataLoader

"""
Main training file. It's this file that has to be launch in order to train the model 
"""

DATA_DIR = join(folder, 'dataset')
 
def my_collate(batch):
    """ custome collate function to create the batch"""
    fas = [item['fas'] for item in batch]
    mat = [item['mat'] for item in batch]
    cmap = [item['cmap'] for item in batch]
    return [fas, mat, cmap]         
        
        

if __name__=='__main__':
    parser = argparse.ArgumentParser(description='Deep Domain training')    
    
    parser.add_argument('-lr', type=float,  default=0.00001 , help='Learning rate')
    parser.add_argument('-lrd', type=float, default=0.000075, help='Learning rate decay')
    parser.add_argument('-wd', type=float, default=0.0, help='Weight decay')
    parser.add_argument('-max_epoch', type=int, default=100, help='Max epoch')
    parser.add_argument('-gpu', default=None, help='Use gpu')
    parser.add_argument('-bs', default=2, help='Batch size')
    parser.add_argument('-gpu_num', type=int, default=1, help='GPU number')
    
    #Retrieving arguments
    args = parser.parse_args()
    
    #If the training is perfomed on the GPU we set the GPU number
    if args.gpu is None:
        gpu = False
    else:
        gpu = True
        torch.cuda.set_device(int(args.gpu_num))
    
    """
    Generation of everything that is needed to perform the training loops
    data: our data
    hidden: thanks to this initialization the hidden matrix from the GRU will be initialized further
    trainer: our trainer (see Training/trainer.py)
    loss_list: to draw the losses plot
    """                          
    training_set = Loader('/home/aria/Stage/deep_domain/dataset')
    hidden = None
    training_loader = DataLoader(training_set, batch_size = args.bs, collate_fn=my_collate)
    trainer = DeepDomainTrainer(batch_size=args.bs, lr=args.lr, weight_decay=args.wd, lr_decay=args.lrd, gpu=gpu)
    loss_list = []
    
    #Training loops

    for epoch in tqdm(range(args.max_epoch)):
        av_loss = 0.0
        for elt in training_loader:
            loss, hidden, guesses = trainer.optimize(elt, hidden)
            av_loss += loss
        loss_list.append(av_loss)
        
    losses = np.asarray(loss_list)
    plt.plot(np.arange(args.max_epoch), losses)
    

# padding
#def my_collate(batch):
#    _fas = [item['fas'] for item in batch]
#    _mat = [item['mat'] for item in batch]
#    _cmap = [item['cmap'] for item in batch]
#    L_max = 0
#    for seq in _fas:
#        L = seq.shape[1]
#        if L > L_max:
#            L_max = L
#    fas, mat, cmap = [], [], []
#    for i, _seq in enumerate(_fas):
#        L = _seq.shape[1]
#        if L != L_max:
#            fill1 = torch.zeros([L, L_max - L])
#            fill2 = torch.zeros([L_max-L, L_max])
#            fas.append(torch.cat([torch.cat([_seq, fill1]), fill2]))
#            mat.append(torch.cat([torch.cat([_mat[i], fill1])], fill2))
#            cmap.append(torch.cat([torch.cat([_cmap[i], fill1])], fill2))
#            print(torch.cat([torch.cat([_seq, fill1]), fill2]).shape)
#        else:
#            fas.append(_seq)
#            mat.append(_mat[i])
#            cmap.append(_cmap[i])
#    del _fas, _mat, _cmap
#    fas = torch.stack(fas, dim=0)
#    mat = torch.stack(mat, dim=0)
#    cmap = torch.stack(cmap, dim=0)         


