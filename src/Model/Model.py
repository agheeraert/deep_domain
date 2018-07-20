import torch
import torch.nn as nn
import numpy as np
from random import randint
from utils.coevolution import Coev

class DeepDomain(nn.Module):
    """
    Our deep learning model class:
        Batch_size = Batch size
        q = number of amino acids + 1 for gaps
        H = number of hidden states in the GRU that influence most layers as well
        p = dropout parameter for the coevolution input
    """
    def __init__(self, batch_size=1, q=21, H=512, p=0.5, gpu=False):
        """
        Initiate layers and hyperparameters
        The first layer is a bidirectionnal GRU that reads the along the dimension of size q and creates H*2 hidden states (*2 because it's bidirectionnal)
        Then a fully connected layer (which is equivalent to a convolutionnal layer with a kernel size of 1)
        The pairwise concatenation isn't a real DL layer so it's not defined here
        A dropout is applied to the coevolution input
        A second fully connected layer (equivalent to a convolutionnal layer with a kernel size of (1,1))
        Finally the softmax layer
        """
        super(DeepDomain, self).__init__()
        self.q = q
        self.H = H
        self.batch_size = batch_size
        
        #Gated Recurrent Unit
        self.gru = nn.GRU(q+1, H, batch_first=True, bidirectional=True)
                
        # Fully connected layer
        self.fully_connected1 = nn.Conv1d(2*H, H, kernel_size=1)
        
        #Dropout layer
        self.dropout = torch.nn.Dropout2d(p=0.5)
        
        #Second fully connected layer
        self.fully_connected2 = nn.Conv2d(2*H+1, 2, kernel_size=(1,1))
        
        #Final softmax
        self.softmax = nn.Softmax(dim=1)
        
    def forward(self, fas, mat, hidden=None):
        """
        the forward motion is performed thanks to four different kinds of input:
            fas: the pytorch sequence representation 
                 Size: (batch size, protein length, 20)
                 Size: (batch size, protein length, protein length)
            hidden: the hidden matrice used in the gru that is optimized at each step 
                 Size: (2, batch size, number of hidden states)
                 (there is a 2 because it is bidirectionnal)
        the final output is:
            The guess, a tensor of size (batch_size, protein length, protein length)
            The new hidden matrix of size (2, batch size, number of hidden states)
        """
        #Generating a random hidden matrix if no one is provided
        if not torch.is_tensor(hidden):
            hidden = torch.randn(2, 1, self.H)
        guesses = []
        for i in range(self.batch_size):
            x = fas[i].unsqueeze(dim=0).transpose(1,2)
            x, hidden = self.gru(x, hidden)
            x.transpose_(1,2)
            x = self.fully_connected1(x)        
            x = torch.cat([torch.einsum('hij,klm->hjmi', (x, x)), torch.einsum('hij,klm->kmjl', (x,x)), self.dropout(mat[i]).unsqueeze(dim=0).unsqueeze(dim=3)], dim=3)
            x = x.transpose(1, 3)
            x = self.fully_connected2(x)
            x = self.softmax(x)
            guesses.append(x)
        return guesses, hidden


