import torch
import torch.optim as optim
from torch.autograd import Variable
import torch.nn as nn   
import atexit
import os
import sys

sys.path.append(os.path.join(os.path.dirname(__file__), ".."))
from Model import DeepDomain


class DeepDomainTrainer():
    """
    The class that perform each training loop
    """
    def __init__(self, batch_size=1, lr=0.001, weight_decay=0.0, lr_decay=0.0001, gpu = False):
        """
        Initiliaze all the optimizer hyperparameters mainly:
            wd: the L2 penalty
            lr: the learning rate
            lrd: the learning rate decay
        Initiliaze the model used
        Initiliaze the optimizer itself (Adam for now)
        Initiliaze the loss function (CrossEntropy)
        """
        self.batch_size = batch_size
        self.wd = weight_decay
        self.lr = lr
        self.lr_decay = lr_decay
        self.gpu = gpu
        self.model = DeepDomain(batch_size=batch_size, gpu=gpu)
        self.optimizer = optim.Adam(self.model.parameters(), lr=self.lr, weight_decay=self.wd)
        if gpu:
            self.model.cuda()
        self.log = None

        atexit.register(self.cleanup)
    
    def new_log(self, log_file_name):
        if not self.log is None:
            self.log.close()
        self.log = open(log_file_name, "w")

    def cleanup(self):
        if not self.log is None:
            self.log.close()

    def optimize(self, data, hidden):
        """
        Optimization step. 
        Input:
            data which contains
                fas: the pytorch sequence representation 
                     Size: (batch size, protein length, 21)
                mat: the pytorch ccmpred output representation
                     Size: (batch size, protein length, protein length)
            hidden: the hidden matrice used in the GRU that is optimized at each step 
                 Size: (2, batch size, number of hidden states)
                 (there is a 2 because the GRU is bidirectionnal)
        Output: loss, hidden data
        """
        #Resetting the gradient and the loss
        self.optimizer.zero_grad()
        loss = torch.zeros(1)
        
        #Retrieving the data from the dict
        fas, mat, cmap = data
        
        #Putting everything on the GPu and wrapping it in Variable to perfom operations on them
        if self.gpu:
            loss = loss.cuda()
            for group in [fas, mat, cmap]:
                for tens in group:
                    tens = tens.cuda()
        loss = Variable(loss)            
        for group in [fas, mat, cmap]:
            for tens in group:
                tens = Variable(tens)
#        fas, mat, cmap = Variable(fas), Variable(mat), Variable(cmap)
        
        #Computing the guess and the new hidden state
        guess, hidden = self.model(fas, mat, hidden)
        
        #Computing the weights for the loss function for each elt of the batch
        for i in range(self.batch_size):
            N_c = torch.sum(cmap[i]).item()
            N_nc = cmap[i].shape[0]**2-N_c
            weight = torch.FloatTensor([N_c, N_nc])
            criterion = nn.CrossEntropyLoss(weight=weight)
            loss = loss + criterion(guess[i], cmap[i].unsqueeze(dim=0))
        #Detaches the hidden data so that the graph can be computed again
        hidden = hidden.detach()
        
        #End of the optimization process
        loss.backward()       
        self.optimizer.step()
        return loss.data, hidden, guess

