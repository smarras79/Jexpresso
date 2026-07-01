import torch
import torch.nn as nn

class FCNN(nn.Module):
    def __init__(self, input_size, output_size):
        super(FCNN, self).__init__()
        self.hidden_size = int(input_size)
        
        self.activation = nn.ReLU()
        #self.activation = nn.Tanh()
        #self.activation = nn.Sigmoid()
        #self.activation = nn.ELU()

        self.fully_connected_nets = nn.Sequential(
            nn.Linear(input_size, self.hidden_size),
            self.activation,
            nn.Linear(self.hidden_size, self.hidden_size),
            self.activation,
            nn.Linear(self.hidden_size, self.hidden_size),
            self.activation,
            nn.Linear(self.hidden_size, self.hidden_size),
            self.activation,
            nn.Linear(self.hidden_size, self.hidden_size),
            self.activation,
            nn.Linear(self.hidden_size, output_size),
        )

    def forward(self, x):
        x = self.fully_connected_nets(x)
        return x
    