import os
import torch
from SLmodel_EL import load_checkpoint
from train_common_EL import get_device, setup_problem, train_and_eval

# Adjust this as needed (env overrides let the pipeline tag per test/grid)
inputfile = os.environ.get('EL_INPUT_TENSOR',  'input_tensor.csv')
outputfile = os.environ.get('EL_OUTPUT_TENSOR', 'output_tensor.csv')
dataname = os.environ.get('EL_DATANAME',      'JX_NN')

b_size     = 100
num_epochs = 1000
ln_rate    = 1e-4

modelname      = dataname + '_model.pth'
training_error = dataname + '_error.mat'

device = get_device()

dataloader_train, dataloader_test, model, criterion = \
    setup_problem(inputfile, outputfile, b_size, device)

# Load from checkpoint
optimizer = torch.optim.Adam(model.parameters(), lr=ln_rate)
model, optimizer, start_epoch, train_loss_history, test_loss_history = \
    load_checkpoint(model, optimizer, modelname)
# Update learning rate or the optimizer
optimizer = torch.optim.Adam(model.parameters(), lr=ln_rate)

model, optimizer, last_epoch, train_loss_history, test_loss_history = train_and_eval(
    model=model,
    optimizer=optimizer,
    dataloader_train=dataloader_train,
    dataloader_test=dataloader_test,
    criterion=criterion,
    num_epochs=num_epochs,
    start_epoch=start_epoch,
    train_loss_history=train_loss_history,
    test_loss_history=test_loss_history,
    modelname=modelname,
    training_error=training_error,
)
