%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Name of the script: myForwardPass.m
% 
% Description: Function that predicts y based on weights and biases previously defined
%  Inputs: The function takes the ground input, as well the weights
%  and biases vector, as inputs.
%  Outputs: The function outputs the values of the hidden and output layers of
%  the neural network.
% 
% Name: Philip Ng
% UID: 305572506
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [h11,h12,y_pred] = myForwardPass(WB, x)
    
    %Define sigmoid activation function
    sigmoid = @(x1,x2,W1,W2,b) 1/ (1 + exp(-(x1*W1 + x2*W2 + b)));

    h11 = sigmoid(x(1),x(2),WB(1),WB(2),WB(7)); %Activate first hidden layer
    h12 = sigmoid(x(1),x(2),WB(3),WB(4),WB(8));

    y_pred = h11*WB(5) + h12*WB(6) + WB(9); %Do not activate output layer
end
