%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Title: update_parameters.m
% Description:
%   This function performs gradient descent and updates the weights and
%   biases of each layer.
% Inputs: batch inputs X and Y, current parameters, and values of each layer
% Outputs: New parameters, gradient-descent optimized
%
% Name: Philip Ng
% UID: 305572506
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parameters = update_parameters(parameters, gradients, learning_rate)
    L = length(parameters);
    for i = 1:L-1
        parameters{i}.W = parameters{i}.W - learning_rate*gradients{i}.dW;
        parameters{i}.b = parameters{i}.b - learning_rate*gradients{i}.db;
    end
end