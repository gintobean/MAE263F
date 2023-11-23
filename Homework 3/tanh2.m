%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Title: tanh2.m
% Description:
%   This function calculates the tanh activation function with the previous
%   layer values as input
% Inputs: matrix X, previous layer values
% Outputs: tanh2 activated values
%
% Name: Philip Ng
% UID: 305572506
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Function to apply the tanh activation function to the input matrix
function Z = tanh2(X)
    % Compute the tanh activation function element-wise for each element in X
    Z = 2 ./ (1 + exp(-2*X)) - 1;
end
