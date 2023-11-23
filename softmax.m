%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Title: softmax.m
% Description:
%   This function calculates the softmax of the output layer
% Inputs: matrix X, output layer values
% Outputs: softmax of output layer
%
% Name: Philip Ng
% UID: 305572506
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Z = softmax(X)
    %Calculate the sum of exponential values for each element in X
    sumExp = sum(exp(X));
    %Divide each element's exponential value by the sum to obtain softmax
    Z = exp(X) ./ sumExp;
end

