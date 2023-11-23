%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Title: compute_cost.m
% Description: Computes the cross-entropy-loss of each sample.
% Inputs: Actual output data and predicted output layer data
% Outputs: Cross-entropy-loss
%
% Name: Philip Ng
% UID: 305572506
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function cost = compute_cost(AL, Y)
     cost = -sum(Y .* log(AL));  % Compute the cross-entropy loss
end
