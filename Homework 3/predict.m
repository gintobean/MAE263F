%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Title: predict.m
% Description: Main code for the feedforward neural network.
% Inputs: Input layer data X
% Outputs: Predicted values of output layer
%
% Name: Philip Ng
% UID: 305572506
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function Y_pred = predict(X, parameters)
    activations = forward_propagation(X, parameters);  % Perform forward propagation to compute activations
    Y_pred = activations{end} == max(activations{end});  % Convert the final activations to a numeric array
end
