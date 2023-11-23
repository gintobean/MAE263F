%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Title: initialize_parameters.m
% Description:
%   This function initializes the parameters (weights and biases) 
%   for a neural network with multiple layers. It takes the dimensions 
%   of each layer as input and returns the initialized parameters 
%   as output.
% Inputs: layer dimensions
% Outputs: Weights and biases for each layer
%
% Name: Philip Ng
% UID: 305572506
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Function to initialize parameters for neural network layers
function parameters = initialize_parameters(layer_dims)

    % Find the total number of layers
    num_layers = length(layer_dims);

    % Initialize the parameters as a cell array
    parameters = cell(1, num_layers-1);
    
    % For each layer
    for i = 1:num_layers - 1
        % Initialize weights randomly for the current layer using standard normal distribution
        parameters{i}.W = randn(layer_dims(i+1), layer_dims(i));
        
        % Initialize biases as zeros for the current layer
        parameters{i}.b = zeros(layer_dims(i+1),1);
    end
end