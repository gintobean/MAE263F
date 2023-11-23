%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Title: forward_propogation.m
% Description:
%   This function performs forward propogation on each layer 
% Inputs: matrix X and the weights and biases from parameters
% Outputs: matrix with activation values of each layer index
%
% Name: Philip Ng
% UID: 305572506
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function activations = forward_propagation(X,parameters)
    num_layers = length(parameters);  % Number of layers in the network
    A = X;
    activations = cell(1, num_layers + 1);  % Cell array to store activations
    activations{1} = X;  % Set input data as the first activation
    
    for i = 1:num_layers  % Iterate through each layer
        Z = parameters{i}.W*A + parameters{i}.b;  % Perform linear transformation
        
        if i == num_layers
            A = softmax(Z);  % Apply tanh2 activation for hidden layers
        else
            A = tanh2(Z);  % Apply softmax activation for the output layer
        end
        
        activations{i+1} = A;  % Store the activation for the current layer
    end
end
