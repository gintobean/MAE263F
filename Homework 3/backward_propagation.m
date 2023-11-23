%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Title: backward_propagation.m
% Description:
%   This function performs backward propagation on each layer.
% Inputs: matrix X, next layer values
% Outputs: The gradients to optimize activation values.
%
% Name: Philip Ng
% UID: 305572506
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function gradients = backward_propagation(X,Y,parameters,activations)
    L = length(parameters);  % Number of layers
    m = size(X,2);  % Number of training examples
    gradients = cell(1,L);  % Cell array to store gradients of parameters

    dZ = activations{end} - Y;  % Compute gradient of last layer
    gradients{L}.dW = (dZ * activations{L}') / m;  % Compute gradient of weights for last layer
    gradients{L}.db = sum(dZ,2)/m;  % Compute gradient of biases for last layer

    for i = L-1:-1:1  % Iterate backwards through layers
        dA = parameters{i+1}.W' * dZ;  % Compute gradient of activations of previous layer
        dZ = dA .* (1 - tanh2(activations{i+1}).^2);  % Compute new gradient of current layer
        
        if i == 1
            Aprev = X;  % If first layer, set previous activations as input data
        else
            Aprev = activations{i};  % Set previous activations as activations of previous layer
        end

        gradients{i}.dW = dZ * Aprev' / m;  % Compute gradient of weights for current layer
        gradients{i}.db = sum(dZ,2) / m;  % Compute gradient of biases for current layer
    end
end
