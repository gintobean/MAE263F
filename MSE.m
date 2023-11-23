%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Name of the script: MSE.m
% 
% Description: This function calculates the mean squared error given an
% input array, output array, and weight and bias arrays.
%  Inputs: The function requires input/output arrays and weights and biases.
%  Outputs: The script outputs the mean squared error.
% 
% Name: Philip Ng
% UID: 305572506
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function L = MSE(X,y,W,b)
    N = size(y,1);
    Lsum = zeros(N,1);
    for i = 1:N
        yhat_i = X(i,:) * W + b;
        Lsum(i) = (y(i) - yhat_i)^2;
    end
    L = 1/N * sum(Lsum);
end