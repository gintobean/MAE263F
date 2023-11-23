%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Name of the script: myWeightsBiases.m
% 
% Description: Function that uses randomization to define weights and
% biases
%  Inputs: The function recieves no inputn parameters
%  Outputs: The function outputs weights and biases.
% 
% Name: Philip Ng
% UID: 305572506
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Function to generate random weights and biases based on table specs
function [W11,W12,W13,W21,W22,W23,W31,W32,b11,b12,b13,b21,b22,b23,b31,b32] = myWeightsBiases()
    W11 = 8 * rand(4,1) - 4;
    W12 = 9 * rand(4,1) - 4;
    W13 = 11 * rand(4,1) - 4;
    W21 = 16 * rand(3,1) - 6;
    W22 = 12 * rand(3,1) - 4;
    W23 = 10 * rand(3,1) - 6;
    W31 = 12 * rand(3,1) - 9;
    W32 = 5 * rand(3,1) + 1;
    b11 = rand(1,1);
    b12 = rand(1,1) - 1;
    b13 = rand(1,1) - 1;
    b21 = rand(1,1) - 1;
    b22 = rand(1,1);
    b23 = rand(1,1);
    b31 = rand(1,1) - 1;
    b32 = rand(1,1);
end

