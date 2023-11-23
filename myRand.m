%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Name of the script: myRand.m
% 
% Description: Function that produces random variable x based on probability
%  density function p(x)
%  Inputs: The function requires no inputs.
%  Outputs: The function outputs a random variable x based on the
%  probability distribution.
% 
% Name: Philip Ng
% UID: 305572506
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function x = myRand()
    y = rand();
    x = 5 - 5 * sqrt(1 - y);
end
