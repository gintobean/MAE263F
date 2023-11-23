%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Name of the script: myInitialization.m
% 
% Description: Function that uses randomization to define weights and
% biases
%  Inputs: The function recieves no input parameters
%  Outputs: The function outputs weights and biases in a 9 X 1 vector.
% 
% Name: Philip Ng
% UID: 305572506
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function WeightsBiases = myInitialization()
    W111 = rand(1,1);
    W121 = rand(1,1);
    W112 = rand(1,1);
    W122 = rand(1,1);
    W211 = rand(1,1);
    W221 = rand(1,1);
    b11 = rand(1,1);
    b12 = rand(1,1);
    b21 = rand(1,1);
    WeightsBiases = [W111;W121;W112;W122;W211;W221;b11;b12;b21];
end
