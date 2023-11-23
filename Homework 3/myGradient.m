%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Name of the script: myGradient.m
% 
% Description: This function carries out gradient descent on a set of
% weights and biases for one iteration.
%  Inputs: The script requires a previous weight and bias, as well as a
%  learning rate, and input and output values.
%  Outputs: The script outputs arrays for weights and biases after one
%  iteration of gradient descent.
% 
% Name: Philip Ng
% UID: 305572506
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function gradient = myGradient(h11,h12,x,y_pred,y_true,WB, lr)

    %Find partial derivatives of error
    dEdy = -2*y_true + 2*y_pred;
    dydh11 = WB(5);
    dydh12 = WB(6);
    dh11dW111 = x(1)*exp(-(WB(1)*x(1) + WB(2)*x(2) + WB(7)))/(1 + exp(-(WB(1)*x(1) + WB(2)*x(2) + WB(7))))^2;
    dh11dW121 = x(2)*exp(-(WB(1)*x(1) + WB(2)*x(2) + WB(7)))/(1 + exp(-(WB(1)*x(1) + WB(2)*x(2) + WB(7))))^2;
    dh12dW112 = x(1)*exp(-(WB(3)*x(1) + WB(4)*x(2) + WB(8)))/(1 + exp(-(WB(3)*x(1) + WB(4)*x(2) + WB(8))))^2;
    dh12dW122 = x(2)*exp(-(WB(3)*x(1) + WB(4)*x(2) + WB(8)))/(1 + exp(-(WB(3)*x(1) + WB(4)*x(2) + WB(8))))^2;
    dydW211 = h11;
    dydW221 = h12;
    dh11b11 = 1;
    dh12b12 = 1;
    dydb21 = 1;
    
    dEdW111 = dEdy * dydh11 * dh11dW111;
    dEdW121 = dEdy * dydh11 * dh11dW121;
    dEdW112 = dEdy * dydh12 * dh12dW112;
    dEdW122 = dEdy * dydh12 * dh12dW122;
    dEdW211 = dEdy * dydW211;
    dEdW221 = dEdy * dydW221;
    dEdb11 = dEdy * dydh11 * dh11b11;
    dEdb12 = dEdy * dydh12 * dh12b12;
    dEdb21 = dEdy * dydb21;
    
    %Update weights and biases
    W111new = WB(1) - lr*dEdW111;
    W121new = WB(2) - lr*dEdW121;
    W112new = WB(3) - lr*dEdW112;
    W122new = WB(4) - lr*dEdW122;
    W211new = WB(5) - lr*dEdW211;
    W221new = WB(6) - lr*dEdW221;
    b11new = WB(7) - lr*dEdb11;
    b12new = WB(8) - lr*dEdb12;
    b21new = WB(9) - lr*dEdb21;

    %Output new weights and biases 
    gradient = [W111new; W121new; W112new; W122new; W211new; W221new; b11new; b12new; b21new];
end
