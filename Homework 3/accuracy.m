%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Title: accuracy.m
% Description: Accuracy calculation for predictions of each image.
% Inputs: Predictions for each image and actual classification for each
% image
% Outputs: Accuracy calculation.
%
% Name: Philip Ng
% UID: 305572506
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function acc = accuracy(Y_pred, Y)
      acc = mean(all(Y_pred == Y, 1));
end
