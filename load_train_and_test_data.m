%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Title: load_train_and_test_data.m
% Description: Loads the train and test data in matrices with the right dimensions.
% Inputs: None
% Outputs: train and test, one hot encoded, normalized, and flattened.
%
% Name: Philip Ng
% UID: 305572506
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to load and preprocess training and testing data
function [X_train, Y_train, X_test, Y_test] = load_train_and_test_data()
    % Load the testing and training image and label data
    test_images = load('test_images.mat');
    test_labels = load('test_labels.mat');
    train_images = load('train_images.mat');
    train_labels = load('train_labels.mat');

    % Get the dimensions of the training images
    s1 = size(train_images.pixel);
    % Get the dimensions of the testing images
    s2 = size(test_images.pixel);
    
    % Reshape and normalize the training images
    X_train = reshape(train_images.pixel, [s1(1)*s1(2) s1(3)]) / 255;
    % Reshape and normalize the testing images
    X_test = reshape(test_images.pixel, [s2(1)*s2(2) s2(3)]) / 255;
    
    % Convert training labels to categorical format
    train_labels.label = categorical(train_labels.label);
    % One-hot encode the training labels
    Y_train = onehotencode(train_labels.label,1);
    
    % Convert testing labels to categorical format
    test_labels.label = categorical(test_labels.label);
    % One-hot encode the testing labels
    Y_test = onehotencode(test_labels.label,1);
end





