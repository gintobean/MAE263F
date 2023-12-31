%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Title: project_305572506.m
% Description: Main code for the feedforward neural network.
% Inputs: The functions and all of their parameters, including the hyperparameters
% Outputs: Figures visualizing loss vs. epochs.
%
% Name: Philip Ng
% UID: 305572506
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all; clc;

%Load Training and Testing Data
[X_train, Y_train, X_test, Y_test] = load_train_and_test_data();

%Define network architecture and hyperparameters

input_size = size(X_train, 1); % input size of the FNN
output_size = size(Y_train, 1);% output size (number of classes) of the FNN
neurons = 64; % neurons each layer
numLayer = 2;

lr = 0.01; % learning rate
epochs = 150; % epochs

layer_dims = zeros(1, numLayer + 2);
layer_dims(1) = input_size;
layer_dims(end) = output_size;

for i = 1:numLayer
    layer_dims(i+1) = neurons;
end

%Train the model
parameters = initialize_parameters(layer_dims);
fprintf('Initialize cost somewhere so that we do not have warnings.\n');

% Train the model using mini-batch gradient descent
m = size(X_train, 2);
batch_size = 64;
num_batches = floor(m / batch_size);
trainLoss = zeros(epochs, 1);
testAccuracy = zeros(epochs, 1);
for i = 1:epochs

    % Shuffle the training data
    indices = randperm(m);
    X_train = X_train(:, indices);
    Y_train = Y_train(:, indices);

    % Initialize cost matrix
    cost = zeros(num_batches, batch_size);

    % Train on mini-batches
    for j = 1:num_batches
        X_batch = X_train(:, (j-1)*batch_size+1:j*batch_size);
        Y_batch = Y_train(:, (j-1)*batch_size+1:j*batch_size);
        forward_pass = forward_propagation(X_batch, parameters);
        cost(j,:) = compute_cost(forward_pass{end}, Y_batch);
        gradients = backward_propagation(X_batch, Y_batch, parameters,forward_pass);
        parameters = update_parameters(parameters, gradients, lr);
    end

    Y_pred = predict(X_test, parameters);
    acc = accuracy(Y_pred, Y_test);
    
    % Print the cost each epoch
    fprintf('Loss after epoch %d: Training: %f\n', i, norm(cost));
    trainLoss(i) = norm(cost);
    testAccuracy(i) = acc;
end

%Evaluate the model on the test set
fprintf('Test accuracy: %f\n', testAccuracy(end));

%Visualize the training progress
visualize_history(epochs, trainLoss, testAccuracy, lr, numLayer);

function predict_single_image(parameters)
% Visualizes the input image and plots a bar graph of probabilities for
% each class
% Inputs:
%   parameters: a struct containing the weights and biases (W1, b1,
%   W2, b2, etc.)
% Output:
%   (None)

    % Load data
    load('test_images.mat');
    index = randi(length(pixel));
    X = reshape(pixel(:,:,index), [size(pixel, 1)*size(pixel, 2), 1]) / 255;
    
    forward = forward_propagation(X, parameters);
    probabilities = forward{end};
    
    figure();
    subplot(1, 2, 1)
    imshow(pixel(:,:,index))
    title('Input Image')
    subplot(1, 2, 2)
    bar(0:9, probabilities)
    xlabel('Classes');
    ylabel('Probability');
    title('Probability Distribution')
end

