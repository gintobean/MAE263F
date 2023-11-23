function visualize_history(epochs, trainLoss, testAccuracy, learning_rate, numLayer)
    % Plotting training loss vs. epochs
    subplot(2, 1, 1);  % Create the first subplot
    plot(1:epochs, trainLoss, 'b-');
    xlabel('Epochs');
    ylabel('Training Loss');
    title(sprintf('Training Loss vs. Epochs\nLearning Rate: %g, Hidden Layers: %d', learning_rate, numLayer));

    % Plotting testing accuracy vs. epochs
    subplot(2, 1, 2);  % Create the second subplot
    plot(1:epochs, testAccuracy, 'r-');
    xlabel('Epochs');
    ylabel('Testing Accuracy');
    title(sprintf('Testing Accuracy vs. Epochs\nLearning Rate: %g, Hidden Layers: %d', learning_rate, numLayer));

    % Save the figure as a PNG file
    filename = sprintf('model %.3f %d %d.png', learning_rate, numLayer, epochs);
    saveas(gcf, filename);
end