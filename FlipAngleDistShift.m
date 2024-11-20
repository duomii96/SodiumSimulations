% Parameters for the normal distribution
mu = 90;    % Mean of the distribution
sigma = 1;  % Standard deviation of the distribution

% Generate 1000 samples from the normal distribution
num_samples = 1000;
normal_samples = normrnd(mu, sigma, num_samples, 1);

% Plot the histogram
figure;
histogram(normal_samples, 30, 'Normalization', 'pdf', 'FaceColor', 'blue');
title('Histogram of Normal Distribution');
xlabel('Value');
ylabel('Density');

% Overlay the theoretical normal distribution curve
hold on;
x = linspace(mu - 4*sigma, mu + 4*sigma, 100); % Values for the x-axis
pdf = normpdf(x, mu, sigma);  % Calculate the normal PDF
plot(x, pdf, 'r-', 'LineWidth', 2); % Plot the PDF
legend('Sampled Data', 'Theoretical PDF');
hold off;
