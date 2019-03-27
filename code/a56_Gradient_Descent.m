close all
clear

% Generate test data
[X, Y] = generate_data(200);
iterations = 200;
rate = 0.01;
% Calculate parameters using gradient descent
[history_slope, history_intercept, history_cost] = gradient_descent(X, Y, iterations, rate);

% Plot scattered points
fig = figure;
scatter(X, Y);
title('Gradient descent line search');
xlabel('X');
ylabel('Y');
axis tight manual;
xlim([-10, 10]);
ylim([-40, 20]);
grid;
hold on;

% Plot real line
x = linspace(-10, 10)';
y2 = 3 * x - 10;
l2 = line(x, y2);
l2.Color = [1,0,0];
hold on;

% Animation
N = size(history_slope, 1);
for i = 1:N
    % Plot current line
    x = linspace(-10, 10)';
    y = history_slope(i) * x + history_intercept(i);
    l1 = line(x, y);
    xlim([-10, 10]);
    ylim([-40, 20]);
    l1.Color = [0,0,0,0.2];
    a_text = strcat('a = ', ' ', num2str(round(history_slope(i),2)));
    b_text = strcat('b = ', ' ', num2str(round(history_intercept(i),2)));
    title_text = {a_text, b_text};
    title(title_text);
    hold on;
    
    % Animation part
    frame = getframe(fig);
	im{i} = frame2im(frame);
end

% Create gif
filename = 'gradient_descent.gif';
gif_speed = 0.05;
for idx = 1:N
    [A,map] = rgb2ind(im{idx},256);
    if idx == 1
        % Create initial gif file
        imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',gif_speed);
    else
        % Append frames to gif
        imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',gif_speed);
    end
end

% Generates linear data without significant noise
function [X, Y] = generate_data(N)
    % Generate 2D linear dataset
    X = 2 * 2 * randn(N, 1);
    Y = 3 * X - 10 + randn(N, 1);
end

% Gradient descent method for calculating slope and intercept of line
function [history_slope, history_intercept, history_cost] = gradient_descent(X, Y, max_steps, rate)
    current_slope = 0;
    current_intercept = 0;
    history_slope = zeros(max_steps, 1);
    history_intercept = zeros(max_steps, 1);
    history_cost = zeros(max_steps, 1);
    % Calculating slope, intercept and cost on each step
    for i=1:max_steps
        history_slope(i,:) = current_slope;
        history_intercept(i,:) = current_intercept;
        history_cost(i,:) = sum((Y - (current_slope * X + current_intercept)).^2) / size(Y, 1);
        [current_slope, current_intercept] = gd_step(X, Y, current_slope, current_intercept, rate);
    end
end

% Single step to calculate slope and intercept
function [slope, intercept] = gd_step(X, Y, slope, intercept, rate)
    N = size(Y, 1);
    % Y-hat = predicted Y
    Y_hat = slope * X + intercept;
    % Error function to be minimized
    error = (Y - Y_hat);
    % Derivatives
    D_slope = (-2/N) * sum(X .* error);
    D_intercept = (-2/N) * sum(error);
    % Calculating new values
    slope = slope - rate * D_slope;
    intercept = intercept - rate * D_intercept;
end

