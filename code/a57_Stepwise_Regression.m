close all
clear
%rng default

% Ex 3
[X, Y] = generate_data(100);
% Calculations
[R_squared_performance, ~, ~, slope, intercept, history_slope, history_intercept] = stepwise_regression(X, Y, 20, 0.95, 0.02);

% Plotting points
fig = figure;
scatter3(X(:,1), X(:,2), Y);
axis equal
xlabel('X1');
ylabel('X2');
zlabel('Y');
set(gca,'XLim',[-10, 10],'YLim',[-10, 10],'ZLim',[-10, 10]);
title('Stepwise regression');
campos([ -166.2451  -39.4836   28.3483])
hold on

% Plane mesh
[x, y] = meshgrid(-10:0.1:10);

% Animation
N = size(history_slope, 1);
for i = 1:N
    % Plot current line
    z = history_slope(i,1)*x + history_slope(i,2)*y + history_intercept(i,:);
    s = surf(x,y,z, 'FaceAlpha',0.3);
    s.EdgeColor = 'none';
    
    % Plot info
    r_squared = strcat('R^2 = ', ' ', num2str(round(R_squared_performance(i,2),2)));
    iter_text = strcat('Iterations = ', ' ', num2str(i));
    title_text = strcat(iter_text, " ", r_squared);
    title(title_text);
    
    % Animation part
    frame = getframe(fig);
	im{i} = frame2im(frame);
    if i < N
        set(s,'visible','off');
    end
end

% Create gif
filename = 'stepwise_regression.gif';
gif_speed = 0.2;
for idx = 1:N
    [A,map] = rgb2ind(im{idx},256);
    if idx == 1
        % Create initial gif file
        imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',gif_speed);
    elseif idx == N
        % Hold last frame
        imwrite(A,map,filename,'gif','WriteMode','append','DelayTime', 3);
    else
        % Append frames to gif
        imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',gif_speed);
    end
end

% Generates linear data with significant noise
function [X, Y] = generate_data(N)
    % Generate 3D linear dataset
    X1 = 0.5 * randn(N, 1);
    X2 = 0.5 * randn(N, 1);
    Y = 1 * X1 + 5 * X2 + randn(N, 1);
    % Add significant noise
    X1 = [X1; -10 + 20*rand(20,1)];
    X2 = [X2; -10 + 20*rand(20,1)];
    Y = [Y; -10 + 20*rand(20,1)];
    % Shuffle
    X = [X1, X2];
    new_idx = randperm(size(X, 1));
    X = X(new_idx,:);
    Y = Y(new_idx,:);
end

% Stepwise regression
function [R_squared_performance, X, Y, slope, intercept, slope_history, intercept_history] = stepwise_regression(X, Y, iterations, R_squared_limit, alpha)
    for i=1:iterations
        % F-test
        [idx] = f_test(X, Y, alpha);
        % Selecting passed test points
        X = X(idx == 0,:);
        Y = Y(idx == 0,:);
        % Calculating R^2
        [slope, intercept] = mean_squares(X, Y);
        R_squared = R_squared_custom(X, Y, slope, intercept);
        % Used for animation
        slope_history(i,:) = slope;
        intercept_history(i,:) = intercept;
        % Used for tracking learning rate
        R_squared_performance(i,:) = [i, R_squared];
        % Loop end condition
        if R_squared > R_squared_limit
            break
        end
    end
end

% Calculates coefficient of determination. Used for quality test.
function R_squared = R_squared_custom(X, Y, slope, intercept) 
    % Total sum of squares
    SS_tot = sum((Y - mean(Y)).^2);
    % Predicted Y value
    Y_hat = (intercept + sum(slope .* X', 1))';
    % Residual sum of squares
    SS_res = sum((Y - Y_hat).^2);
    % Regression sum of squares = explained sum of squares
    SS_reg = sum((Y_hat - mean(Y)).^2);
    % R^2
    R_squared = 1 - SS_res / SS_tot;
end

% Mean squares method
function [slope, intercept] = mean_squares(X, Y)
    slope = pinv(X' *  X) * X' * Y;
    intercept = sum(mean(Y) - slope .* mean(X)');
end

% Returns array, where 0 - passed test, 1 - failed test
function [idx] = f_test(X, Y, alpha)
    N = size(X, 1);
    % Iterating over N points and returning 'passed'/'not passed' for each point 
    idx = arrayfun(@(x) sum(get_Ftest(X, Y, x, alpha)), (1:N)');
end

% F-test implementation
function F_result = get_Ftest(X_combined, Y_combined, idx, alpha)
    % Remove i-th point from combined array
    X_simple = X_combined(setdiff(1:end, idx), :);
    Y_simple = Y_combined(setdiff(1:end, idx), :);
    % Calculate residual sum of squares
    RSS1 = get_RSS(X_simple, Y_simple);
    RSS2 = get_RSS(X_combined, Y_combined);
    % Calculate F score and compare to distribution F score
    F_calculated = ((RSS2 - RSS1)/2) / (RSS2/1);
    F_distribution = finv(alpha, 2, 1);
    % Return result, 0 - passed test, 1 - failed test
    F_result = F_calculated > F_distribution;
end

% Calculating residual sum of squares
function RSS = get_RSS(X, Y)
    [slope, intercept] = mean_squares(X, Y);
    Y_hat = (intercept + sum(slope .* X', 1))';
    RSS = sum((Y - Y_hat).^2);
end
