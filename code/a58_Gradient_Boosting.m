close all
clear

[X, Y] = create_data();
fig = figure;
scatter(X, Y);
set(gca,'XLim',[0, 150],'YLim',[0, 35]);
hold on;

N = 50;

for i=1:N
    [error_history, Y_hat] = gradient_boosting(X, Y, i);
    p = plot(X, Y_hat, 'r');
    title('Plotting data and boosted model');
    xlabel('X');
    ylabel('Y');
    grid on;
    
    % Plot info
    iter_text = strcat('Epochs = ', ' ', num2str(i));
    title_text = strcat(iter_text);
    title(title_text);
    
    % Animation part
    frame = getframe(fig);
	im{i} = frame2im(frame);
    if i < N
        set(p,'visible','off');
    end
end

% Create gif
filename = 'gradient_boosting.gif';
gif_speed = 0.1;
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

function [x, y] = create_data()
    x = (1:150)';
    
    y1 = randi([0,20],10,1) / 10;
    y2 = randi([60,110],10,1) / 10;
    y3 = randi([120,150],10,1) / 10;
    y4 = randi([60,110],10,1) / 10;
    y5 = randi([0,20],10,1) / 10;
    
    y6 = randi([120,150],10,1) / 10;
    y7 = randi([170,210],10,1) / 10;
    y8 = randi([290,300],10,1) / 10;
    y9 = randi([170,210],10,1) / 10;
    y10 = randi([120,150],10,1) / 10;
    
    y11 = randi([0,20],10,1) / 10;
    y12 = randi([60,110],10,1) / 10;
    y13 = randi([120,150],10,1) / 10;
    y14 = randi([60,110],10,1) / 10;
    y15 = randi([0,20],10,1) / 10;
    
    y = [y1; y2; y3; y4; y5; y6; y7; y8; y9; y10; y11; y12; y13; y14; y15];
end

% Gradient boosting algorithm
function [error_history, y_hat] = gradient_boosting(X, Y, epochs)
    Y_temp = Y;
    N = size(Y_temp, 1);
    y_hat = zeros(N, 1);
    error_history = zeros(epochs, 1);
    % Training model
    for i=1:epochs
        % Creating stump
        stump = fitrtree(X,Y_temp,'minparent',size(X,1),'prune','off','mergeleaves','off');
        % Splitting data in two parts
        cut_point = stump.CutPoint(1);
        left_idx = find(X <= cut_point);
        right_idx = find(X > cut_point);
        % Calculating residuals for each part
        residuals = zeros(N, 1);
        residuals(left_idx) = mean(Y_temp(left_idx));
        residuals(right_idx) = mean(Y_temp(right_idx));
        % Predicting new Y value
        y_hat = y_hat + residuals;
        error = Y - y_hat;
        Y_temp = error;
        % Saving error for plotting
        error_history(i,:) = sum(abs(error));
    end
end
