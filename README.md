### 51. Monte Carlo method [Wikipedia](https://en.wikipedia.org/wiki/Monte_Carlo_method)

Monte Carlo method to calculate value of Pi using ratio of square and circle.  
Plotting rectangle, circle and points up to n = 50 000.
```mathematica
getPiCircle[points_, n_] := N[4*Mean[If[Norm[#] < 1, 1, 0] & /@ points], 5];
gPlot2D[points_, n_] := Module[{},
   pi = getPiCircle[points, n];
   gRectangle = {White, EdgeForm[Black], Rectangle[{-1, -1}, {1, 1}]};
   gCircle =  Circle[];
   gPoints = {If[Norm[#] < 1, Green, Red], PointSize[.002], Point[#]} & /@ points;
   gPi = Text[StringForm["Pi = ``", pi], {0, 1.1}];
   gN = Text[StringForm["n = ``", n], {0, 1.2}];
   Show[Graphics[{gRectangle, gCircle, gPoints, gPi, gN}], ImageSize -> Large]
];
n = 50000;
points = RandomReal[{-1, 1}, {n, 2}];
frames = Table[gPlot2D[Take[points, {1, i}], i], {i, 500, 50000, 500}];
ListAnimate@frames;
```
![Image of Monte Carlo 2D](/img/51.1.MonteCarlo.gif)

Second method calculates Pi in 3D coordinates using ratio of cube and sphere.  
Plotting cube, sphere and points up to n = 50 000.
```mathematica
getPiSphere[points_, n_] := N[6*Mean[If[Norm[#] < 1, 1, 0] & /@ points], 5];
gPlot3D[points_, n_] := Module[{},
   pi = getPiSphere[points, n];
   gSphere =  {Opacity[.1], Sphere[{0, 0, 0}]};
   gPoints = {If[Norm[#] < 1, Green, Red], PointSize[.004], Point[#]} & /@ points;
   gPi = "Pi = " <> ToString@pi;
   gN = "n = " <> ToString@n;
   Show[Graphics3D[{gSphere, gPoints}], ImageSize -> Large, Epilog -> {Inset[gN, {.02, 0.99}], Inset[gPi, {.02, 0.96}]}]
];
n = 50000;
points = RandomReal[{-1, 1}, {n, 3}];
frames = Table[gPlot3D[Take[points, {1, i}], i], {i, 500, 50000, 500}];
ListAnimate@frames;
```
![Image of Monte Carlo 3D](/img/51.2.MonteCarlo.gif)

### 52. Lagrange interpolation [Wikipedia](https://en.wikipedia.org/wiki/Lagrange_polynomial)

Defining Lagrange interpolation function. Built-in function [Interpolation](http://reference.wolfram.com/language/ref/Interpolation.html) also implements Lagrange polynomials.  
```mathematica
interpolate[points_] := Module[{y = 0},
   {xvalues, yvalues} = points;
   Table[{
     Table[
      If[i != j, lx *= (x - xvalues[[i]])/(xvalues[[j]] - xvalues[[i]])], {i, 1, Length@xvalues}];
      y += yvalues[[j]]*lx;
      lx = 1;
     }, {j, 1, Length@xvalues}];
   Simplify@y
   ];
```

Creating points on Sine wave.
```mathematica
points = Table[{i, 5*N@Sin[i]}, {i, 0, 32, 1}];
```

Drawing points, interpolating Sine function and plotting result.
```mathematica
gFunction[points_] := Plot[interpolate[Transpose[points]], {x, 0, First@Last@points}, AspectRatio -> Automatic, PlotRange -> {{0, 10*Pi}, {-6, 6}}, Ticks -> {Table[t, {t, 0, 10*Pi, Pi}], {-10, -5, 0, 5, 10}}];
gPoints[points_] := ListPlot[points, PlotStyle -> Red];
gPlot[points_] := Show[gFunction[points], gPoints[points], ImageSize -> Large];

frames = Table[gPlot[Take[points, {1, i}]], {i, 2, Length@points}];
ListAnimate[frames];
```

![Image of Lagrange interpolation](/img/52.LagrangeInterpolation.gif)

### 53. Least Squares [Wikipedia](https://en.wikipedia.org/wiki/Least_squares) [Linear Regression](https://en.wikipedia.org/wiki/Linear_regression) [MIT 18.02SC](https://www.youtube.com/watch?v=YwZYSTQs-Hk)

Simulating data.
```mathematica
data = {#, 2 + 0.5*# + RandomReal[{-5, 5}]} & /@ RandomReal[{-10, 10}, 100];
```

Defining Least squares function explicitly. Built-in function [LeastSquares](http://reference.wolfram.com/language/ref/LeastSquares.html).  
Plotting result
```mathematica
draw[data_] := Module[{},
   {xi, yi} = Transpose[data];
   {a1, b1} = {a, b} /. First@Solve[{a*Total[xi^2] + b*Total[xi] == Total[xi*yi], a*Total[xi] + b*Length[xi] == Total[yi]}, {a, b}];
   eq = b1 + a1*x;
   gN = Graphics@Text[StringForm["n = ``", Length@data], {-8.5, 9.5}];
   gA = Graphics@Text[StringForm["a = ``", a1], {-8.5, 9}];
   gB = Graphics@Text[StringForm["b = ``", b1], {-8.5, 8.5}];
   plot1 = ListPlot[Transpose[{xi, yi}], PlotStyle -> Blue, AspectRatio -> Automatic, PlotRange -> {{-10, 10}, {-10, 10}}];
   plot2 = Plot[eq, {x, -10, 10}, PlotStyle -> {Red, Thin}];
   Show[plot1, plot2, gN, gA, gB, ImageSize -> Large]
   ];
```

Folding points and creating frames.
```mathematica
points = Table[Take[data, {1, i}], {i, 2, Length@data}];
frames = Table[draw[points[[i]]], {i, 1, Length@points}];
ListAnimate@frames;
```

![Image of Least Squares](/img/53.LeastSquares.gif)

### 54. k-NN classification [Wikipedia](https://en.wikipedia.org/wiki/K-nearest_neighbors_algorithm) [GeeksforGeeks](https://www.youtube.com/watch?v=odqIu23OSbs)

Simulating data. Creating 10 green points, 10 blue points and 50 unclassified points.
```mathematica
data = Join[{#, Blue} & /@ RandomReal[{-5, 1}, {10, 2}], {#, Green} & /@ RandomReal[{-1, 5}, {10, 2}]];
unclassified = RandomReal[{-5, 5}, {50, 2}];
```

Calculating Euclidean distance between given point and data. Classifying given point with most common color.  
Explicitly defined, built-in methods [NearestNeighbors](http://reference.wolfram.com/language/ref/method/NearestNeighbors.html)
```mathematica
distance[p1_, p2_] := N@Norm[(p1 - p2)^2];
classify[point_, data_, k_] := Module[{},
   sorted = SortBy[{distance[point, #[[1]]], #[[2]]} & /@ data, First];
   First@First@GatherBy[Take[sorted[[All, 2]], {1, k}]]
   ];
```

Visualizing data.
```mathematica
visualizePoints[data_] := Table[{data[[i, 2]], PointSize[.02], Point[data[[i, 1]]]}, {i, 1, Length@data}];
draw[data_] := Show[VoronoiMesh[data[[All, 1]], {{-5, 5}, {-5, 5}}, PlotTheme -> "Monochrome"], Graphics[visualizePoints[data]], ImageSize -> Large];
```

Classifying data, k=1.
```mathematica
k = 1;
classified = Table[{unclassified[[i]], classify[unclassified[[i]], data, k]}, {i,1, Length@unclassified}];
```

![Image of kNN](/img/54.kNN.gif)

### 55. Naive Bayes classifier [Wikipedia](https://en.wikipedia.org/wiki/Naive_Bayes_classifier) [Saed Sayad](http://www.saedsayad.com/naive_bayesian.htm)

Simulating data. Creating 10 of each red/green/blue points on (x,y) plane. Distributed according to normal distribution.
```mathematica
createData[mu_, class_, n_] :=  Table[{RandomVariate[MultinormalDistribution[mu, IdentityMatrix[2]]], class}, {i, 1, n}];
data = {createData[{0, -2}, Red, 10], createData[{-2, 2}, Green, 10], createData[{2, 2}, Blue, 10]};
unclassified = RandomReal[{-4, 4}, {100, 2}];
```

Defining Normal distrubution PDF formula, we use it to calculate probability.
```mathematica
normalDist[t_, std_, mean_] := 1/(Sqrt[2*Pi]*std)*E^-((t - mean)^2/(2*std^2));
probability[point_, class_] := normalDist[First@point, StandardDeviation[class[[All, 1]]], Mean[class[[All, 1]]]]*normalDist[Last@point, StandardDeviation[class[[All, 2]]], Mean[class[[All, 2]]]];
```

Classifying data, sorting results according to probabilities and choosing most popular result.  
Built-in method [NaiveBayes](http://reference.wolfram.com/language/ref/method/NaiveBayes.html)
```mathematica
classify[data_, point_] := Module[{},
   {red, green, blue} = data;
   pRed = {Red, probability[point, red[[All, 1]]]};
   pGreen = {Green, probability[point, green[[All, 1]]]};
   pBlue = {Blue, probability[point, blue[[All, 1]]]};
   color = First@Last@SortBy[{pRed, pGreen, pBlue}, Last];
   {point, color}
   ];
classified = Table[classify[data, unclassified[[i]]], {i, 1, Length@unclassified}];
```

![Image of NaiveBayes](/img/55.NaiveBayes.gif)

### 56. Gradient Descent [Wikipedia](https://en.wikipedia.org/wiki/Gradient_descent)

Simulating data. Generating points with normally distributed noise. Original slope a = 3 and intercept b = -10.
```matlab
function [X, Y] = generate_data(N)
    % Generate 2D linear dataset
    X = 2 * 2 * randn(N, 1);
    Y = 3 * X - 10 + randn(N, 1);
end
```

Single step calculates and minimizes error function
```matlab
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
```

Iterating until max_step is reached. Saving slope and intercept for animation.
```matlab
function [history_slope, history_intercept, history_cost] = gradient_descent(X, Y, max_steps, rate)
    % Calculating slope, intercept and cost on each step
    for i=1:max_steps
        history_slope(i,:) = current_slope;
        history_intercept(i,:) = current_intercept;
        history_cost(i,:) = sum((Y - (current_slope * X + current_intercept)).^2) / size(Y, 1);
        [current_slope, current_intercept] = gd_step(X, Y, current_slope, current_intercept, rate);
    end
end
```

![Image of GradientDescent](/img/56.GradientDescent.gif)

### 57. Stepwise Regression [Wikipedia](https://en.wikipedia.org/wiki/Stepwise_regression)

Simulating data. Generating 3D points with significant noise. Original slope a1 = 1 and a2 = 5.
```matlab
function [X, Y] = generate_data(N)
    % Generate 3D linear dataset
    X1 = 0.5 * randn(N, 1);
    X2 = 0.5 * randn(N, 1);
    Y = 1 * X1 + 5 * X2 + randn(N, 1);
    % Add significant noise
    X1 = [X1; -10 + 20*rand(20,1)];
    X2 = [X2; -10 + 20*rand(20,1)];
    Y = [Y; -10 + 20*rand(20,1)];
end
```

Defining stepwise regression running function. Running F-test on data to remove points with significant variance. Using mean squares method to calculate parameters. Using R^2 to validate model and if model not good enough then going for next iteration.
```matlab
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
```

Running implemented stepwise model builder with initial parameters: max_iterations = 20, minimum R^2 value = 0.96 and alpha = 0.02.
```matlab
[R_squared_performance, ~, ~, slope, intercept, history_slope, history_intercept] = stepwise_regression(X, Y, 20, 0.95, 0.02);
```

![Image of StepwiseRegression](/img/57.StepwiseRegression.gif)

### 58. Gradient boosting [Wikipedia](https://en.wikipedia.org/wiki/Gradient_boosting) [Tutorial]https://medium.com/mlreview/gradient-boosting-from-scratch-1e317ae4587d

Defining function for gradient boosting. Using decision stumps as weak models.
```matlab
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
```

Creating data and running model for N epochs.
```matlab
[X, Y] = create_data();
for i=1:N
    [error_history, Y_hat] = gradient_boosting(X, Y, i);
end
```

![Image of GradientBoosting](/img/58.GradientBoosting.gif)
