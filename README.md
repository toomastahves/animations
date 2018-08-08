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

Visualizing data. Using VoroinoiMesh to visualize boundaries.
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
