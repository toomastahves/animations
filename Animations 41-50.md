
### 41. Conway's Game of Life [Wikipedia](https://en.wikipedia.org/wiki/Conway%27s_Game_of_Life)

Defining function to check 8 neighbors of the cell.
```mathematica
isValid[x_, y_, n_] := If[(x < 1 || y < 1 || x > n || y > n), False, True];
checkNeighbors[matrix_, x_, y_] := Module[{},
   n = Length@matrix;
   up = If[isValid[x, y + 1, n], matrix[[x, y + 1]], 0];
   down = If[isValid[x, y - 1, n], matrix[[x, y - 1]], 0];
   right = If[isValid[x + 1, y, n], matrix[[x + 1, y]], 0];
   rightup = If[isValid[x + 1, y + 1, n], matrix[[x + 1, y + 1]], 0];
   rightdown = If[isValid[x + 1, y - 1, n], matrix[[x + 1, y - 1]], 0];
   left = If[isValid[x - 1, y, n], matrix[[x - 1, y]], 0];
   leftup = If[isValid[x - 1, y + 1, n], matrix[[x - 1, y + 1]], 0];
   leftdown = If[isValid[x - 1, y - 1, n], matrix[[x - 1, y - 1]], 0];
   Total[{up, down, right, rightup, rightdown, left, leftup, leftdown}]
   ];
```

Defining GameOfLife function and rules how to act depending on neighbors.  
Rewrites entire table. Faster solution would be to update only part that has changed.
```mathematica
GameOfLife[matrix_] := Module[{},
   neighbors = Flatten[Table[{checkNeighbors[matrix, i, j], matrix[[i, j]], {i, j}}, {i, 1, Length@matrix}, {j, 1, Length@matrix}], 1];
   updates = Table[{
      Switch[neighbors[[i, 1]],
       0, neighbors[[i, 3]] -> If[neighbors[[i, 2]] == 0, 0, 0],
       1, neighbors[[i, 3]] -> If[neighbors[[i, 2]] == 0, 0, 0],
       2, neighbors[[i, 3]] -> If[neighbors[[i, 2]] == 0, 0, 1],
       3, neighbors[[i, 3]] -> If[neighbors[[i, 2]] == 0, 1, 1],
       4, neighbors[[i, 3]] -> If[neighbors[[i, 2]] == 0, 0, 0],
       5, neighbors[[i, 3]] -> If[neighbors[[i, 2]] == 0, 0, 0],
       6, neighbors[[i, 3]] -> If[neighbors[[i, 2]] == 0, 0, 0],
       7, neighbors[[i, 3]] -> If[neighbors[[i, 2]] == 0, 0, 0],
       8, neighbors[[i, 3]] -> If[neighbors[[i, 2]] == 0, 0, 0]
       ]
      },
     {i, 1, Length@neighbors}];
   updates = Flatten@updates;
   ReplacePart[matrix, updates]
   ];
```

Creating matrix, setting up initial position and running game.  
On the picture glider crashes into block, destroys beehive+blinker and ends with 4 blinkers.
```mathematica
matrix = createMatrix[20];
seed = {
   {5, 14}, {6, 14}, {5, 16}, {6, 16}, {4, 15}, {7, 15},
   {10, 10}, {11, 11}, {10, 11}, {11, 10},
   {15, 11}, {15, 12}, {15, 13},
   {4, 5}, {5, 6}, {5, 7}, {4, 7}, {3, 7}
   };
matrix = ReplacePart[matrix, seed -> 1];
result  = NestList[GameOfLife, matrix, 100];
Animate[Show[draw[result[[i]]]], {i, 1, Length@result, 1}, AnimationRunning -> False];
```

![Image of Game of Life](/img/41.GameOfLife.gif)

### 42. Mandelbrot set [Wikipedia](https://en.wikipedia.org/wiki/Mandelbrot_set) [Rosettacode](http://rosettacode.org/wiki/Mandelbrot_set)

Calculating Mandelbrot set.  
Input c = y + x*I. Formula z = z^2 + c. Iteration count n. If Abs[z] > 2, then return current iteration number.
```mathematica
MandelBrotSet = With[
   {color = Module[{c = #1, z = #1, n = #2},
       Do[If[Abs[z] > 2, Return@i]; z = z^2 + c, {i, n}]
       ] &
    },
   Compile[{{n, _Integer}},
    Table[color[y + x*I, n], {x, -2, 2, 0.005}, {y, -2, 2, 0.005}]
    ]
   ];
```

Plotting result upto n = 30.
```mathematica
frames = ParallelTable[ArrayPlot[MandelBrotSet[i], ImageSize -> Large], {i, 0, 30}];
ListAnimate[frames];
```

![Image of Mandelbrot set](/img/42.MandelbrotSet.gif)

### 43. Julia set [Wikipedia](https://en.wikipedia.org/wiki/Julia_set) [Stackexchange](https://mathematica.stackexchange.com/questions/21714/why-is-this-mandelbrot-sets-implementation-infeasible-takes-a-massive-amount-o)

Calculating Julia set.
```mathematica
JuliaSet = With[
   {color = Module[{z = #1, c = #2, n = #3},
       Do[If[Abs[z = z^2 + c] > 2, Return[i]], {i, n}]
       ] &
    },
   Compile[{{c, _Complex}, {n, _Integer}},
    Table[color[y + x*I, c, n], {x, -2, 2, 0.005}, {y, -2, 2, 0.005}]
    ]
   ];
```

Plotting result for c = 0.7885 * E^(a * I), where a=[0,2Pi] and n = 50.
```mathematica
frames = ParallelTable[ArrayPlot[JuliaSet[0.7885*E^(a*I), 50], ImageSize -> Large, ColorFunction -> "TemperatureMap"], {a, Pi/90, 2*Pi, Pi/90}];
ListAnimate@frames;
```

![Image of Julia set](/img/43.JuliaSet.gif)

### 44. Pentaflake [N-flake](https://en.wikipedia.org/wiki/N-flake) [Golden ratio](https://en.wikipedia.org/wiki/Golden_ratio) [Pentaflake](http://mathworld.wolfram.com/Pentaflake.html) 

First version of pentaflake.  
1. Defining 5 points on inner circle of pentagon (in the middle of the each edge of pentagon). Using GoldenRatio to calculate next radius.  
2. Rotating original pentagon around that point by Pi degree. Making 5 rotations, thus creating 5 new pentagons.  
3. Flipping result by Pi degrees, not necessary for creating pentaflake, but visually makes image stable.
```mathematica
PentaFlake1[{polygon0_, n_: 1}] := Module[{},
   points = CirclePoints[{0.5*GoldenRatio^(2*n - 1), 3*Pi/10}, 5];
   rotatedPolygons = GeometricTransformation[polygon0, RotationTransform[Pi, #] & /@ points];
   flippedPolygon = GeometricTransformation[{polygon0, rotatedPolygons}, RotationTransform[Pi, {0, 0}]];
   {flippedPolygon, n + 1}
   ];
```
![Image of Pentaflake1](/img/44.1.Pentaflake.gif)

Second version of pentaflake.
1. Defining 5 points on inner circle of pentagon (in the middle of the each edge of pentagon).  
2. Translating pentagon towards each point. Direction vector is built using pentagon middle point {0,0} and previously created points.
3. Flipping for visual purpose not necessary, looks good.
```mathematica
PentaFlake2[{polygon0_, n_: 1}] := Module[{},
   points = CirclePoints[{GoldenRatio^(2*n - 1), Pi/10}, 5];
   polygon1 = GeometricTransformation[polygon0, TranslationTransform /@ points];
   {polygon1, n + 1}
   ];
```
![Image of Pentaflake2](/img/44.2.Pentaflake.gif)

### 45. Menger Sponge [Wikipedia](https://en.wikipedia.org/wiki/Menger_sponge)

Input is array of cube coordinates.  
Splitting each cube into 27 smaller cubes and removing selected 7 cubes.
```mathematica
MengerSponge[cubes_] := Module[{newcubes},
   c = Table[{
      newcuboid = cubes[[i]];
      length = N@Abs[newcuboid[[2, 1]] - newcuboid[[1, 1]]]/3;
      newcubes = Table[{
         {x, y, z}, {x + length, y + length, z + length}},
        {x, newcuboid[[1, 1]], newcuboid[[1, 1]] + 2*length, length},
        {y, newcuboid[[1, 2]], newcuboid[[1, 2]] + 2*length, length},
        {z, newcuboid[[1, 3]], newcuboid[[1, 3]] + 2*length, length}];
      newcubes = Flatten[newcubes, 2];
      newcubes = Delete[newcubes, {{5}, {11}, {13}, {14}, {15}, {17}, {23}}];
      newcubes
      }, {i, 1, Length@cubes}];
   c = Flatten[c, 2]
   ];
```

Defining first cube, splitting n = 4 and animating.
```mathematica
init = {{{0, 0, 0}, {3, 3, 3}}};
result = NestList[MengerSponge, init, 4];
frames = Table[{sponge = result[[i]]; draw = Table[Graphics3D@{Yellow, EdgeForm[Black], Cuboid[sponge[[i, 1]], sponge[[i, 2]]]}, {i, 1, Length@sponge}]; Show[draw, ImageSize -> Large]}, {i, 1, Length@result}];
ListAnimate@frames;
```

![Image of MengerSponge](/img/45.MengerSponge.gif)

### 46. Barnsley fern [Wikipedia](https://en.wikipedia.org/wiki/Barnsley_fern)

Function that returns random point according to rules.
```mathematica
BarnsleyFern[{x_, y_}] := Module[{},
   i = RandomInteger[{1, 100}];
   If[i <= 1, {xt = 0, yt = 0.16*y},
    If[i <= 8, {xt = 0.2*x - 0.26*y, yt = 0.23*x + 0.22*y + 1.6},
     If[i <= 15, {xt = -0.15*x + 0.28*y, yt = 0.26*x + 0.24*y + 0.44},
      {xt = 0.85*x + 0.04*y, yt = -0.04*x + 0.85*y + 1.6}]]];
   {xt, yt}];
```

Animating upto 200 000 points.
```mathematica
points = NestList[BarnsleyFern, {0, 0}, 200000];
gPoints = {Hue[.35, 1, .7], PointSize[.001], Point[#]} & /@ points;
partitioned = Partition[gPoints, 1000];
folded = FoldList[List, partitioned];
frames = Show[Graphics[#], ImageSize -> Large] & /@ folded;
ListAnimate@frames;
```

![Image of Barnsley fern](/img/46.BarnsleyFern.gif)

### 47. Tesseract [Wikipedia](https://en.wikipedia.org/wiki/Tesseract) [Math Stackexchange](https://mathematica.stackexchange.com/questions/9580/how-to-create-this-four-dimensional-cube-animation) [4D Visualization](http://eusebeia.dyndns.org/4d/vis/vis)

Defining vertices (points), edges (lines) and faces (planes) for tesseract.
```mathematica
vertices = Tuples[{-1, 1}, 4];
edges = Select[Subsets[Range[Length[vertices]], {2}], Count[Subtract @@ vertices[[#]], 0] == 3 &];
faces = Select[Union[Flatten[#]] & /@ Subsets[edges, {4}], Length@# == 4 &] /. {a_, b_, c_, d_} -> {a, b, d, c};
```

Defining rotation function. Double rotation in two hyperplanes. Calculating for vertices.  
Projection function projects 4D object into 3D coordinates, so we see "3D shadow" of 4D object.
```mathematica
rotation[t_] := RotationMatrix[t, {{0, 0, 1, 0}, {0, 1, 0, 0}}].RotationMatrix[t, {{1, 0, 0, 0}, {0, 0, 0, 1}}].# & /@ vertices;
projection[t_] := Most[#]/(3 - Last[#]) & /@ rotation[t];
```

Building visual components for tesseract.
```mathematica
tesseract[t_] := GraphicsComplex[projection[t], {
   {Red, Sphere[Range[16], 0.01]},
   {Black, Line[edges]},
   {Hue[.7, 1, .8, .2], Polygon[faces]}
   }]
Animate[Show[Graphics3D[tesseract[t], Boxed -> False], PlotRange -> 1], {t, Pi/90, Pi/2}]
```

![Image of Tesseract](/img/47.Tesseract.gif)

### 48. GeoGraphics [Mapping GPS Data](http://blog.wolfram.com/2009/04/17/mapping-gps-data/) [GPS Mountainbike analysis](http://community.wolfram.com/groups/-/m/t/567757)

Defining function for parsing data. Input is location of .gpx file.  
Parsing out {latitude, longitude} points, calculating distance (meters) and time (seconds), scaling so they start from 0.
```mathematica
getData[filename_] := Module[{},
   SetDirectory[NotebookDirectory[]];
   xml = Import[filename, "XML"];
   data = Cases[xml, XMLElement["trkpt", _, _], All];
   points = {ToExpression["lat" /. #[[2]]], ToExpression["lon" /. #[[2]]]} & /@ data;
   distance = FoldList[Plus, 0, GeoDistance @@ # & /@ Partition[points, 2, 1]];
   time = Map[FirstCase[#, XMLElement["time", {}, {time_}] :> AbsoluteTime[time] + $TimeZone*3600, Missing[], \[Infinity]] &, data]; time -= First[time];
   {points, distance, time}
   ];
```

Defining plotting function.
1. Reading in data.  
2. Defining map bounds by fining Min and Max values for latitude and longitude.  
3. Multiple filtering methods for interactive animation.  
4. Using ListLinePlot to draw plot.  
5. Drawing GeoPath from points, using Epilog to insert plot.
```mathematica
gPlot[{paths_, distances_, times_}, i_] := Module[{},
   {p1, p2} = paths; {d1, d2} = distances/1000; {t1, t2} = times/60;
   range = {{Min@paths[[All, All, 1]], Max@paths[[All, All, 1]]}, {Min@paths[[All, All, 2]], Max@paths[[All, All, 2]]}};
   
   d1f = If[i <= Length@d1, Take[d1, {1, i}], d1];
   d2f = If[i <= Length@d2, Take[d2, {1, i}], d2];
   t1f = If[i <= Length@t1, Take[t1, {1, i}], t1];
   t2f = If[i <= Length@t2, Take[t2, {1, i}], t2];
   p1f = If[i <= Length@p1, Take[p1, {1, i}], p1];
   p2f = If[i <= Length@p2, Take[p2, {1, i}], p2];
   p1pos = If[i <= Length@p1, p1[[i]], Last@p1];
   p2pos = If[i <= Length@p2, p2[[i]], Last@p2];
   
   stats = ListLinePlot[
     {Transpose[{d1f, t1f}], Transpose[{d2f, t2f}]},
     AxesLabel -> {"(km)", "(min)"}, ImageSize -> 300, 
     PlotStyle -> {Red, Blue}, PlotRange -> {{0, 10}, {0, 120}}
     ];
   
    GeoGraphics[
    {Red, GeoPath[p1f], PointSize[.02], Point[GeoPosition[p1pos]],
     Blue, GeoPath[p2f], PointSize[.02], Point[GeoPosition[p2pos]]
     },
    GeoRange -> range, GeoRangePadding -> 100, ImageSize -> Large, 
    Epilog -> {Inset[stats, {Right, Bottom}, {Right, Bottom}]}
    ]
   ];
```

Using functions and animating.
```mathematica
{p1, d1, t1} = getData["48.path1.gpx"];
{p2, d2, t2} = getData["48.path2.gpx"];
data = {{p1, p2}, {d1, d2}, {t1, t2}};
frames = Table[gPlot[data, i], {i, 2, 1200, 10}];
ListAnimate@frames;
```

![Image of GeoGraphics](/img/48.GeoGraphics.gif)

### 49. Rolling shutter [Wikipedia](https://en.wikipedia.org/wiki/Rolling_shutter) 

Defining propeller (Rose curve), shutter (line) and plotting functions.
```mathematica
propeller[t_, a_] := {Cos[5*t]*Cos[(t + a/5)], Cos[5*t]*Sin[(t + a/5)]};
shutter[s_, a_] := {s, a/Pi};
plotPropeller[a_] := ParametricPlot[propeller[t, a], {t, -Pi, Pi}, PlotRange -> 1, PlotStyle -> Black, Axes -> False];
plotShutter[a_] := ParametricPlot[shutter[s, a], {s, -Pi, Pi}, PlotStyle -> Blue];
```

Defining function to calculate intersection points between propeller and shutter.
```mathematica
calculateIntersection[a_] := Module[{t, s},
   solve = NSolve[{propeller[t, a][[1]] == shutter[s, a][[1]], propeller[t, a][[2]] == shutter[s, a][[2]]}, {t, s}, Reals] /. C[1] -> 1;
   coordinates = If[Length@solve != 0, shutter[#, a] & /@ ({t, s} /. solve)[[All, 2]]];
   If[Length@coordinates != 0, Graphics[{PointSize[0.005], Red, Point[#]}] & /@ coordinates, Graphics[{}]]
   ];
```

Calculating 90 lines where intersection points are laid.  
Using Association (key-value map) to hold rotation angle+points. Using Take to fold points, so they appear on screen after line has passed.
```mathematica
n = 90;
points = Table[{a, calculateIntersection[a]}, {a, -Pi, Pi, Pi/n}];
assoc = Association[Table[points[[i, 1]] -> Take[points[[All, 2]], {1, i}], {i, 1, Length@points}]];

frames = Table[Show[plotPropeller[a], plotShutter[a], assoc[a], ImageSize -> Large], {a, -Pi, Pi, Pi/n}];
ListAnimate@frames;
```

![Image of Rolling shutter](/img/49.RollingShutter.gif)

### 50. Bessel function [Vibrations of a circular membrane](https://en.wikipedia.org/wiki/Vibrations_of_a_circular_membrane) 

Defining 3D membrane, plotting function, projections and animating.
```mathematica
f[n_, k_, time_] := {s* Cos[t], s* Sin[t], Cos[time]*Cos[n*t]*N@BesselJ[n, s*N@BesselJZero[n, k]]};
gDrum[n_, k_, time_] := ParametricPlot3D[Evaluate@{f[n, k, time], ReplacePart[f[n, k, time], 2 -> 3], ReplacePart[f[n, k, time], 1 -> -3]}, {s, 0, 2}, {t, 0, 2*Pi}, PlotRange -> {{-3, 3}, {-3, 3}, {-2, 2}}, Mesh -> True, Boxed -> True, Axes -> True, PlotStyle -> Lighter@Pink];

frames = Table[Show[gDrum[3, 1, time], ImageSize -> Large], {time, -Pi, Pi, Pi/45}];
ListAnimate@frames;
```

![Image of Bessel function](/img/50.BesselFunction.gif)