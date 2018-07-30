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
