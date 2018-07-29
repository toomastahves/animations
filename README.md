### 51. Monte Carlo method [Wikipedia](https://en.wikipedia.org/wiki/Monte_Carlo_method)

Defining amount of points and creating random points between {-1,1} on {x,y} plane.
```mathematica
n = 50000;
points = RandomReal[{-1, 1}, {n, 2}];
```

Calculating potential Pi value and drawing plot.
```mathematica
gPlot[points_, n_] := Module[{},
   pi = N[4*Count[If[#[[1]]^2 + #[[2]]^2 < 1, True, False] & /@ points, True]/n, 5];
   gRectangle = {White, EdgeForm[Black], Rectangle[{-1, -1}, {1, 1}]};
   gCircle =  Circle[];
   gPoints = {If[#[[1]]^2 + #[[2]]^2 < 1, Green, Red], PointSize[.001], Point[#]} & /@ points;
   gN = Text[StringForm["Pi = ``", pi], {0, 1.1}];
   gPi = Text[StringForm["n = ``", n - 1], {0, 1.2}];
   Show[Graphics[{gRectangle, gCircle, gPoints, gN, gPi}], ImageSize -> Large]
   ];
```

Animating.
```mathematica
folded = Table[{n, Take[points, {1, n}]}, {n, 1, 50000, 500}];
frames = Table[gPlot[folded[[n, 2]], folded[[n, 1]]], {n, 1, 100}];
ListAnimate@frames;
```

![Image of Monte Carlo](/img/51.MonteCarlo.gif)