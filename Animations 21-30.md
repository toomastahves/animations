### 21. Lorenz attractor [Wikipedia](https://en.wikipedia.org/wiki/Lorenz_system)

Setting up curve parameters and equations.
```mathematica
sigma = 10; rho = 28; beta = 8/3;
equations = {
   x'[t] == sigma*(y[t] - x[t]), 
   y'[t] == x[t]*(rho - z[t]) - y[t], 
   z'[t] == x[t]*y[t] - beta*z[t],
   x[0] == 10,
   y[0] == 10,
   z[0] == 10
   };
```

Solving differential equations.
```mathematica
{px, py, pz} = NDSolveValue[equations, {x, y, z}, {t, 0, 50}];
```

Plotting points.
```mathematica
gPlot[u_] := ParametricPlot3D[{px[t], py[t], pz[t]}, {t, 0, u}];
```

![Image of Lorenz attractor](/img/21.LorenzAttractor.gif)

### 22. Klein Bottle [Wikipedia](https://en.wikipedia.org/wiki/Klein_bottle)

Defining Klein Bottle shape equation in parametric form.
```mathematica
kbx = -2/15*Cos[u]*(3*Cos[v] - 30 Sin[u] + 90*Cos[u]^4*Sin[u] - 60*Cos[u]^6*Sin[u] + 5*Cos[u]*Cos[v]*Sin[u]);
kby = -1/15*Sin[u]*(3*Cos[v] - 3*Cos[u]^2*Cos[v] - 48*Cos[u]^4*Cos[v] + 48*Cos[u]^6*Cos[v] - 60*Sin[u] + 5*Cos[u]*Cos[v]*Sin[u] - 5*Cos[u]^3*Cos[v]*Sin[u] - 80*Cos[u]^5*Cos[v]*Sin[u] + 80*Cos[u]^7*Cos[v]*Sin[u]);
kbz = 2/15*(3 + 5*Cos[u]*Sin[u])*Sin[v];
```

Plotting Bottle shape. 0 < u < Pi and 0 < v < 2Pi.
```mathematica
gKleinBottle[u1_, u2_, v1_, v2_] := ParametricPlot3D[{kbx, kby, kbz}, {u, u1, u2}, {v, v1, v2}, 
   Boxed -> False, Axes -> False, MeshStyle -> Directive[Red], PlotRange -> {{-1.6, 2}, {-0.2, 4.4}, {-0.8, 0.8}}, Mesh -> False, 
   ColorFunction -> Function[{x, y, z, t}, Hue[x/4 + 0.55, 0.7, 0.9, 0.5]]];
```

Plotting neck. 0 < v < 2Pi. Moving animation is created by changing parameter u: 0 < u < Pi
```mathematica
gNeck[t_] := ParametricPlot3D[{kbx /. u -> t, kby /. u -> t, kbz /. u -> t}, {v, 0, 2*Pi}, Boxed -> False, Axes -> False, 
  PlotRange -> {{-1.6, 2}, {-0.2, 4.4}, {-0.8, 0.8}}, Mesh -> False, PlotStyle -> Hue[0.9, 1, 1, 0.7]]
```

![Image of Klein Bottle](/img/22.KleinBottle.gif)

### 23. 3D Lissajous [MathCurve](https://www.mathcurve.com/courbes3d.gb/lissajous3d/noeudlissajous.shtml)

Defining 3D Lissajous curve and projections.
```mathematica
fxyz[a_, b_, c_, q_, p_, r_] := {a*Cos[q*t], b*Sin[p*t], c*Sin[r*t]};
fyz[a_, b_, c_, q_, p_, r_] := {-3, b*Sin[p*t], c*Sin[r*t]};
fxz[a_, b_, c_, q_, p_, r_] := {a*Cos[q*t], 3, c*Sin[r*t]};
fxy[a_, b_, c_, q_, p_, r_] := {a*Cos[q*t], b*Sin[p*t], -3};
```

Plotting curve and projections.
```mathematica
gPlot[a_, b_, c_, q_, p_, r_, u_] := ParametricPlot3D[{fxyz[a, b, c, q, p, r], fyz[a, b, c, q, p, r], fxz[a, b, c, q, p, r], fxy[a, b, c, q, p, r]}, {t, 0, u}];
```

Animating curve: a=1,b=1,c=1,q=3,p=6,r=8, 0<u<2Pi
```mathematica
Animate[Show[gPlot[1, 1, 1, 3, 6, 8, u], {u, 0, 2*Pi, Pi/45}]
```

![Image of 3D Lissajous](/img/23.3DLissajous.gif)

### 24. Circle Around Curve

Defining 3D curve and projections.
```mathematica
fxyz = {Cos[t] + 2*Cos[2*t], Sin[t] - 2*Sin[2*t], 2*Sin[3*t]};
fyz = ReplacePart[fxyz, 1 -> -range];
fxz = ReplacePart[fxyz, 2 -> range];
fxy = ReplacePart[fxyz, 3 -> -range];
```

Calculating parameter t values.
```mathematica
tvalues = Table[t, {t, Pi/45, 2*Pi, Pi/45}];
```

Defining function that calculates circle equation around any curve.  
1. Depending on t value, calculating point on curve.  
2. Finding derivative of curve.  
3. Finding tangent unit vector to curve.  
4. Normalizing found vector.
5. Finding vector perpendicular to tangent.  
6. Normalizing found vector.  
7. Finding second vector perpendicular to both previous vectors.  
8. Defining 3D circle using two perpendicular vectors (both on the same plane) and one point (center of circle).  
9. Creating projection equations.
```mathematica
createCircleEquations[curve_, radius_, tpoint_, range_] := Module[{circle, circleYZ, circleXZ, circleXY},
   point = N[curve /. t -> tpoint, 5];
   curveDerivative =  D[curve, t];
   curveTangent = N[curveDerivative /. t -> tpoint, 5];
   curveTangentNormalized = curveTangent/Norm@curveTangent;
   normalUnitVector1 = {1, 0, a3 /. First@Solve[curveTangentNormalized[[1]]*a1 + curveTangentNormalized[[2]]*a2 + curveTangentNormalized[[3]]*a3 == 0 /. {a1 -> 1, a2 -> 0},a3]};
   normalUnitVector1Normalized = normalUnitVector1/Norm@normalUnitVector1;
   normalUnitVector2 = Cross[curveTangentNormalized, normalUnitVector1Normalized];
   circle = First@{point + radius*Cos[t]*normalUnitVector1Normalized + radius*Sin[t]*normalUnitVector2};
   circleYZ = ReplacePart[circle, 1 -> -range];
   circleXZ = ReplacePart[circle, 2 -> range];
   circleXY = ReplacePart[circle, 3 -> -range];
   {circle, circleYZ, circleXZ, circleXY}
   ];
```

Using createCircleEquations method and plotting. 
```mathematica
circle[k_] := createCircleEquations[fxyz, 0.5, tvalues[[k]], range];
gPlot[k_] := ParametricPlot3D[{Evaluate@Flatten@circle[k], fxyz, fxy, fxz, fyz}, {t, 0, 2*Pi}];
Animate[Show[gPlot[k], box], {k, 1, Length@tvalues, 1}];
```

![Image of Circle Around Curve](/img/24.CircleAroundCurve.gif)

### 25. Conical Rose [MathCurve](https://www.mathcurve.com/courbes3d.gb/rosaceconique/rosaceconique.shtml)

Defining Conical Rose curve and plotting function.
```mathematica
conicalrose[a_, b_, n_, c_] := {a*Cos[n*t + c]*Cos[t], a*Cos[n*t + c]*Sin[t], b*Cos[n*t + c]};
gPlot[f_, style_] := ParametricPlot3D[f, {t, 0, 2*Pi}, PlotRange -> {{-3, 3}, {-3, 3}, {0, 1}}, style];
```

Animating rose.  
Changing parameter a in rose curve will control rose opening.  
Changing value of ViewPoint (Pi*(1 - j)) from Pi to 0 moves camera.
```mathematica
Animate[Show[
  gPlot[Table[conicalrose[1*j, 0.5, 7, i], {i, -Pi, Pi, Pi/4}], ColorFunction -> Function[{x, y, z}, Hue[Abs[z]/7]]],
  gPlot[Table[conicalrose[2*j, 1, 4, i], {i, -Pi, Pi, Pi/4}], ColorFunction -> Function[{x, y, z}, Hue[Abs[z]/4]]],
  gPlot[Table[conicalrose[3*j, 1, 4, i], {i, -Pi, Pi, Pi/4}], ColorFunction -> Function[{x, y, z}, Hue[Abs[z]/3]]],
  ViewPoint -> {0, Pi*(1 - j), Pi}], {j, 0, 1, 0.01}
];
```

![Image of Conical Rose](/img/25.ConicalRose.gif)

### 26. 3D Basin [MathCurve](https://www.mathcurve.com/courbes3d.gb/vasque3d/vasque3d.shtml)

Defining Basin curve and plotting function.  
Pot function is just basin XY divided by Cos[n*t].
```mathematica
basin[a_, b_, n_, z0_] := {a*Cos[t]*Cos[n*t], a*Sin[t]*Cos[n*t], z0 + b*Cos[n*t]^2};
pot[a_, b_, n_, z0_] := {a*Cos[t], a*Sin[t], z0 + b*Cos[n*t]^2};
gPlot[f_, style_, u_] := ParametricPlot3D[f, {t, 0, u}, PlotStyle -> style, PlotRange -> {{-2, 2}, {-2, 2}, {0, 2.4}}];
```

Defining four parts of object.
```mathematica
bottom[u_] := gPlot[ReplacePart[basin[2, -0.25, 20/7, 1.75], 3 -> 0], Hue[0.1, 0.8, 0.3, 1], u];
wall[u_] := gPlot[pot[2, 1.5, 20/7, 0], Hue[0.4, 0.8, 0.5, 1], u];
lid[u_] := gPlot[basin[2, -0.3, 20/7, 1.8], Hue[0.6, 0.8, 0.6, 1], u];
knob[u_] := gPlot[basin[0.3, -0.3, 10/7, 2.1], Hue[0.8, 0.8, 0.7, 1], u];
```

Animating basin creation.
```mathematica
Animate[Show[bottom[u], wall[u], lid[u], knob[u], {u, Pi/15, 14*Pi, Pi/15}];
```

![Image of 3D Basin](/img/26.3DBasin.gif)

### 27. Langtons ant [Wikipedia](https://en.wikipedia.org/wiki/Langton%27s_ant)

Function for creating nxn matrix.
```mathematica
createMatrix[n_] := Table[0, {n}, {n}];
```

Set of movement rules.  
> At a white square, turn 90° right, flip the color of the square, move forward one unit  
> At a black square, turn 90° left, flip the color of the square, move forward one unit
```mathematica
turnRight[x_, y_, direction_] := Switch[direction,
   "N", {x + 1, y, "E"},
   "E", {x, y - 1, "S"},
   "S", {x - 1, y, "W"},
   "W", {x, y + 1, "N"}];
turnLeft[x_, y_, direction_] := Switch[direction,
   "N", {x - 1, y, "W"},
   "W", {x, y - 1, "S"},
   "S", {x + 1, y, "E"},
   "E", {x, y + 1, "N"}];
```

Choosing move according to cell color.
```mathematica
makeMove[{x_, y_, direction_}, matrix_] := Switch[matrix[[#1]][[#2]], 0, turnRight[#1, #2, #3], 1, turnLeft[#1, #2, #3]] &[x, y, direction];
```

Toggling cell color.
```mathematica
toggleValue[{x_, y_, direction_}, matrix_] := ReplacePart[matrix, {#1, #2} -> If[matrix[[#1]][[#2]] == 0, 1, 0]] &[x, y];
```

Defining function to recurse.
```mathematica
LangtonsAnt[{ant_, matrix_}] := {makeMove[ant, matrix], toggleValue[ant, matrix]};
```

Visualizing matrix.
```mathematica
visualizeMatrix[{{x0_, y0_, direction_}, matrix_}] :=
  Graphics[Table[{EdgeForm[Thin], If[{x, y} == {x0, y0}, Red, If[matrix[[x]][[y]] == 1, Black, White]], 
     Rectangle[{x, y}, {x + 1, y + 1}]}, {x, 1, Length@matrix}, {y, 1,Length@matrix}]];
```

Using code. Setting initial X,Y position, ant direction, creating matrix and defining steps count.  
Animated image shows 100x100 matrix and 11500 steps (skipping over 10 to make faster).
```mathematica
init = {50, 50, "N"};
matrix = createMatrix[100];
steps = 11500;
result = Nest[LangtonsAnt, {init, matrix}, steps];
visualizeMatrix[result]
```

![Image of Langtons ant](/img/27.LangtonsAnt.gif)

### 28. 2D Rotation matrix [Wikipedia](https://en.wikipedia.org/wiki/Rotation_matrix)

Building clock using rotation matrix. Matrix is modified so movement is clock-wise.  
Using dot product of rotation matrix and direction vector to build line. Line can be used as hand of clock.  
n - rotation speed (used for hours, minutes, seconds).  
v - line length.
```mathematica
line[t_, n_] := Dot[{{Cos[t/n], Sin[t/n]}, {-Sin[t/n], Cos[t/n]}}, {0, v}];
gLine[f_, u_, style_] := ParametricPlot[f, {v, 0, u}, PlotStyle -> style];
```

Calculating angles and coordinates for clock numbers using parametric circle equation.
```mathematica
coordinates = Table[{1*Cos[i], 1*Sin[i]}, {i, Pi/3, -3*Pi/2, -Pi/6}];
gNumbers = Table[Graphics@Text[RomanNumeral@i, coordinates[[i]]], {i, 1, 12}];
```

Animating clock hands.
```mathematica
Animate[Show[gLine[line[t, 1], 0.75, Black], gLine[line[t, 12], 0.5, Red], gNumbers, ImageSize -> Large], {t, 0, 24*Pi}]
```

![Image of 2D Rotation matrix](/img/28.2DRotationMatrix.gif)

### 29. 3D Rotation matrix [Wikipedia](https://en.wikipedia.org/wiki/Rotation_matrix)

Defining function that rotates lines.  
1. Defining rotation matrix around z-axis.  
2. Defining perpendicular line to xy-plane, that is going to rotate around z-axis.  
3. Defining circle equation, our line moves on circlular trajectory, using it to generate points.  
4. Defining point.  
5. Creating tangent vector.  
6. Creating normal vector.  
7. Replacing v=0, we just need direction.  
8. Using RotationTransform to build second rotation matrix.  
9. Using rotation matrix to rotate line around its own axis too.
```mathematica
rotatingLine[t_] := Module[{},
   zRotationMatrix = {{Cos[t], -Sin[t], 0}, {Sin[t], Cos[t], 0}, {0, 0, 1}};
   perpendicularLine = Dot[zRotationMatrix, {1, 0, v}];
   circle[u_] := {Cos[u], Sin[u], 0};
   point = N[circle[t], 3];
   tangent = point + N[circle'[t], 3]*v;
   normal = Cross[tangent, {0, 0, 1}];
   axis = normal /. v -> 0;
   pRotationMatrix = RotationTransform[t, axis, point];
   line = Simplify[pRotationMatrix[perpendicularLine]];
   line
   ];
```

Creating lines and 3 projections.
```mathematica
lines = {
     rotatingLine[#],
     ReplacePart[rotatingLine[#], 1 -> -range],
     ReplacePart[rotatingLine[#], 2 -> range],
     ReplacePart[rotatingLine[#], 3 -> -range]
     } & /@ Range[Pi/72, 2*Pi, Pi/72];
```

Plotting and animating.
```mathematica
gLines[i_, style_] := ParametricPlot3D[lines[[i]], {v, -0.5, 0.5}];
gMesh = Show[gLines[#, GrayLevel[.5, .7]] & /@ Range[1, Length@lines]];
Animate[Show[box, gMesh, gLines[i, Red]], {i, 1, Length@lines, 1}];
```
Created surface is known as [Möbius strip](https://en.wikipedia.org/wiki/M%C3%B6bius_strip)  

![Image of 3D Rotation matrix](/img/29.3DRotationMatrix.gif)

### 30. Sierpinski triangle [Wikipedia](https://en.wikipedia.org/wiki/Sierpinski_triangle) [Oftenpaper](http://www.oftenpaper.net/sierpinski.htm)

Defining function for splitting triangles. Input is list of polygons (triangles).
```mathematica
SierpinskiTriangle[polygons_] := Module[{},
   newpolygons = {};
    Table[{
     {p1, p2, p3} = polygons[[i]];
     polygon1 = {p1, (p1 + p2)/2, (p1 + p3)/2};
     polygon2 = {p2, (p2 + p3)/2, (p1 + p2)/2};
     polygon3 = {p3, (p1 + p3)/2, (p2 + p3)/2};
     newpolygons = Join[newpolygons, {polygon1, polygon2, polygon3}]
     }, {i, 1, Length@polygons}];
   newpolygons
   ];
```

Define list of initial polygons (just 1 triangle), go through it recursively n = 8 and plot result.
```mathematica
polygons = {{{0, 0}, {1, Sqrt[3]}/2, {1, 0}}};
result = NestList[SierpinskiTriangle, polygons, 8];
Graphics[{White, EdgeForm[Directive[Thin, Blue]], Polygon[Last@result]}]
```

![Image of Sierpinski triangle](/img/30.SierpinskiTriangle.gif)