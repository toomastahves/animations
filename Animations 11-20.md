### 11. "A Bird in Flight (2014)," by Hamid Naderi Yeganeh [AMS](http://ams.org/mathimagery/displayimage.php?album=40&pid=570#top_display_media)

Creating two sets of points.
```mathematica
n = 2000;
points1 = Table[{3*Sin[2*Pi*i/n]^3, -Cos[8*Pi*i/n]}, {i, 1, n}];
points2 = Table[{3/2*Sin[2*Pi*i/n]^3, -1/2*Cos[6*Pi*i/n]}, {i, 1, n}];
```

Creating graphical points and lines.
```mathematica
gPoints1 = Table[Graphics[{PointSize[0.001], Point[points1[[i]]]}], {i, 1, n}];
gPoints2 = Table[Graphics[{PointSize[0.001], Point[points2[[i]]]}], {i, 1, n}];
gLines = Table[Graphics[{GrayLevel[.2, .2], Thin, Line[{points1[[i]], points2[[i]]}]}], {i, 1, n}];
```

Unfolding lists for animation.
```mathematica
gPoint1Frames = FoldList[List, gPoints1];
gPoint2Frames = FoldList[List, gPoints2];
gLineFrames = FoldList[List, gLines];
```

![Image of A Bird in Flight 2014](/img/11.abirdinflight2014.gif)


### 12. "A Bird in Flight (2015)," by Hamid Naderi Yeganeh [AMS](http://ams.org/mathimagery/displayimage.php?album=40&pid=616#top_display_media)

Creating two sets of points.
```mathematica
n = 500;
points1 = Table[{3/2*Sin[2*Pi*i/n + Pi/3]^7, 1/4*Cos[6*Pi*i/n]^2}, {i, 1, n}];
points2 = Table[{1/5*Sin[6*Pi*i/n + Pi/5], -2/3*Sin[2*Pi*i/n - Pi/3]^2}, {i, 1, n}];
```

Creating graphical points and lines.
```mathematica
gPoints1 = Table[Graphics[{PointSize[0.001], Point[points1[[i]]]}], {i, 1, n}];
gPoints2 = Table[Graphics[{PointSize[0.001], Point[points2[[i]]]}], {i, 1, n}];
gLines = Table[Graphics[{GrayLevel[.1, .6], Thin, Line[{points1[[i]], points2[[i]]}]}], {i, 1, n}];
```

Unfolding lists for animation.
```mathematica
gPoint1Frames = FoldList[List, gPoints1];
gPoint2Frames = FoldList[List, gPoints2];
gLineFrames = FoldList[List, gLines];
```

![Image of A Bird in Flight 2015](/img/12.abirdinflight2015.gif)

### 13. Starr Rose [Wolfram MathWorld](http://mathworld.wolfram.com/StarrRose.html)

Defining Starr Rose in parametric form.  
a,b,c are curve parameters. n controls scale of Rose.
```mathematica
starr[a_, b_, c_, n_] := {n*(2 + Sin[a *t]/2)*Cos[(t + Sin[b *t]/c)], n*(2 + Sin[a *t]/2)*Sin[(t + Sin[b *t]/c)]};
gPlot[f_, u_] := ParametricPlot[f, {t, 0, u}, PerformanceGoal -> "Speed", PlotStyle -> Hue[0.85, 0.9, 0.5, 0.6], Axes -> False, PlotRange -> {{-range, range}, {-range, range}}];
```

Creating roses. a=8 is constant. t values (0; 2*Pi).  
Inner loop {i, 0.05, 2.9, 0.05} creates multiple Roses on top of each other.  
Outer loop {j, 0.05, 16, 0.05} changes parameters b and c.  
Last Rose has parameters a=8, b=16, c=16.
```mathematica
frames = Table[Show[Table[gPlot[N[starr[8, j, j, i], 4], 2*Pi], {i, 0.05, 2.9, 0.05}]], {j, 0.05, 16, 0.05}]; // Timing
```

![Image of Starr Rose](/img/13.starrrose.gif)

### 14. Maurer Rose [Wikipedia](https://en.wikipedia.org/wiki/Maurer_rose)

Maurer rose function just creates line between two points on a given curve.  
Function inputs: n - Rose curve parameter, d - step between two points, k - amount of points on a curve.
```mathematica
MauerRose[n_, d_, k_] := Table[Graphics[{Hue[0.5, 0.8, 0.5, 0.8], Line[
      {{Sin[n *2*Pi*t/k]*Cos[2*Pi*t/k],Sin[n*2*Pi*t/k]*Sin[2*Pi*t/k]},
       {Sin[n *2*Pi*(t + d)/k]*Cos[2*Pi*(t + d)/k],Sin[n*2*Pi*(t + d)/k]*Sin[2*Pi*(t + d)/k]}}
      ]}], {t, 1, k}];
```

Animation below has n = 10, k = 360 and d changes dynamically from 1 to 360.
```mathematica
Show[MauerRose[10, d, 360]
```

![Image of Maurer Rose](/img/14.maurerrose.gif)

### 15. Truchet Tiling [Wolfram MathWorld](http://mathworld.wolfram.com/TruchetTiling.html)

Define function to generate tiles.
```mathematica
generateGrid[n_] := Table[
   Rotate[{Circle[{x, y}, 1, {0, Pi/2}], Circle[{x + 2, y + 2}, 1, {Pi, 3*Pi/2}]}, RandomChoice[{0, Pi/2}]],
   {x, 0, 2*n - 1, 2},
   {y, 0, 2*n - 1, 2}];
grid = generateGrid[10];
```

Define function to rotate each tile for animation.
```mathematica
rotateGrid[grid_, n_, u_] := 
  Graphics[Table[Rotate[grid[[i]][[j]], u], {i, 1, n}, {j, 1, n}], AspectRatio -> Automatic, PlotRange -> {{0, 2*n - 1}, {0, 2*n - 1}}];
Animate[rotateGrid[grid, 10, u], {u, 0, Pi}];
```

![Image of Truchet Tiling](/img/15.truchettiling.gif)

### 16. Pacman Curve [Mathematics Stack Exchange](https://math.stackexchange.com/questions/418641/smooth-pac-man-curve)

Define Pacman curve in polar coordinates. Parameter 'i' controls mouth openness. Parameter 'side' controls which side it opens.
```mathematica
pacman[i_, side_] := PolarPlot[1 + Tanh[10000*(i + side*Cos[t])], {t, 0, 2*Pi}, Axes -> False, PlotStyle -> Hue[0.16, 1, 1, 1]];
```

Define friend function in parametric form. Parameter 'style' controls what clothes it wears.
```mathematica
friend[style_] := ParametricPlot[{{t - Sin[t], 3.5 - 3*Cos[t]}, {t, 0.5*Cos[4*t]}}, {t, 0, 2*Pi}, Axes -> False, AspectRatio -> Automatic, PlotStyle -> style];
```

Scale function is necessary because we are using one animation frame, but two different movements with different speed.  
Pacman moves left/right (slower) and opens mouth (faster).
```mathematica
scale[i_, min_, max_] := (0.95 - 0.3)*(i - min)/(max - min) + 0.3;
```

Defining empty plot to draw objects on. Using Inset to place object onto empty plot. 
```mathematica
plot[pacman_, friend_, x0_, y0_, side_] := ParametricPlot[, {t, 0, 2*Pi}, Axes -> False, PlotRange -> {{-range, range}, {-range, range}}, 
   Epilog -> {Inset[pacman, {x0, y0}, {0, 0}, {1, 1}], Inset[friend, {x0 - side*1.5, y0 - 0.2}, {0, 0}, {0.5, 0.5}]}];
```

Animating movement. Floor function rounds 'i' to closest Integer, so scaling numbers between (-5,-4), (-4,-3) etc.
```mathematica
Animate[Show[plot[pacman[scale[i, Floor[i], Floor[i + 1]], -1], friend[Blue], i, -1, -1]], {i, -5, 4, 0.1}];
```

![Image of Pacman Curve](/img/16.pacman.gif)

### 17. Superellipse [Wolfram MathWorld](http://mathworld.wolfram.com/Superellipse.html) [Wikipedia](https://en.wikipedia.org/wiki/Superellipse)

17.1. Defining superellipse in parametric form.
```mathematica
gPlot1[a_, b_, c_] := ParametricPlot[{Abs[Cos[t]]^(2/c)*a*Sign[Cos[t]], Abs[Sin[t]]^(2/c)*b*Sign[Sin[t]]}, {t, Pi/180, 2*Pi + Pi/180}];
```

Plotting multiple superellipses by changing parameter c in range [-5;5] \ {0}
```mathematica
frames1 = Table[gPlot1[2, 2, i], {i, Join[Range[-5, -0.1, 0.1], Range[0.1, 5, 0.1]]}];
```

![Image of Superellipse 1](/img/17.1.Superellipse.gif)

17.2. Defining superellipse in polar coordinates.
```mathematica
gPlot2[n1_, n2_, n3_, a_, b_, m_, u_] := PolarPlot[(Abs[Cos[m*t/4]/a]^n2 + Abs[Sin[m*t/4]/b]^n3)^-(1/n1), {t,Pi/18, u}];
```

Plotting multiple superellipses by changing parameters n1=n2=n3 in range [-5;5] \ {0}
```mathematica
frames2 = Table[gPlot2[i, i, i, 1, 1, 8, 2*Pi + Pi/18], {i, Join[Range[-5, -0.1, 0.1], Range[0.1, 5, 0.1]]}];
```

![Image of Superellipse 2](/img/17.2.Superellipse.gif)

### 18. Meander curve [Mathcurve](https://www.mathcurve.com/courbes2d.gb/giration/motifs.shtml)

 Defining curve and parametric plotting function.  
 a - size of curve  
 rotate - rotates curve around P(0,0), changes [0;2*Pi]  
 v - length of curve  
 i - curve parameter
```mathematica
gPlot[a_, i_, rotate_, v_] := ParametricPlot[
   {a*NIntegrate[Cos[i*Sin[u] + rotate], {u, 0, t}], a*NIntegrate[Sin[i*Sin[u] + rotate], {u, 0, t}]}, 
   {t, -v, v}, Axes -> False, PlotStyle -> Hue[0.8, 1, 0.5, 0.6], PlotRange -> {{-range, range}, {-range, range}}];
```

Plotting four curves, each rotated Pi/2. Parameter i changes [-4*Pi;4*Pi]
```mathematica
Show[{gPlot[2, i, 0, 2*Pi], gPlot[2, i, Pi/2, 2*Pi], gPlot[2, i, Pi, 2*Pi], gPlot[2, i, 3*Pi/2, 2*Pi]}]
```

![Image of Meander curve](/img/18.Meander.gif)

### 19. Pendulum [Wikipedia](https://en.wikipedia.org/wiki/Pendulum_(mathematics))

Pendulum parameters. l - length, m - mass, t2 - max oscillation time, b - damping ratio.
```mathematica
g = 9.81; l = 1; b = 0.1; m = 1; t2 = 30;
```

Defining and solving differential equations with damping. Table solves 18 equations. a is starting position of pendulum.
```mathematica
interpols = Table[y /. NDSolve[{y''[t] + g/l*y[t] + b/m*y'[t] == 0, y'[0] == 0, y[0] == a}, y, {t, 0, t2}][[1]], {a, Pi/18, Pi, Pi/18}];
```

Describing pendulums as line with a circle at its end. Rotating line around point P(0,0) and angle comes from dif. eq. solution.
```mathematica
pendulum[t_, i_] := Graphics@Rotate[{Hue[i/10], Line[{{0, 0}, {0, -l}}], Circle[{0, -l}, 0.05]}, interpols[[i]][t], {0, 0}];
```

Animating pendulums.
```mathematica
Animate[Show[box, Table[pendulum[t, i], {i, 1, 18}]], {t, 0, t2, 0.05}]
```

![Image of Pendulum](/img/19.Pendulum.gif)

### 20. Double Pendulum [Math24](https://www.math24.net/double-pendulum/)

Pendulum parameters. l1,l2 - pendulum length, m1,m2 - mass, t2 - max oscillation time, a1,a2 - starting angle
```mathematica
g = 9.81; m1 = 2; m2 = 1; l1 = 1; l2 = 1; t2 = 30; a1 = Pi; a2 = Pi/2;
```

Defining and solving Lagrange differential equations.
```mathematica
lagrangian = {(m1 + m2)*l1*y1''[t] + m2*l2*y2''[t]*Cos[y1[t] - y2[t]] + m2*l2*y2'[t]^2*Sin[y1[t] - y2[t]] + (m1 + m2)*g*Sin[y1[t]] == 0,
   l2*y2''[t] + l1*y1''[t]*Cos[y1[t] - y2[t]] - l1*y1'[t]^2*Sin[y1[t] - y2[t]] + g*Sin[y2[t]] == 0,
   y1[0] == a1, y2[0] == a2, y1'[0] == 0,  y2'[0] == 0};
```

Getting results from solution. theta1 is angle of first pendulum. theta2 is angle of second pendulum.
```mathematica
intpol = NDSolve[lagrangian, {y1, y2}, {t, 0, t2}][[1]];
theta1 = y1 /. intpol;
theta2 = y2 /. intpol;
```

Defining points depending on time.
```mathematica
p0 = {0, 0};
p1[t_] := {l1*Sin[theta1[t]], -l1*Cos[theta1[t]]};
p2[t_] := {l1*Sin[theta1[t]] + l2*Sin[theta2[t]], -l1*Cos[theta1[t]] - l2*Cos[theta2[t]]};
```

Visualizing pendulums and trajectory.
```mathematica
pendulum1[t_] := Graphics@{Red, Line[{p0, p1[t]}], Circle[p1[t], 0.05]};
pendulum2[t_] := Graphics@{Blue, {Line[{p1[t], p2[t]}], Circle[p2[t], 0.05]}};
trajectory[u_] = ParametricPlot[Evaluate[p2[t]], {t, 0, u}, PlotRange -> {{-(l1 + l2), (l1 + l2)}, {-(l1 + l2), (l1 + l2)}}, Axes -> False, PlotStyle -> {GrayLevel[0.1, 0.2]}];
```

![Image of Double Pendulum](/img/20.DoublePendulum.gif)