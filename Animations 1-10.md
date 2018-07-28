### 1. Archimedean spirals [Wikipedia](https://en.wikipedia.org/wiki/Archimedean_spiral)

Defining function of the spiral. Spiral parameters: a (inclination in xy-plane), b (distance between intersception points of the spiral on a straight line) and theta (length of the spiral).
```mathematica
spiral[theta_, a_, b_] := a + b*theta;
```

Drawing circles using polar plot, where i*b*Pi is circle formula in polar coordinates. Multiplying by i increases radius.
```mathematica
gCircles[b_] := Table[PolarPlot[i*b*Pi, {theta, 0, 2*Pi}, Ticks -> False, PlotStyle -> {Dashed, Black, Thin}], {i, 0, count}];
```

Drawing spiral using polar plot. Angle theta changes from 0 to arg, Changing theta is necessary when animating.
```mathematica
gSpiral[arg_, a_, b_, style_] := PolarPlot[spiral[theta, a, b], {theta, 0, arg}, PlotStyle -> {style}];
```

Animating two spirals. One spiral has positive {a,b} parameters, second one negative {-a,-b} parameters.
```mathematica
Animate[Show[gCircles[b], gSpiral[arg, a, b, Blue], gSpiral[arg, -a, -b, Red]], {arg, 0, count*Pi}, AnimationRate -> 0.025];
```

![Image of Spirals](/img/1.spirals.gif)

### 2. Lissajous curve [Wikipedia](https://en.wikipedia.org/wiki/Lissajous_curve)

Defining function of the curve in parametric form.
```mathematica
lissajous[a_, b_, phasex_, phasey_, kx_, ky_] := {a*Sin[kx*t + phasex], b*Sin[ky*t + phasey]};
```

Plotting curve.
```mathematica
gPlot[theta_, a_, b_, phasex_, phasey_, kx_, ky_] := ParametricPlot[lissajous[a, b, phasex, phasey, kx, ky], {t, 0, theta},  PlotRange -> {{-5, 5}, {-5, 5}},  Axes -> False];
```

When code is run inside Wolfram Mathematica, then parameters can be changed during running animation.
```mathematica
Manipulate[Animate[Show[gPlot[theta, a, b, phasex, phasey, kx, ky], gPoint[a*Sin[kx*theta + phasex], b*Sin[ky*theta + phasey]]], 
{theta, 0, 2*Pi}], {a, -5, 5}, {b, -5, 5}, {phasex, 0, 2*Pi}, {phasey, 0, 2*Pi}, {kx, -50, 50}, {ky, -50, 50}];
```

![Image of Lissajous curve](/img/2.lissajous.gif)

### 3. Hypocycloid [Wikipedia](https://en.wikipedia.org/wiki/Hypocycloid)

Defining hypocycloid formula.
```mathematica
hypocycloid[ t_, R_, r_] := {(R - r)*Cos[t] + r*Cos[(R - r)/r*t], (R - r)*Sin[t] - r*Sin[(R - r)/r*t]};
gHypocycloid[u_, R_, r_, style_] := ParametricPlot[hypocycloid[t, R, r], {t, 0, u}, Axes -> False, PlotRange -> {{-22, 22}, {-22, 22}}, PlotStyle -> {style}];
```

Generating multiple hypocycloids with different parameters and colors.
```mathematica
gShow2[u_, R_, r_] := Show[
   gHypocycloid[u, 2, r, Red],
   gHypocycloid[u, 3, r, Green],
   gHypocycloid[u, 4, r, Blue],
   gHypocycloid[u, 5, r, Cyan],
   gHypocycloid[u, 6, r, Magenta],
   gHypocycloid[u, 7, r, Yellow],
   gHypocycloid[u, 8, r, Brown],
   gHypocycloid[u, 9, r, Orange],
   gHypocycloid[u, 10, r, Pink],
   gHypocycloid[u, 11, r, Red],
   gHypocycloid[u, 12, r, Green],
   gHypocycloid[u, 13, r, Blue],
   gHypocycloid[u, 14, r, Cyan],
   gHypocycloid[u, 15, r, Magenta],
   gHypocycloid[u, 16, r, Yellow],
   gHypocycloid[u, 17, r, Brown],
   gHypocycloid[u, 18, r, Orange],
   gHypocycloid[u, 19, r, Pink],
   gHypocycloid[u, 20, r, Purple],
   gHypocycloid[u, 21, r, Gray],
   gCircle[0, 0, R], ImageSize -> Large];
Animate[gShow2[u, 21, 1], {u, Pi/180, 2*Pi, Pi/180}];
```

![Image of Hypocycloids](/img/3.hypocycloids.gif)

### 4. Hypotrochoid [Wikipedia](https://en.wikipedia.org/wiki/Hypotrochoid)

Defining hypotrochoid.
```mathematica
hypotrochoid[ t_, R_, r_, d_] := {(R - r)*Cos[t] + d*Cos[(R - r)/r*t], (R - r)*Sin[t] - d*Sin[(R - r)/r*t]};
gHypotrochoid[u_, R_, r_, d_, style_] := ParametricPlot[hypotrochoid[t, R, r, d], {t, 0, u}, Axes -> False, PlotRange -> {{-70, 70}, {-70, 70}}, PlotStyle -> {style}];
```

Generating hypotrochoid with parameters R=73, r=37, d=17. Interesting patterns emerge when using prime numbers.
```mathematica
gShow[u_, R_, r_, d_] := Show[gHypotrochoid[u, R, r, d, Pink], ImageSize -> Large];
Animate[gShow[u, 73, 37, 17], {u, Pi/180, 74*Pi}];
```

![Image of Hypotrochoid](/img/4.hypotrochoid.gif)


### 5. Epicycloid [Wikipedia](https://en.wikipedia.org/wiki/Epicycloid)

Defining epicycloid.
```mathematica
epicycloid[ t_, R_, r_] := {(R + r)*Cos[t] - r*Cos[(R + r)/r*t], (R + r)*Sin[t] - r*Sin[(R + r)/r*t]};
gEpicycloid[u_, R_, r_, style_] := ParametricPlot[epicycloid[t, R, r], {t, 0, u},  Axes -> False, PlotRange -> {{-90, 90}, {-90, 90}}, PlotStyle -> {style}];
```

Generating 4 epicycloids on top of each other. Parameters {1,21,41,61} are radiuses of big circle. Parameter {10} is radius of small circle (not visible).
```mathematica
gShow[u_] := Show[
   gEpicycloid[u, 1, 10, Purple],
   gEpicycloid[u, 21, 10, Pink],
   gEpicycloid[u, 41, 10, Orange],
   gEpicycloid[u, 61, 10, Magenta],
   ImageSize -> Large
   ];  
Animate[gShow[u ], {u, Pi/72, 20*Pi, Pi/72}];
```

![Image of Epicycloid](/img/5.epicycloids.gif)

### 6. Rotating Circle

Defining and plotting curve. Sophisticated function with multiple parameters. 
```mathematica
gCircle[u_, xa_, ya_, xb_, yb_, xc_, yc_, xd_, yd_, e_] := 
  ParametricPlot[{xa*(xb*Cos[xc*t] + Cos[xd*t]*Cos[e*t]), ya*(yb*Sin[yc*t] + Sin[yd*t] )}, {t, 0, u},
   PlotStyle -> {Thin},  PlotPoints -> 1000,   Axes -> False, PlotRange -> {{-range, range}, {-range, range}}];
```

Playing with circle, rotating around own axis and putting on the move on one of the Lissajous curves.
```mathematica
gShow[u_] := Show[
   gCircle[u,
    1, 1, (* xa, ya - radius of circle *)
    8, 8, (* xb, yb - range on plot *)
    3, 4, (* xc, yc - Lissajous curve parameters *)
    720, 720, (* xd, yd - amount of rotations *)
    1],(* e - rotating circle around own axis *)
   ImageSize -> Large];
```


![Image of Rotating circle](/img/6.circle.gif)

### 7. Heart Curve [Wolfram MathWorld](http://mathworld.wolfram.com/HeartCurve.html)

Defining and plotting heart curve. x0 and y0 are heart position parameters on xy-plane. 
```mathematica
gHeart[x0_, y0_] := 
  ParametricPlot[{x0 + 16*Sin[t]^3, y0 + 13*Cos[t] - 5*Cos[2*t] - 2*Cos[3*t] - Cos[4*t]},
  {t, 0, 2*Pi}, PlotStyle -> {Thin, Pink},  Axes -> False, PlotRange -> {{-range, range}, {-range, range}}];
```

Creating hearts in different places in xy-plane on heart trajectory.
```mathematica
hearts = Table[gHeart[16*Sin[u]^3, 13*Cos[u] - 5*Cos[2*u] - 2*Cos[3*u] - Cos[4*u]], {u, 0, 2*Pi, Pi/90}];
```

Heart itself is continuous function. 
But moving heart on heart trajectory is not continuous. To make it look continuous we unfold and animate frames.
```mathematica
frames = FoldList[List, hearts];
```

![Image of Heart Curve](/img/7.heart.gif)

### 8. Evolute [Wikipedia](https://en.wikipedia.org/wiki/Evolute)

Defining function. In this case Rose curve with parameter k=2. 
```mathematica
f = {Cos[t] + Cos[3*t], Sin[t] - Sin[3*t]};
```

Defining evolute formula
```mathematica
evolute[f_] := {
  f[[1]] - (D[f[[2]],t]*(D[f[[1]],t]^2 + D[f[[2]],t]^2))/(D[f[[1]], t]*D[D[f[[2]], t], t] - D[f[[2]], t]*D[D[f[[1]], t], t]),
  f[[2]] + (D[f[[1]],t]*(D[f[[1]],t]^2 + D[f[[2]],t]^2))/(D[f[[1]], t]*D[D[f[[2]], t], t] - D[f[[2]], t]*D[D[f[[1]], t], t])
};
```

Finding dy/dx
```mathematica
dydx[f_] := D[f[[2]], t]/D[f[[1]], t];
```

Creating points on Rose curve.
```mathematica
points = Table[f /. t -> u, {u, 0, 2*Pi, Pi/90}];
```

Calculating slopes of Rose curve.
```mathematica
slopes = Quiet[Table[dydx[f] /. Solve[{f[[1]] == points[[i]][[1]], f[[2]] == points[[i]][[2]]},t][[1]] /. C[1] -> 0, 
{i, 1, Length[points]}]];
```

Building normals of Rose curve.  
When animating normals, then their movement generates hypocycloid with parameter k = 8.  
Normals of Rose (k=2) are tangents of hypocycloid (k=8).
```mathematica
normals = Quiet[Table[-1/slopes[[i]]*(x - points[[i]][[1]]) + points[[i]][[2]], {i, 1, Length[slopes]}]];
```

![Image of Evolute](/img/8.evolute.gif)


### 9. Epitrochoid [Wikipedia](https://en.wikipedia.org/wiki/Epitrochoid)

Defining and plotting epitrochoid. 
R - radius of big circle, r - radius of small circle, d - distance from center of small circle
```mathematica
f[R_, r_, d_] := {(R + r)*Cos[t] - d*Cos[(R + r)/r*t], (R + r)*Sin[t] - d*Sin[(R + r)/r*t]};
gPlot[f_, u_, style_] := ParametricPlot[f, {t, 0, u}, PlotRange -> {{-range, range}, {-range, range}}, 
   PlotStyle -> Hue[style], AspectRatio -> Automatic, Axes -> False];
```

Creating 30 epitrochoids with different parameters.  
Big circle radiuses vary from 1 to 30.  
Small circle radiuses vary from 1/7 to 30/7.  
Distance from center of small circle is constant 3.
```mathematica
Animate[Show[Table[gPlot[f[i, i/7, 3], u, 0.75], {i, 1, 30}]], {u, 0, 2*Pi, Pi/45}]
```

![Image of Epitrochoid](/img/9.epitrochoid.gif)

### 10. Cornu Spiral [Wolfram MathWorld](http://mathworld.wolfram.com/CornuSpiral.html) [Wikipedia](https://en.wikipedia.org/wiki/Euler_spiral)

Defining Cornu spiral (clothoid). Tuning NIntegrate parameters to calculate faster.
```mathematica
clothoid[fi_, t_?NumericQ, x0_, y0_, a_] := {
   a*NIntegrate[Sin[fi], {s, 0, t}, Method -> {GlobalAdaptive, SymbolicProcessing -> 0, SingularityHandler -> None, Method -> {GaussKronrodRule, Points -> 9}}] + x0,
   a*NIntegrate[Cos[fi], {s, 0, t}, Method -> {GlobalAdaptive, SymbolicProcessing -> 0, SingularityHandler -> None, Method -> {GaussKronrodRule, Points -> 9}}] + y0
   };
```

Creating 30 clothoids with different parameters. Base angle is s^3/3 - 2.19*s. Curvature (derivative of angle) s^2 - 2.19.
```mathematica
gPlot[clothoid[s^3/3 - 2.19*s + Pi/2, t, 0, 0, 1.13], 2*Pi],
```

![Image of Clothoid](/img/10.clothoid.gif)
