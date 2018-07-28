### 31. Gastropod shell [Wikipedia](https://en.wikipedia.org/wiki/Gastropod_shell) [MathCurve](https://www.mathcurve.com/surfaces.gb/coquillage/coquillage.shtml)

Defining equation.
```mathematica
gastropod[u_, v_, m_, b_] := {E^(m*u)*Cos[u]*(1 + b*Cos[v]), E^(m*u)*Sin[u]*(1 + b*Cos[v]), E^(m*u)*(1 + b*Sin[v])};
```

Applying texture while plotting.
```mathematica
gPlot[t_] := 
  ParametricPlot3D[
   gastropod[u, v, .1, .5], {u, -30, t}, {v, 0, 2*Pi},  
   PlotPoints -> 100, ViewPoint -> {0, -Pi/2, -Pi/2}, 
   ViewVertical -> {0, Pi, 0},
   PlotRange -> {{-range, range}, {-range, range}, {0, 2.3}}, 
   Axes -> False, Boxed -> False, Mesh -> None,
   PlotStyle -> {Directive[
      Hue[.2, 1, 1, .8],
      Specularity[White, 100],
      Texture[ExampleData[{"Texture", "Sand3"}]]
      ]},
   TextureCoordinateFunction -> ({10*#4} &)];
```

![Image of Gastropod shell](/img/31.GastropodShell.gif)

### 32. Koch snowflake [Wikipedia](https://en.wikipedia.org/wiki/Koch_snowflake)

Defining drawing function.  
1. Defining generator line path. Not necessary for default Koch curve.  
2. Defining 3 curves and 2D rotation matrices to build equilateral triangle and symmetric snowflake.
```mathematica
draw[n_] := Module[{},
   path = {{1, 0}, {1, Pi/3}, {1, -2*Pi/3}, {1, Pi/3}};
   Graphics[{
     GeometricTransformation[KochCurve[n, path], RotationTransform[Pi, {0.5, 0}]],
     GeometricTransformation[KochCurve[n, path], RotationTransform[-Pi/3, {1, 0}]],
     GeometricTransformation[KochCurve[n, path], RotationTransform[Pi/3, {0, 0}]]}]
   ];
```

Drawing Koch snowflake.
```mathematica
Show[draw[5]]
```

![Image of Koch snowflake](/img/32.KochSnowflake.gif)

### 33. LSystem [Wikipedia](https://en.wikipedia.org/wiki/L-system)

Defining function to build axiom according to rules.
```mathematica
buildAxiom[{axiom_, rules_}] := {StringReplace[axiom, rules], rules};
```

Defining basic L-system function.  
1. Building axiom string, recursion depth n.  
2. Alphabet is (F,f,+,-).  
> "F" : Move forward a step of length d. The state of the turtle changes to (x, y, α), where x = x0 + d cos α and y = y0 + d sin α. A line segment between points (x0, y0) and (x, y) is drawn.  
> "f" : Move forward a step of length d without drawing a line.  
> "+" : Turn left by angle δ. The next state of the turtle is (x, y, α+δ). The positive orientation of angles is counterclockwise.  
> "−" : Turn right by angle δ. The next state of the turtle is (x, y, α − δ).
```mathematica
LSystem[axiom_, rules_, angle_, n_] := Module[{},
   chars = First@Characters@Nest[buildAxiom, {axiom, rules}, n];
   x0 = y0 = a = 0;
   lines = {};
   previouspoint = {0, 0};
   For[i = 1, i <= Length@chars, i++, {Switch[chars[[i]],
       "F", {
        newpoint = {x0 = x0 + Cos[a], y0 = y0 + Sin[a]}; 
        lines = Append[lines, {previouspoint, newpoint}];
        previouspoint = newpoint},
       "f", {
        previouspoint = {x0 = x0 + Cos[a], y0 = y0 + Sin[a]}
        },
       "+", a = a + angle,
       "-", a = a - angle];
     }];
   lines
   ];
```

Trivia: Previously drawn Koch snowflake can be drawn using these parameters:
```mathematica
axiom="F--F--F"; 
rule="F+F--F+F"; 
angle=Pi/3;
```

Drawing island and lakes using LSystem.
```mathematica
axiom = "F+F+F+F";
rules = {"F" -> "F+f-FF+F+FF+Ff+FF-f+FF-F-FF-Ff-FFF", "f" -> "ffffff"};
angle = Pi/2;
lines = LSystem[axiom, rules, angle, 2];
Graphics[Line[lines]]
```

![Image of LSystem](/img/33.LSystem.gif)

### 34. 2D Hilbert curve [Wikipedia](https://en.wikipedia.org/wiki/Hilbert_curve)

Using previously written LSystem function to calculate hilbert curve.
```mathematica
axiom = "-L";
rules = {"L" -> "LF+RFR+FL-F-LFLFL-FRFR+", "R" -> "-LFLF+RFRFR+F+RF-LFL-FR"};
angle = Pi/2;
lines = LSystem[axiom, rules, angle, 3];
gLines = Table[Line[lines[[i]]], {i, 1, Length@lines}];
Graphics[gLines];
```

![Image of 2D Hilbert Curve](/img/34.2DHilbertCurve.gif)

### 35. 3D Hilbert curve [Wikipedia](https://en.wikipedia.org/wiki/Hilbert_curve)

Defining LSystem3D function.  
1. Building axiom according to rules.  
2. Defining initial point {0,0,0} and initial direction vector as dv = IdentityMatrix[3].  
3. Angle is Pi/2, and rotations are around X,Y,Z axes, clockwise and counter-clockwise.  
4. Movement is done by "F". Calculating next point using newpoint = previouspoint + Transpose[dv][[1]].
```mathematica
LSystem3D[axiom_, rules_, angle_, n_] := Module[{},
   chars = First@Characters@Nest[buildAxiom, {axiom, rules}, n];
   previouspoint = {0, 0, 0};
   dv = IdentityMatrix[3];
   lines = {};
   For[i = 1, i <= Length@chars, i++, {
     Switch[chars[[i]],
       "F", {
        newpoint = previouspoint + Transpose[dv][[1]];
        lines = Append[lines, {previouspoint, newpoint}];
        previouspoint = newpoint;},
       "+", dv = dv.mz[angle],
       "-", dv = dv.mz[-angle],
       "&", dv = dv.my[angle],
       "∧", dv = dv.my[-angle],
       "\\", dv = dv.mx[angle],
       "/", dv = dv.mx[-angle],
       "|", dv = dv.mz[Pi]
       ];
     }];
   lines
   ];
```

Defining 3D Hilbert curve path ( from the book [Algorithmic Beauty of Plants](http://algorithmicbotany.org/papers/#abop) ) and drawing lines.
```mathematica
axiom = "A"; angle = Pi/2;
rules = {
    "A" -> "B-F+CFC+F-D&F∧D-F+&&CFC+F+B//",
    "B" -> "A&F∧CFB∧F∧D∧∧-F-D∧|F∧B|FC\∧F∧A//",
    "C" -> "|D∧|F∧B-F+C∧F∧A&&FA&F∧C+F+B∧F\∧D//",
    "D" -> "|CFB-F+B|FA&F∧A&&FB-F+B|FC//"
   };
lines = LSystem3D[axiom, rules, angle, 3];
Graphics3D[Line[lines]]
```

![Image of 3D Hilbert Curve](/img/35.3DHilbertCurve.gif)


### 36. 2D Bracketed OL-systems [L-system](https://en.wikipedia.org/wiki/L-system) [Strahler number](https://en.wikipedia.org/wiki/Strahler_number)

Expanding LSystem function with brackets.  
"[": Push current state to stack.  
"]": Pop current state from stack.
```mathematica
Switch[chars[[i]],
  "[", {
   depth += 1;
   stack = Append[stack, {currentpoint, currentAngle, depth}];
   },
  "]", {
   {currentpoint, currentAngle, depth} = stack[[Length@stack]];
   stack = Drop[stack, -1];
   depth -= 1;
   }
  ];
```

Defining drawing function. Choosing style according to stack depth.
```mathematica
draw[{axiom_, rules_, angle_, startingpoint_, n_}] := Module[{},
   lines = LSystem[axiom, rules, angle, startingpoint, n];
   gLines = Table[
     {points, index} = lines[[i]];
     Graphics[Append[style[[index]], Line[points]]],
     {i, 1, Length@lines}]
   ];
   style = {
   {Thickness[0.01], Hue[.3, 1, .4]},
   {Thickness[0.005], Hue[.3, 1, .5]},
   {Thickness[0.004], Hue[.3, 1, .6]},
   {Thickness[0.003], Hue[.3, 1, .8]},
   {Thickness[0.002], Hue[.3, 1, .9]},
   {Thickness[0.001], Hue[.3, 1, 1]}
   };
```

Drawing weeds from book [Algorithmic Beauty of Plants](http://algorithmicbotany.org/papers/#abop).  
Recursion depth 6, Total 15625 lines, skipping over 40 to draw faster.
```mathematica
gLines = draw[{"F", {"F" -> "F[+F]F[-F][F]", "F" -> "FF-[-F+F+F]+[+F-F-F]", "F" -> "F[+F]F[-F]F"}, 22.5 Degree, {0, 0}, 6}];
Show[gLines]
```

![Image of 2D Bracketed OL-systems](/img/36.2DBracketedOLSystems.gif)

### 37. 3D Bracketed OL-systems [L-system](https://en.wikipedia.org/wiki/L-system) [Strahler number](https://en.wikipedia.org/wiki/Strahler_number)

Updating LSystem3D function to support brackets. Stack holds current point, 3D direction vector and stack depth.
```mathematica
"[", {
 depth += 1;
 stack = Append[stack, {currentpoint, dv, depth}];
 },
"]", {
 {currentpoint, dv, depth} = stack[[Length@stack]];
 stack = Drop[stack, -1];
 depth -= 1;
```

Defining some basic axiom and rule, angle Pi/4, recursion depth 4.
```mathematica
gLines = draw["F", {"F" -> "FF[+F][-F][/F][∧F][&F]"}, Pi/4, {0, 0, 0}, 4];
Show[gLines]
```

![Image of 3D Bracketed OL-systems](/img/37.3DBracketedOLSystems.gif)

### 38. Fourier Series [Mathologer](https://www.youtube.com/watch?v=qS4H6PEcCCA) [Stackexchange](https://mathematica.stackexchange.com/questions/171755/how-can-i-draw-a-homer-with-epicycloids) [Discrete Fourier Transform](http://mathworld.wolfram.com/DiscreteFourierTransform.html)

Defining function to turn text into points.
```mathematica
textToPoints[text_, skip_] := Module[{},
   img = Rasterize[Style[text, FontFamily -> "Lucida Handwriting"], ImageSize -> 500, RasterPadding -> 1];
   pts = {#2, -#1} & @@@ Position[ImageData[EdgeDetect[img]], 1, {2}];
   pts = Map[{0, 100} + # &, pts];
   shortest = Last@FindShortestTour@pts;
   pts = pts[[shortest]][[;; ;; skip]]
   ];
```

Defining function to create animation frames. Explanation in Mathologer video and Stackexchange post.
```mathematica
createFrames[pts_, m_] := Module[{},
   SetAttributes[toPt, Listable];
   toPt[z_] := ComplexExpand[{Re@z, Im@z}] // Chop;
   cf = Compile[{{z, _Complex, 1}}, Module[{n = Length@z}, 1/n*Table[Sum[z[[k]]*Exp[-I*i*k*2 Pi/n], {k, 1, n}], {i, -m, m}]]];
   z = pts[[All, 1]] + I*pts[[All, 2]];
   cn = cf[z];
   {f[t_], g[t_]} = Sum[cn[[j]]*Exp[I*(j - m - 1)*t], {j, 1, 2 m + 1}] // toPt;
   r = Abs[cn];
   theta = Arg[cn];
   index = {m + 1}~Join~Riffle[Range[m + 2, 2 m + 1], Reverse[Range[1, m]]];
   
   points[t_] = Accumulate@Table[cn[[j]]*Exp[I*(j - m - 1)*t], {j, index}] // toPt;
   circles[t_] = Table[Circle[points[t][[i - 1]], r[[index[[i]]]]], {i, 1, 2 m + 1}];
   
   range = {{0, 500}, {-150, 250}};
   ParallelTable[
    ParametricPlot[{f[s], g[s]}, {s, 0, t}, AspectRatio -> Automatic, 
     PlotStyle -> Red, PlotPoints -> 200, PlotRange -> range,
     Epilog -> {circles[t][[2 ;;]], Line[points[t]], 
       Point[points[t]]}, ImageSize -> Large], {t, Pi/36, 4*Pi, Pi/36}]
   ];
```

Using functions.  
1. Writing text to animate. '1' defines over how many points we skip.  
2. Points go to createFrames function. 300 is iteration count.
```mathematica
pts = textToPoints["FourierSeries", 1];
frames = createFrames[pts, 300];
ListAnimate[frames];
```

![Image of Fourier Series](/img/38.FourierSeries.gif)

### 39. Bezier curve [Wikipedia](https://en.wikipedia.org/wiki/B%C3%A9zier_curve)

Defining function to turn text into points.
```mathematica
textToPoints[text_, skip_] := Module[{},
   img = Rasterize[Style[text, FontFamily -> "Lucida Handwriting"], ImageSize -> 500];
   pts = {#2, -#1} & @@@ Position[ImageData[EdgeDetect[img]], 1, {2}];
   shortest = Last@FindShortestTour@pts;
   pts[[shortest]][[;; ;; skip]]
   ];
points = textToPoints["Bezier", 5];
```

Partitioning points in 3, offset 2. Calculating Bezier curve between them. 
```mathematica
pp = Partition[points, 3, 2];
curves = Table[{(1 - t)^2*pp[[i, 1]] + 2*t*(1 - t)*pp[[i, 2]] + t^2*pp[[i, 3]]}, {i, 1, Length@pp}];
```

Folding and animating.
```mathematica
folded = FoldList[List, plots];
frames = ParallelTable[Show[folded[[i]]], {i, 1, Length@folded}];
ListAnimate[frames]
```

![Image of Bezier curve](/img/39.BezierCurve.gif)

### 40. Turmite [Wikipedia](https://en.wikipedia.org/wiki/Turmite)

Defining Turmite function. New state and new color are defined according to rules matrix.
```mathematica
Turmite[{x_, y_, angle_, matrix_}] := Module[{},
   {color, state} = matrix[[x, y]];
   {newcolor, ruleangle, newstate} = rules[[state + 1, color + 1]];
   m = ReplacePart[matrix, {x, y} -> {newcolor, newstate}];
   {x + Cos[angle + ruleangle], y + Sin[angle + ruleangle], angle + ruleangle, m}
   ];
```

Grid drawing function.
```mathematica
draw[result_] := {
   {x0, y0, angle, matrix} = result;
   gGrid = Table[{
      {color, state} = matrix[[x, y]];
      color = Switch[color, 0, White, 1, Black];
      If[x == x0 && y == y0,
       Graphics[{EdgeForm[Thin], Red, Rectangle[{x0, y0}, {x0 + 1, y0 + 1}]}],
       Graphics[{EdgeForm[Thin], color, Rectangle[{x, y}, {x + 1, y + 1}]}]]
      }, {x, 1, Length@matrix}, {y, 1, Length@matrix}]
   };
```

Defining rules and running Turmite function recursively, step count 10000.
```mathematica
rules = {
   {{1, -Pi/2, 1}, {0, 0, 1}},
   {{1, 0, 0}, {0, Pi/2, 1}}};
result = NestList[Turmite, {50, 50, Pi/2, createMatrix[100]}, 10000];
Show[draw[result[[Length@result]]]];
```

![Image of Turmite](/img/40.Turmite.gif)
