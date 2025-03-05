F = nfinit(y^2 - 11);/*The field Q(sqrt(11))*/
A = alginit(F, [y - 3, -3]);/*Quaternion algebra (y-3,-3/F), which is split at one real place and has a>0 at that split place*/
X = afuchinit(A);/*Initialize and compute the fundamental domain*/
afuchfdom_latex(X, "example_fdom", 1, 1, 1, 1);/*Create (and compile) the fundamental domain in Latex, stored in ./plots/. Requires pdftex to compile it. Make the second last input a 0 if you don't want auto-compilation.*/
S = afuchsignature(X);/*Signature*/
P = afuchpresentation(X);/*Retrieve presentation*/
elts = afuchelts(X);
elt = algmulvec(A, elts, [1, 2, 7, 1, 2, 2, 1, 5, 7]);/*Generating an element of the group (happens to be hyperbolic)*/
word = afuchword(X, elt);/*Write it as a word in P[1].*/
geod = afuchgeodesic(X, elt);/*Finding the geodesic in the fundamental domain*/
afuchfdom_python(X, "example_fdom");/*Store the domain for fdomviewer.py*/
afuchgeodesic_python(X, geod, "example_geodesic");/*Store the geodesic for fdomviewer.py*/
fdomviewer("example_fdom example_geodesic");/*Open it (calls "python3 fdomviewer.py example_fdom example_geodesic" from the command line on Linux/Mac, and "cmd.exe /C start py fdomviewer.py example_fdom example_geodesic" on Windows. Alternatively, open up another terminal and type these commands in from this folder.).*/