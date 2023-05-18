print("\n\nType '?fdom' for help.\n\n");
addhelp(fdom, "This package can be used to compute fundamental domains for Arithmetic Fuchsian groups.\n For each subtopic ``P (p)'', call ?p to access a basic description and list of methods. Subtopics:");
parigp_version=version();
fdom_library=strprintf("./libfdom-%d-%d.so", parigp_version[1], parigp_version[2]);

\\fdom.c

\\SECTION 1: GEOMETRIC METHODS
	\\1: MATRIX ACTION ON GEOMETRY
	install(klein_act,"GGp",,fdom_library);
	addhelp(klein_act,"klein_act(M, z): returns the action of M on z, where we are working in the Klein model. M=[A, B] with |A|^2-|B|^2=1 acts on the unit disc model via the normal Mobius action of [A, B;conj(B), conj(A)].");

\\SECTION 3: QUATERNION ALGEBRA METHODS
	
	\\3: INITIALIZE ARITHMETIC FUCHSIAN GROUPS
	install(afuchinit,"GDGDGD1,L,p",,fdom_library);
	addhelp(afuchinit,"afuchinit(al, {O}, {type}, {flag=1}): initializes the arithmetic Fuchsian group in the algebra al with respect to the order O and of the given type. We work in the Klein model where p=Pi/8+0.5*I is sent to 0. The default order O is the stored maximal order in A, and the default type is 0. If flag = 1, also computes the fundamental domain. flag = 2 also computes the presentation and signature.");
	install(afuch_changep,"vGG",,fdom_library);
	addhelp(afuch_changep,"afuch_changep(X, p): changes the value of p, i.e. what is sent to 0 under the map from the upper half plane to the unit disc/Klein model. We recompute the fundamental domain and presentation if they were initialized.");

	\\3: ALGEBRA FUNDAMENTAL DOMAIN METHODS
	install(afucharea,"G",,fdom_library);
	addhelp(afucharea,"afucharea(X): returns the area of the computed fundamental domain of X.");
	install(afuchelts,"G",,fdom_library);
	addhelp(afuchelts,"afuchelts(X): returns the vector of elements giving the sides of the fundamental domain of X, which generate the group.");
	install(afuchfdom,"G",,fdom_library);
	addhelp(afuchfdom,"afuchfdom(X): returns the fundamental domain of X. The elements returned are with respect to the basis of O, so you must convert back if you want to use them in A.");
	install(afuchgeodesic,"GG",,fdom_library);
	addhelp(afuchgeodesic,"afuchgeodesic(X, g): computes the image of the closed geodesic associated to the hyperbolic element g in the fundamental domain (g belongs to X as well). The return is a vector with entries being [g, s1, s2, v1, v2, [a, b, c]], where each component runs from vertex v1 on side s1 to vertex v2 on side s2, which has equation ax+by=c=0 or 1. The components are listed in order.");
	install(afuchlist,"GGDGD1,L,",,fdom_library);
	addhelp(afuchlist,"afuchlist(F, Amin, {Amax}, {split=1}: given a totally real number field F (with variable not x), we find all possible quaternion algebras over F that are split at the unique real place given by split, for which the area of the fundamental domain is between Amin and Amax. If Amax is not passed, we go from 0 to Amin. The return is [[[a, b], area, rprimes]], where A=alginit(F, [a, b]) has area area, and rprimes is the multiset of primes lying above the finite ramified primes of A.");
	install(afuchpresentation,"G",,fdom_library);
	addhelp(afuchpresentation,"afuchpresentation(X): returns the presentation P of X. P[1] is the vector of generators, and P[2] is the vector of relations, where [1, -4, 3, 3] corresponds to P[1][1]*P[1][4]^-1*P[1][3]*P[3][3] being in the centre of A. P[3] tracks the sides of the fundamental domain in terms of the generators here, used to write an element as a word in these generators.");
	install(afuchsignature,"G",,fdom_library);
	addhelp(afuchsignature,"afuchsignature(X): returns the signature of X, i.e. [genus, [lengths of elliptic cycles], # of parabolic cycles].");
	install(afuchspair,"G",,fdom_library);
	addhelp(afuchspair,"afuchspair(X): returns S, the side pairing of X. The format is a Vecsmall, where side i is paired with side S[i]. Note that under some conventions, S[i]=i corresponds to there being an extra vertex at the midpoint of the edge (which is fixed), and the pairing is between the two new distinct sides formed.");
	install(afuchword,"GG",,fdom_library);
	addhelp(afuchword,"afuchword(X, g): writes g as a word in terms of the presentation of X. The format is a Vecsmall v, corresponding to the product of P[1][|v[i]|]^{sign(v[i])}, where P is the presentation. We do not check the relations of the presentation and eliminate their occurences, though I suspect that this non-trivial behaviour will occur rarely/never.");

	\\3: FINDING ELEMENTS
	install(afuchfindelts,"GD1,G,D1,L,DG",,fdom_library);
	addhelp(afuchfindelts,"afuchfindelts(X, {nm=1}, {N=1}, {C=default}: find N non-trivial elements of X of norm nm in the normalizer of the order O by solving Q_{z, 0}^nm(g)<=C for random points z. Elements found may be equal.");

	\\3: ALGEBRA HELPER METHODS
	install(algab,"G",,fdom_library);
	addhelp(algab,"algab(A): returns [a, b] where A=(a, b/F).");
	install(alg1ijktoalg,"GG",,fdom_library);
	addhelp(alg1ijktoalg,"alg1ijktoalg(A, g): returns what g=[e, f, g, h]=e+fi+gj+hk is in the algebraic representation.");
	install(alg1ijktobasis,"GG",,fdom_library);
	addhelp(alg1ijktobasis,"alg1ijktobasis(A, g): returns what g=[e, f, g, h]=e+fi+gj+hk is in the basis representation.");
	install(algalgto1ijk,"GG","algalgto1ijk",fdom_library);
	addhelp(algalgto1ijk,"algalgto1ijk(A, g): returns what g is in the 1ijk representation.");
	install(algbasisto1ijk,"GG",,fdom_library);
	addhelp(algbasisto1ijk,"algbasisto1ijk(A, g): returns what g is in the 1ijk representation.");
	install(algmulvec,"GGG",,fdom_library);
	addhelp(algmulvec,"algmulvec(A, G, L): returns G[L[1]]*G[L[2]]*...*G[L[n]]. If an index is negative, we take the inverse of that element.");
	install(algisorder,"iGG",,fdom_library);
	addhelp(algisorder,"algisorder(A, O): returns 1 if the lattice generated by the columns of O is an order, 0 if not.");
	install(algorderlevel,"GGD0,L,",,fdom_library);
	addhelp(algorderlevel,"algorderlevel(A, O, {flag=0}: returns the level of the order O in A as an ideal in the centre of A. If flag=1, this is in factored form.");
	
	\\4: FINCKE POHST FOR FLAG=2 WITH PRUNING
	install(qfminim_prune,"GGD1,L,p",,fdom_library);
	addhelp(qfminim_prune,"qfminim_prune(M, C, {prunetype=1}): Finds a set of vectors x solving x~*M*x<=C, and returns them as columns of a matrix. If prunetype=1, we use linear pruning (c.f. Schnorr-Euchner).");


\\fdom_extra.c

\\SECTION 1: VISUALIZATION
	install(afuchfdom_latex,"vGrD1,L,D1,L,D1,L,D1,L,",,fdom_library);
	addhelp(afuchfdom_latex,"afuchfdom_latex(X, filename, {model=1}, {boundcircle=1}, {compile=1}, {open=1}): writes the fundamental domain of X to a LaTeX document in plots/build/filename.tex. If model=0, use the upper half plane model, if model=1, use the unit disc model, and if model=2, use the Klein model. If boundcircle=0, does not print the bounding circle. If compile=1, compiles the document and moves it up to plots/build. If open=1, also opens the file (WSL only). Requires standalone, which can be found in texlive-latex-extra. NOTE: displaying in the Klein model is not suggested, as points are closer to the unit disc, and it does not show up very well.");
	install(afuchfdom_python,"vGr",,fdom_library);
	addhelp(afuchfdom_python,"afuchfdom_python(X, filename): writes the fundamental domain of X to a file ./fdoms/filename.dat, which can be read by the Python program fdomviewer to visualize the domain. Call 'fdomviewer.py filename' to open the Python application. Requires having both numpy and mathplotlib installed.");
	install(afuchgeodesic_python,"vGGr",,fdom_library);
	addhelp(afuchgeodesic_python,"afuchgeodesic_python(X, g, filename): writes the geodesic given by g in the fundamental domain of X to the file ./fdoms/filename.dat, which can be read by the Python program fdomviewer to visualize the domain. We can pass either g as an element of the algebra, or as the output of afuchgeodesic. call 'fdomviewer.py fdomname geodname' to see the geodesic in the fundamental domain.");
	install(fdomviewer,"vr",,fdom_library);
	addhelp(fdomviewer,"fdomviewer(files): assuming the user is using WSL with Python installed in Windows, calls 'fdomviewer.py files' in Python to launch the viewer. Files should be a string of space separated files, the first being the fundamental domain you want to visualize, and the rest being the geodesics.'");

\\SECTION 2: TUNING

	\\2: BEST C
	install(tune_Cn,"LGGD4,L,D20,L,p",,fdom_library);
	addhelp(tune_Cn,"tune_Cn(n, Cmin, Cmax, {testsperalg=4}, {tests=20}): For the degree n (between 1 and 9), we compute the fundamental domains for a range of algebras with C_n between Cmin and Cmax. We return the values of C_n and the total time taken for each one. This is used to determine the best value of C_n.");

default(parisize, "4096M");\\Must come at the end

