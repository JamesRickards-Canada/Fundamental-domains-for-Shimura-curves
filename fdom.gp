print("\n\nType '?fdom' for help.\n\n");
addhelp(fdom, "This package can be used to compute fundamental domains for congruence Arithmetic Fuchsian groups.\n Installed methods:\nklein_act.\nafuchinit, afuchnewp, afuchnewtype, afuchmoreprec.\nafuchalg, afucharea, afuchelliptic, afuchelts, afuchelttype, afuchgeodesic, afuchlist, afuchmakefdom, afuchnormalizernorms, afuchorder, afuchpresentation, afuchsides, afuchsignature, afuchspair, afuchvertices, afuchword.\nafuchfindoneelt\nalgab, alg1ijktoalg, alg1ijktobasis, algalgto1ijk, algbasisto1ijk, algmulvec, algisorder, algorderalgtoorder, algordertoalgorder, algorderlevel, algreduceddisc.\nqfminim_prune\nafuchfdom_latex, afuchfdom_python, afuchgeodesic_python, fdomviewer\nalgeichlerorder\ntune_Cn");
parigp_version = version();
fdom_library = strprintf("./libfdom-%d-%d-%d.so", parigp_version[1], parigp_version[2], parigp_version[3]);

/*fdom.c*/

/*SECTION 1: GEOMETRIC METHODS*/
  /*1: MATRIX ACTION ON GEOMETRY*/
  install(klein_act,"GGp",,fdom_library);
  addhelp(klein_act,"klein_act(M, z): returns the action of M on z, where we are working in the Klein model. M=[A, B] with |A|^2-|B|^2=1 acts on the unit disc model via the normal Mobius action of [A, B;conj(B), conj(A)].");

/*SECTION 3: QUATERNION ALGEBRA METHODS*/
    
  /*3: INITIALIZE ARITHMETIC FUCHSIAN GROUPS*/
  install(afuchinit,"GDGDGD1,L,p");
  addhelp(afuchinit,"afuchinit(al, {O}, {type}, {flag=1}): initializes the arithmetic Fuchsian group in the algebra al with respect to the order O and of the given type. We work in the Klein model where p=Pi/8+0.5*I is sent to 0. The default order O is the stored maximal order in A, and the default type is 0. type=0 means O^1, type=1 means totally positive unit norm, type=2 is AL(O)^+, and type=3 is the whole positive normalizer. If flag = 1, also computes the fundamental domain and presentation. To compute a fundamental domain for a general congruence arithmetic Fuchsian group, initialize with type=3 and then use afuch_newtype.");
  install(afuchnewp,"GG");
  addhelp(afuchnewp,"afuchnewp(X, p): returns the Fuchsian group with the changed the value of p, i.e. what is sent to 0 under the map from the upper half plane to the unit disc/Klein model. We also recompute the fundamental domain and presentation if they were initialized.");
  install(afuchnewtype,"GG");
  addhelp(afuchnewtype,"afuchnewtype(X, type): returns a the Fuchsian group X but we change the type. If the fundamental domain was already computed with type=3, this is very efficient, as we do not have to search for new generators. In order to input a subgroup Gamma between O^1 and N_{A^x}^+(O), consider S=concat(afuchnormalizernorms(X)). We know that N_{A^x}^+(O)/O^1=(Z/2Z)^#S, with norms of generators corresponding to the entries of S. A subgroup of this can be specified by a matrix with #S rows, where each column consists of 0's/1's, indicating the generating set of the subgroup. Thus the 0 matrix gets you O^1, and any matrix with rank #S gets the full positive normalizer.");
  install(afuchmoreprec,"GD1,L,");
  addhelp(afuchmoreprec,"afuchmoreprec(X, inc=1): returns the Fuchsian group X with precision increased by inc steps, i.e. inc*DEFAULTPREC more precision.");

  /*3: ALGEBRA FUNDAMENTAL DOMAIN METHODS*/
  install(afuchalg,"G");
  addhelp(afuchalg,"afuchalg(X): retrieves the stored algebra in X.");
  install(afucharea,"G");
  addhelp(afucharea,"afucharea(X): retrieves the area of the computed fundamental domain of X.");
  install(afuchelliptic,"G");
  addhelp(afuchelliptic,"afuchelliptic(X): retrieves the elliptic elements of X.");
  install(afuchelts,"G");
  addhelp(afuchelts,"afuchelts(X): retrieves the vector of elements giving the sides of the stored fundamental domain of X, which generate the group.");
  install(afuchelttype,"iGG");
  addhelp(afuchelttype,"afuchelttype(X, g): returns 1 if g is hyperbolic, 0 if g is parabolic, and -1 if g is elliptic.");
  install(afuchgeodesic,"GG");
  addhelp(afuchgeodesic,"afuchgeodesic(X, g): computes the image of the closed geodesic associated to the hyperbolic element g in the fundamental domain (g belongs to X as well). The return is a vector with entries being [g, s1, s2, v1, v2, [a, b, c]], where each component runs from vertex v1 on side s1 to vertex v2 on side s2, which has equation ax+by=c=0 or 1. The components are listed in order.");
  install(afuchlist,"GGDGD1,L,");
  addhelp(afuchlist,"afuchlist(F, Amin, {Amax}, {split=1}: given a totally real number field F (with variable not x), we find all possible quaternion algebras over F that are split at the unique real place given by split, for which the area of the fundamental domain is between Amin and Amax. If Amax is not passed, we go from 0 to Amin. The return is [[[a, b], area, rprimes]], where A=alginit(F, [a, b]) has area area, and rprimes is the multiset of primes lying above the finite ramified primes of A.");
  install(afuchmakefdom,"G");
  addhelp(afuchmakefdom,"afuchmakefdom(X): returns X with the fundamental domain and presentation computed.");
  install(afuchminimalcycles,"G");
  addhelp(afuchminimalcycles,"afuchminimalcycles(X): computes the minimal cycles in X, returning [cycles, types], where cycles[i] has type types[i]. Type 0=parabolic, 1=accidental, m>=2=elliptic of order m. It is returned with the types sorted, i.e. parabolic cycles first, then accidental, then elliptic.");
  install(afuchnormalizernorms,"G");
  addhelp(afuchnormalizernorms,"afuchnormalizernorms(X): let N_{B^x}(O)^+ be the norm totally positive normalizer of the order O. The quotient of this group by F^xO^1 is a 2-group. This method returns the norms of a generating set for this group. The format is [unit, AL, rest], where unit is the set of norms that are units, AL is the set of norms appearing in the Atkin-Lehner group, and rest is the remaining norms in the normalizer, coming from the 2-torsion of the class group. We expect the norms found to be linearly independant with respect to the normalizer quotient, but this is not necessarily guaranteed.");
  install(afuchorder,"G");
  addhelp(afuchorder,"afuchorder(X): retrieves the stored order in X.");
  install(afuchpresentation,"G");
  addhelp(afuchpresentation,"afuchpresentation(X): retrieves the stored presentation P in X. P[1] is the vector of generators, and P[2] is the vector of relations, where [1, -4, 3, 3] corresponds to P[1][1]*P[1][4]^-1*P[1][3]*P[3][3] being in the centre of A.");
  install(afuchsides,"G");
  addhelp(afuchsides,"afuchsides(X): retrieves the sides of the fundamental domain. The format of a side is [a, b, r], where the equation for the side is ax+by=1 in the Klein model, (x-a)^2+(y-b)^2=r^2=a^2+b^2-1 in the unit disc model.");
  install(afuchsignature,"G");
  addhelp(afuchsignature,"afuchsignature(X): returns the signature of X, i.e. [genus, [lengths of elliptic cycles], # of parabolic cycles].");
  install(afuchspair,"G");
  addhelp(afuchspair,"afuchspair(X): retrieves S, the side pairing of X. The format is a Vecsmall, where side i is paired with side S[i]. Note that under some conventions, S[i]=i corresponds to there being an extra vertex at the midpoint of the edge (which is fixed), and the pairing is between the two new distinct sides formed.");
  install(afuchvertices,"GD0,L,");
  addhelp(afuchvertices,"afuchvertices(X, {model=0}): returns the vertices of the fundamental domain, where side i corresponds to going from vertex i-1 to i counterclockwise around the unit disc. If model=0 they are in the Klein model, and if model=1 they are in the unit disc model.");
  install(afuchword,"GG");
  addhelp(afuchword,"afuchword(X, g): writes g as a word in terms of the presentation of X. The format is a Vecsmall v, corresponding to the product of P[1][|v[i]|]^{sign(v[i])}, where P is the presentation. We do not check the relations of the presentation and eliminate their occurences, though I suspect that this non-trivial behaviour will occur rarely. If your element g is ''too big'', this will internally recompute X to more and more precision. If you do this with multiple elements, it is more efficient to increase X's precision first.");

  /*3: FINDING ELEMENTS*/
  install(afuchfindoneelt,"GD1,G,DG");
  addhelp(afuchfindoneelt,"afuchfindoneelt(X, {nm=1}, {C=default}: returns one non-trivial element of X of norm nm in the normalizer of the order O by solving Q_{z, 0}^nm(g)<=C for random points z. The element is primitive in O, so if no such element exists, this will be an infinite loop.");

  /*3: ALGEBRA HELPER METHODS*/
  install(algab,"G");
  addhelp(algab,"algab(A): returns [a, b] where A=(a, b/F).");
  install(alg1ijktoalg,"GG");
  addhelp(alg1ijktoalg,"alg1ijktoalg(A, g): returns what g=[e, f, g, h]=e+fi+gj+hk is in the algebraic representation.");
  install(alg1ijktobasis,"GG");
  addhelp(alg1ijktobasis,"alg1ijktobasis(A, g): returns what g=[e, f, g, h]=e+fi+gj+hk is in the basis representation.");
  install(algalgto1ijk,"GG");
  addhelp(algalgto1ijk,"algalgto1ijk(A, g): returns what g is in the 1ijk representation.");
  install(algbasisto1ijk,"GG");
  addhelp(algbasisto1ijk,"algbasisto1ijk(A, g): returns what g is in the 1ijk representation.");
  install(algmulvec,"GGG");
  addhelp(algmulvec,"algmulvec(A, G, L): returns G[L[1]]*G[L[2]]*...*G[L[n]]. If an index is negative, we take the inverse of that element.");
  install(algisorder,"iGG",);
  addhelp(algisorder,"algisorder(A, O): returns 1 if the lattice generated by the columns of O is an order, 0 if not.");
  install(algorderalgtoorder,"GG");
  addhelp(algorderalgtoorder,"algorderalgtoorder(A, Oalg): given an order of A expressed as Oalg, a vector of elements of A in algebraic form, this returns the matrix whose columns are the basis in basis form.");
  install(algordertoorderalg,"GG");
  addhelp(algordertoorderalg,"algordertoorderalg(A, O): given an order O of A expressed as a matrix, whose columns form the basis of the order, we convert it to algebraic form, i.e. return a vector whose entries are the basis in algebraic form.");
  install(algorderlevel,"GGD0,L,");
  addhelp(algorderlevel,"algorderlevel(A, O, {flag=0}: returns the level of the order O in A as an ideal in the centre of A. If flag=1, this is in factored form.");
  install(algreduceddisc,"G");
  addhelp(algreduceddisc,"algreduceddisc(A): returns the reduced discriminant of A as an ideal in its centre.");
  
  /*4: FINCKE POHST FOR FLAG=2 WITH PRUNING*/
  install(qfminim_prune,"GGD1,L,p");
  addhelp(qfminim_prune,"qfminim_prune(M, C, {prunetype=1}): Finds a set of vectors x solving x~*M*x<=C, and returns them as columns of a matrix. If prunetype=1, we use linear pruning (c.f. Schnorr-Euchner).");


/*fdom_extra.c*/

/*SECTION 1: VISUALIZATION*/
  install(afuchfdom_latex,"vGrD1,L,D1,L,D1,L,D1,L,");
  addhelp(afuchfdom_latex,"afuchfdom_latex(X, filename, {model=1}, {boundcircle=1}, {compile=1}, {open=1}): writes the fundamental domain of X to a LaTeX document in ./plots/build/filename.tex. If model=0, use the Klein model, if model=1, use the unit disc model, and if model=2, use the upper half plane model. If boundcircle=0, does not print the bounding circle. If compile=1, compiles the document and moves it up to ./plots/filename.pdf. If open=1, also opens the file (WSL only). Requires standalone, which can be found in texlive-latex-extra. NOTE: displaying in the Klein model is not suggested, as points are closer to the unit disc, and it does not show up very well.");
  install(afuchfdom_python,"vGr");
  addhelp(afuchfdom_python,"afuchfdom_python(X, filename): writes the fundamental domain of X to a file ./fdoms/filename.dat, which can be read by the Python program fdomviewer to visualize the domain. Call 'fdomviewer.py filename' to open the Python application. Requires having both numpy and mathplotlib installed.");
  install(afuchgeodesic_python,"vGGr");
  addhelp(afuchgeodesic_python,"afuchgeodesic_python(X, g, filename): writes the geodesic given by g in the fundamental domain of X to the file ./fdoms/filename.dat, which can be read by the Python program fdomviewer to visualize the domain. We can pass either g as an element of the algebra, or as the output of afuchgeodesic. call 'fdomviewer.py fdomname geodname' to see the geodesic in the fundamental domain.");
  install(fdomviewer,"vr");
  addhelp(fdomviewer,"fdomviewer(files): assuming the user is using WSL with Python installed in Windows, calls 'fdomviewer.py files' in Python to launch the viewer. Files should be a string of space separated files, the first being the fundamental domain you want to visualize, and the rest being the geodesics.'");

/*SECTION 2: EICHLER ORDERS*/
  install(algeichlerorder,"GG");
  addhelp(algeichlerorder,"algeichlerorder(A, I): returns an Eichler order of level I in A the stored maximal order of A. This uses Magma, so requires Magma to be installed. Creates and modifies the files 'fdom_make_Eichler.m' and 'fdom_make_Eichler_output.dat', and the algebra A cannot use the variables i, j, or k anywhere.");

/*SECTION 3: TESTING AND TUNING*/
  install(afuchcheck,"lG");
  addhelp(afuchcheck,"afuchcheck(X):runs a series of checks on X with the fundamental domain and presentation initialized. Returns 0 if all passed, and something non-zero else. These return codes are:\n\t1: signature area formula does not match computed area;\n\t2: presentation has wrong number of generators;\n\t3: presentation has wrong number of relations;\n\t4: one of the relations fails;\n\t5: one of the side pairing element relations fails;\n\t6: afuchfdomword fails on a random element (15 random elements tested).");
  install(tune_Cn,"LGGD4,L,D20,L,p");
  addhelp(tune_Cn,"tune_Cn(n, Cmin, Cmax, {testsperalg=4}, {tests=20}): For the degree n (between 1 and 9), we compute the fundamental domains for a range of algebras with C_n between Cmin and Cmax. We return the values of C_n and the total time taken for each one. This is used to determine the best value of C_n.");

default(parisize, "4096M");