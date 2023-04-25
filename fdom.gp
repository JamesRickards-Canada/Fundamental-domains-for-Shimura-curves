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
	install(afuchinit,"GDGDGDGD1,L,p",,fdom_library);
	addhelp(afuchinit,"afuchinit(al, {O}, {type}, {p}, {flag=2}): initializes the arithmetic Fuchsian group in the algebra al with respect to the order O and of the given type. We work in the Klein model where p is an upper half plane point that is sent to 0. The default order O is the stored maximal order in A, the default type is 0, and the default value of p is Pi/8+0.5*I. If flag = 1, also computes the fundamental domain. flag = 2 also computes the presentation and signature.");

	\\3: ALGEBRA FUNDAMENTAL DOMAIN METHODS
	install(afuchfdom,"G",,fdom_library);
	addhelp(afuchfdom,"afuchfdom(X): returns the fundamental domain of X. The elements returned are with respect to the basis of O, so you must convert back if you want to use them in A.");
	install(afuchpresentation,"G",,fdom_library);
	addhelp(afuchpresentation,"afuchpresentation(X): returns the presentation P of X. P[1] is the vector of generators, and P[2] is the vector of relations, where [1, -4, 3, 3] corresponds to P[1][1]*P[1][4]^-1*P[1][3]*P[3][3] being in the centre of A. P[3] tracks the sides of the fundamental domain in terms of the generators here, used to write an element as a word in these generators.");
	install(afuchsignature,"G",,fdom_library);
	addhelp(afuchsignature,"afuchsignature(X): returns the signature of X, i.e. [genus, [lengths of elliptic cycles], # of parabolic cycles].");
	install(afuchredelt,"GDGD0,G,",,fdom_library);
	addhelp(afuchredelt,"afuchredelt(X, {g=id}, {z=0}: reduces gz to the fundamental domain, returning [g'g, g'gz, decomp], where g'gz is reduced, and g'=algmulvec(A, U[1], decomp). Except for a set of z of area 0, g'g should be trivial.");
	install(afuchword,"GG",,fdom_library);
	addhelp(afuchword,"afuchword(X, g): writes g as a word in terms of the presentation of X. We do not check the relations of the presentation and eliminate their occurences, though I suspect that this non-trivial behaviour will occur rarely/never.");

	\\3: FINDING ELEMENTS
	install(afuchfindelts,"GDGDGD1,L,",,fdom_library);
	addhelp(afuchfindelts,"afuchfindelts(X, {z}, {C}, {maxelts=1}: finds the elements of norm 1 for which Q_{z, 0}(g)<=C, which will happen if gz is close to 0. Returns at most maxelts, unless this is set to 0 when it returns all elements found. If C is not passed, we use the default value which attemps to minimize the expected time to find an element (when we run it with random points, stopping once the first non-trivial element is found). If z is not passed we pick a uniformly random point from a ball of radius large enough so that we expect it to be random with respect to the fundamental domain.");

	\\3: ALGEBRA HELPER METHODS
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
	addhelp(afuchfdom_latex,"Inputs: U, filename, {model=1}, {boundcircle=1}, {compile=1}, {open=1}.\n Writes the fundamental domain U to a LaTeX document in plots/build/filename.tex. If model=0, use the upper half plane model, if model=1, use the unit disc model, and if model=2, use the Klein model. If boundcircle=0, does not print the bounding circle. If compile=1, compiles the document and moves it up to plots/build. If open=1, also opens the file (WSL only). Requires standalone, which can be found in texlive-latex-extra. NOTE: displaying in the Klein model is not suggested, as points are closer to the unit disc, and it does not show up very well.");

\\SECTION 2: TUNING

	\\2: BEST C
	install(tune_Cn,"LGGD4,L,D20,L,p",,fdom_library);
	addhelp(tune_Cn,"tune_Cn(n, Cmin, Cmax, {testsperalg=4}, {tests=20}): For the degree n (between 1 and 9), we compute the fundamental domains for a range of algebras with C_n between Cmin and Cmax. We return the values of C_n and the total time taken for each one. This is used to determine the best value of C_n.");
	install(tune_Nelts,"lGGD100,L,p",,fdom_library);
	addhelp(tune_Nelts,"tune_Nelts(X, C, {nelts=100}): Computes the time taken to find nelts elements with the given value of C.");

default(parisize, "4096M");\\Must come at the end

