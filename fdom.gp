print("\n\nType '?fdom' for help.\n\n");
addhelp(fdom, "This package can be used to compute fundamental domains for Arithmetic Fuchsian groups.\n For each subtopic ``P (p)'', call ?p to access a basic description and list of methods. Subtopics:");
parigp_version=version();
fdom_library=strprintf("./libfdom-%d-%d.so", parigp_version[1], parigp_version[2]);

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
	install(afuchsignature,"G",,fdom_library);
	addhelp(afuchsignature,"afuchsignature(X): returns the signature of X, i.e. [genus, [lengths of elliptic cycles], # of parabolic cycles].");
	install(afuchredelt,"GGDGD0,G,",,fdom_library);
	addhelp(afuchredelt,"afuchredelt(X, U, {g=id}, {z=0}: reduces gz to the normalized boundary U, returning [g'g, g'gz, decomp], where g'gz is reduced, and g'=algmulvec(A, U[1], decomp).");

	\\3: FINDING ELEMENTS
	install(afuchfindelts,"GDGDGD1,L,",,fdom_library);
	addhelp(afuchfindelts,"afuchfindelts(X, {z}, {C}, {maxelts=1}: Finds the elements of norm 1 for which Q_{z, 0}(g)<=C, which will happen if gz is close to 0. Returns at most maxelts, unless this is set to 0 when it returns all elements found. If C is not passed, we use the default value which attemps to minimize the expected time to find an element (when we run it with random points, stopping once the first non-trivial element is found). If z is not passed we pick a uniformly random point from a ball of radius large enough so that we expect it to be random with respect to the fundamental domain.");

	\\3: ALGEBRA HELPER METHODS
	install(algmulvec,"GGG",,fdom_library);
	addhelp(algmulvec,"algmulvec(A, G, L): Returns G[L[1]]*G[L[2]]*...*G[L[n]]. If an index is negative, we take the inverse of that element.");

default(parisize, "4096M");\\Must come at the end


