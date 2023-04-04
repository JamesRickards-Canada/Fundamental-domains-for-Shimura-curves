print("\n\nType '?fdom' for help.\n\n");
addhelp(fdom, "This package can be used to compute fundamental domains for Arithmetic Fuchsian groups.\n For each subtopic ``P (p)'', call ?p to access a basic description and list of methods. Subtopics:");
parigp_version=version();
fdom_library=strprintf("./libfdom-%d-%d.so", parigp_version[1], parigp_version[2]);

\\SECTION 1: GEOMETRIC METHODS
	\\1: MATRIX ACTION ON GEOMETRY
	install(klein_act,"GGp",,fdom_library);
	addhelp(klein_act,"klein_act(M, z): returns the action of M on z, where we are working in the Klein model. M=[A, B] with |A|^2-|B|^2=1 acts on the unit disc model via the normal Mobius action of [A, B;conj(B), conj(A)].");

	\\1: DISTANCES/AREAS
	install(hdiscrandom,"Gp",,fdom_library);
	addhelp(hdiscrandom,"Input R, a positive real number.\n Returns a random point in the ball of radius R centred at 0 in the unit disc model of the hyperbolic plane.");
	install(hdist,"GGD2,L,p",,fdom_library);
	addhelp(hdist,"hdist(z1, z2, {flag=2}): returns the hyperbolic distance between z1 and z2. If flag=0, uses the upper half plane model. If flag=1, uses the unit disc model. If flag=2, uses the Klein model.");


\\SECTION 3: QUATERNION ALGEBRA METHODS
	
	\\3: INITIALIZE ARITHMETIC FUCHSIAN GROUPS
	install(afuchinit,"GGGDGp",,fdom_library);
	addhelp(afuchinit,"afuchinit(al, O, type, {p}): initializes the arithmetic Fuchsian group in the algebra al with respect to the order O and of the given type. We work in the Klein model where p is an upper half plane point that is sent to 0. The default value of p is Pi/8+0.5*I.");

	\\3: ALGEBRA FUNDAMENTAL DOMAIN METHODS
	install(afuchicirc,"GG",,fdom_library);
	addhelp(afuchicirc,"afuchicirc(X, g): returns the isometric circle of g, an element of non-zero norm.");
	install(afuchklein,"GG",,fdom_library);
	addhelp(afuchklein,"afuchklein(X, g): returns the vector giving the action of g on the Klein model, which can be supplied to klein_act.");
	install(afuchnormbasis,"GG",,fdom_library);
	addhelp(afuchnormbasis,"afuchnormbasis(X, G): computes the normalized basis of the set of elements of G, where X is an initialized arithmetic Fuchsian group. Elements of G may not have norm 0.");
	install(afuchnormbound,"GG",,fdom_library);
	addhelp(afuchnormbound,"afuchnormbound(X, G): computes the normalized boundary of the set of elements in G, where X is an initialized arithmetic Fuchsian group. Elements of G may not have norm 0.");
	install(afuchnormbound_append,"GGG",,fdom_library);
	addhelp(afuchnormbound_append,"afuchnormbound_append(X, U, G): given a non-trivial normalized boundary U, this computes the new normalized boundary obtained by adding in the isometric circles for all elements in G. More efficient than calling afuchnormbound(X, concat(U[1], G)), as we use the already computed U to optimize the computation.");
	install(afuchredelt,"GGDGD0,G,",,fdom_library);
	addhelp(afuchredelt,"afuchredelt(X, U, {g=id}, {z=0}: reduces gz to the normalized boundary U, returning [g'g, g'gz, decomp], where g'gz is reduced, and g'=algmulvec(A, U[1], decomp).");

	\\3: ALGEBRA HELPER METHODS
	install(algmulvec,"GGG",,fdom_library);
	addhelp(algmulvec,"algmulvec(A, G, L): Returns G[L[1]]*G[L[2]]*...*G[L[n]]. If an index is negative, we take the inverse of that element.");
	install(algramifiedplacesf,"G",,fdom_library);
	addhelp(algramifiedplacesf,"algramifiedplacesf(al): vector of the finite places of the center of al that ramify in al. Each place is described as a prime ideal.");

default(parisize, "4096M");\\Must come at the end


