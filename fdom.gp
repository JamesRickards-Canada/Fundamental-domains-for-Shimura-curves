print("\n\nType '?fdom' for help.\n\n");
addhelp(fdom, "This package can be used to compute fundamental domains for Arithmetic Fuchsian groups.\n For each subtopic ``P (p)'', call ?p to access a basic description and list of methods. Subtopics:");
parigp_version=version();
fdom_library=strprintf("./libfdom-%d-%d.so", parigp_version[1], parigp_version[2]);

\\SECTION 1: GEOMETRIC METHODS
	\\addhelp(geo,"These methods deal with geometry. Available methods:\n hdiscrandom, hdist, mat_eval");

	install(hdiscrandom,"Gp",,fdom_library);
	addhelp(hdiscrandom,"Input R, a positive real number.\n Returns a random point in the ball of radius R centred at 0 in the unit disc model of the hyperbolic plane.");
	install(hdist,"GGD1,L,p",,fdom_library);
	addhelp(hdist,"Inputs z1, z2, {flag=1}: complex numbers in the upper half plane z1 and z2, and flag=0, 1.\n Returns the hyperbolic distance between z1 and z2. If flag=0 we use the upper half plane model, and if flag=1 we use the unit disc model.");


\\SECTION 3: QUATERNION ALGEBRA METHODS
	
	\\3: INITIALIZE ARITHMETIC FUCHSIAN GROUPS
	install(afuchinit,"GGGDGp",,fdom_library);
	addhelp(afuchinit,"afuchinit(al, O, type, {p}): initializes the arithmetic Fuchsian group in the algebra al with respect to the order O and of the given type. We work in the Klein model where p is an upper half plane point that is sent to 0. The default value of p is 0.5+Pi/8*I.");

	\\3: ALGEBRA FUNDAMENTAL DOMAIN METHODS
	install(afuchnormbound,"GGp",,fdom_library);
	addhelp(afuchnormbound,"afuchnormbound(X, G): computes the normalized boundary of the set of elements in G, where X is an initialized arithmetic Fuchsian group. Elements of G may not have norm 0.");

	\\3: ALGEBRA HELPER METHODS
	install(algramifiedplacesf,"G",,fdom_library);
	addhelp(algramifiedplacesf,"algramifiedplacesf(al): vector of the finite places of the center of al that ramify in al. Each place is described as a prime ideal.");

default(parisize, "4096M");\\Must come at the end


