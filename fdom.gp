print("\n\nType '?fdom' for help.\n\n");
addhelp(fdom, "This package can be used to compute fundamental domains for Shimura curves (the PARI code can be easily adapted to compute fundamental domains for any discrete subgroup of PSL(2, R)).\n For each subtopic ``P (p)'', call ?p to access a basic description and list of methods. Subtopics:\n Geometry (geo)\n Visualizing fundamental domains with Python (vfd)\n Quaternion methods (quat)");
parigp_version=version();
fdom_library=strprintf("./libfdom-%d-%d.so", parigp_version[1], parigp_version[2]);

\\GEOMETRY
	\\addhelp(geo,"These methods deal with geometry. Available methods:\n hdiscrandom, hdist, mat_eval");

	install("hdiscrandom","Gp",,fdom_library);
	addhelp(hdiscrandom,"Input R, a positive real number.\n Returns a random point in the ball of radius R centred at 0 in the unit disc model of the hyperbolic plane.");
	install("hdist","GGD1,L,p",,fdom_library);
	addhelp(hdist,"Inputs z1, z2, {flag=1}: complex numbers in the upper half plane z1 and z2, and flag=0, 1.\n Returns the hyperbolic distance between z1 and z2. If flag=0 we use the upper half plane model, and if flag=1 we use the unit disc model.");

default(parisize, "4096M");\\Must come at the end