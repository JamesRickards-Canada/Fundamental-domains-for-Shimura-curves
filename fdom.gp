print("\n\nType '?fdom' for help.\n\n");
addhelp(fdom, "This package can be used to compute fundamental domains for Shimura curves (the PARI code can be easily adapted to compute fundamental domains for any discrete subgroup of PSL(2, R)).\n For each subtopic ``P (p)'', call ?p to access a basic description and list of methods. Subtopics:\n Geometry (geo)\n Visualizing fundamental domains with Python (vfd)\n Quaternion methods (quat)");
default(help, "gphelp -detex");

\\GEOMETRY
	addhelp(geo,"These methods deal with geometry. Available methods:\n hdist, hdist_ud, randompoint_ud");
	
	install("hdist","GGp","hdist","./libfdom.so");
	addhelp(hdist,"Inputs z1, z2 complex numbers in the upper half plane.\n Returns the hyperbolic distance between z1 and z2.");
	install("hdist_ud","GGp","hdist_ud","./libfdom.so");
	addhelp(hdist_ud,"Inputs z1, z2 complex numbers inside the unit disc.\n Returns the hyperbolic distance between z1 and z2 in the unit disc model.");
	install("mat_eval","GG","mat_eval","./libfdom.so");
	addhelp(mat_eval, "Inputs M, x; M a matrix, and x number.\n Returns Mx with M acting via Mobius transformation. x=+/-oo is allowed.");
	install("randompoint_ud","Gp","randompoint_ud","./libfdom.so");
	addhelp(randompoint_ud,"Input R, a positive real number.\n Returns a random point in the ball of radius R centred at 0 in the unit disc model of the hyperbolic plane.");

\\Visualization
	addhelp(vfd,"These methods allow one to save fundamental domains and geodesics, and view them with a Python program. Available methods:\n python_plotviewer, python_printarcs, python_printfdom.");
	
	install("python_plotviewer","vr","python_plotviewer","./libfdom.so");
	addhelp(python_plotviewer,"Input S: string denoting the file names of data fundamental domains/geodesics.\n Launches the python file fdviewer.py to view the domain/geodesics. Enter the files separated by spaces (they must be stored in the sub-folder 'fdoms').");
	install("python_printarcs","vGrD0,L,Drp","python_printarcs","./libfdom.so");
	addhelp(python_printarcs,"Inputs arcs, filename, {view=0}, {extrainput=NULL}: a set of arcs arcs, string filename, view=0, 1, and extrainput=NULL or a string.\n Prints the arcs specified by arcs to the file fdoms/filename.dat, ready for plotviewer. If view=1, calls plotviewer with the additional input of extrainput if you want to include other arcs/fundamental domains.");
	install("python_printfdom","vGrp","python_printfdom","./libfdom.so");
	addhelp(python_printfdom,"Input U, filename: fundamental domain U, string filename.\n Prints the fundamental domain U to the file fdoms/filename.dat, ready for the plot viewer. The filename must start with 'fd' to work properly.")


\\Quaternion methods
	addhelp(quat,"These methods allow for the computation of fundamental domains for Eichler orders in quaternion algebras split at one real place. Available methods:\n algfdomarea, algfdom, algfdomminimalcycles, algfdompresentation, algfdomreduce, algfdomrootgeodesic, algfdomsignature, algfromnormdisc, algmulvec, algramifiedplacesf, algnormalizedbasis, algnormalizedboundary, algshimura, algsmallnorm1elts, algswapab.");

	install("algfdomarea","Gp","algfdomarea","./libfdom.so");
	addhelp(algfdomarea,"Input A, a quaternion algebra split at one real place.\n Returns the area of the fundamental domain associated to the group of units of norm 1 in the order stored in A. Requires the computation of the zeta_K(2) (Dedekind zeta function for the centre), which may require some calls to allocatemem() if K is ``large''.");
	install("algfdom","GGD1,L,D0,G,D0,G,p","algfdom","./libfdom.so");
	addhelp(algfdom,"Inputs A, p, {dispprogress=1}, {area=0}, {ANRdata=0}: quaternion algebra A split at one real place, upper half plane point p, dispprogress=0,1, ANRdata=0 or a length 5 vector, area=0 or the area of the fundamental domain.\n Computes and returns the fundamental domain for the group of units of norm 1 in A. We use the unit disc model, which the upper half plane is mapped to via p->0. p need to NOT be a fixed point of this group under the standard embedding into PSL(2, R) (p=I/2 is safe for quaternion algebras over Q). If area is passed in, this method will not re-compute it, which can save a little bit of time. ANRdata is a technical entry, and can be passed in to specify some/all constants used in the enumeration of Page (they greatly affect the running time, but the optimal choices are not totally clear).");
	install("algfdomminimalcycles","GGp","algfdomminimalcycles","./libfdom.so");
	addhelp(algfdomminimalcycles,"Inputs A, U: quaternion algebra A split at one real place with fundamental domain U.\n Returns the set of minimal cycles of the side pairing. The format is [cycles, types], where an element of cycle is a vecsmall [i1,i2,...,in] so that the cycle is v_i1, v_i2, ..., v_in. cycle[i] has type types[i], where type 0=parabolic, 1=accidental, m>=2=elliptic of order m. The vecsmall types is sorted from small to large.");
	install("algfdompresentation","GGp","algfdompresentation","./libfdom.so");
	addhelp(algfdompresentation,"Inputs A, U: quaternion algebra A split at one real place with fundamental domain U.\n Returns the group presentation of the fundamental domain U. The return is a vector, where the 1st element is the list of indices of the generators (with respect to U[1]), 2nd element is the vector of relations, whose ith element is a relation of the form [indices, powers], where indices and powers are vecsmall's. If indices=[i1,i2,...,ik] and powers=[p1,p2,...,pk], then this corresponds to g_{i1}^p1*...*g_{ik}^{pk}=1.");
	install("algfdomreduce","GGGD0,G,p","algfdomreduce","./libfdom.so");
	addhelp(algfdomreduce,"Inputs A, U, g, {z=0}: quaternion algebra A split at one real place, the fundamental domain U of the group of units of norm 1 of the order in A, an element g of this group, and a point z in the unit disc.\n Returns the triple [gammabar, delta, decomp], where gammabar=delta*g is (G,z)-reduced (i.e. distance between gammabar*z and 0 is less than or equal to the distance between g'*gammabar*z for all g' in G), and decomp is the vecsmall [i1, i2, ..., in] with delta=G[i1]*G[i2]*...*G[in]. If z=0, then gammabar=+/-1.");
	install("algfdomrootgeodesic","GGGp","algfdomrootgeodesic","./libfdom.so");
	addhelp(algfdomrootgeodesic,"Inputs A, U, g: quaternion algebra A split at one real place, the fundamental domain U of the group of units of norm 1 of the order in A, an element g of this group.\n Returns the root geodesic of g in the fundamental domain. The format is [g's, circle arcs, vecsmall(sides hit), vecsmall(sides emenating from)].");
	install("algfdomsignature","GGp","algfdomsignature","./libfdom.so");
	addhelp(algfdomsignature,"Inputs A, U: quaternion algebra A split at one real place with fundamental domain U.\n Returns the signature of the algebra A. The format is [g, V, s], where g is the genus, V=[m1,m2,...,mt] (vecsmall) are the orders of the elliptic cycles (all >=2), and s is the number of parabolic cycles. The signature is normally written as (g;m1,m2,...,mt;s), and the group is generated by elements a_1, ..., a_g, b_1, ..., b_g, g_1, ..., g_{t+s} satisfying the relations g_i^m_i=1 for 1<=i<=t and [a_1,b_1]*...*[a_g,b_g]*g_1*...*g_{t+s}=1, where [x, y]=x*y*y^(-1)*x^(-1) is the commutator.");
	install("algfromnormdisc","GGG","algfromnormdisc","./libfdom.so");
	addhelp(algfromnormdisc,"Inputs F, D, infram: number field F, positive real D, vector infram of 0/1's of length deg(F).\n Returns a quaternion algebra over F with infinite ramification infram and discriminant disc, where |N_{F/Q}(disc)|=D, if it exists. If it does not exist, returns 0.");	
	install("algmulvec","GGG","algmulvec","./libfdom.so");
	addhelp(algmulvec,"Inputs A, G, L: algebra A, G=vector of elements of A, L a vecsmall or vector of indices.\n Returns G[L[1]]*G[L[2]]*...*G[L[n]].");
	install("algramifiedplacesf","G","algramifiedplacesf","./libfdom.so");
	addhelp(algramifiedplacesf,"Input A, an algebra.\n Returns the vector of finite places that ramify.");
	install("algnormalizedbasis","GGGp","algnormalizedbasis","./libfdom.so");
	addhelp(algnormalizedbasis, "Inputs A, G, p: quaternion algebra A split at one real place, set G of elements of norm 1 in the order in A, upper half plane point p.\n Returns the normalized basis associated to G.");
	install("algnormalizedboundary","GGGp","algnormalizedboundary","./libfdom.so");
	addhelp(algnormalizedboundary, "Inputs A, G, p: quaternion algebra A split at one real place, set G of elements of norm 1 in the order in A, upper half plane point p.\n Returns the normalized boundary associated to G. The format of the output is [elements, icircs, vertices, vertex angles, matrices, area, 0, mats]. The circle corresponding to elements[i] is icircs[i], and the vertices are vertices[i-1] and vertices[i]. matrices[i] is the image in PSU(1,1) of elements[i]. The element 1 corresponds to a section on the unit circle, which also corresponds to a circle of 0. Vertex angles stores the radial angle to the ith vertex (with base angle being the first one). The area is the area, and the 0 stores the side pairing when we have a fundamental domain (so a priori stores nothing).");
	install("algshimura","GGD1,L,","algshimura","./libfdom.so");
	addhelp(algshimura,"Inputs F, D, {place=1}: totally real number field F, positive integer D, integer place between 1 and deg(F).\n Returns a quaternion algebra over F that is split at the infinite place place only, and has discriminant disc, where |N_{F/Q}(disc)|=D, if it exists. If it does not exist, returns 0. This also guarantees that a>0 at the split infinite place, hence the output is suitable for fundamental domain methods.");
	install("algsmallnorm1elts","GGGD0,G,D0,G,p","algsmallnorm1elts","./libfdom.so");
	addhelp(algsmallnorm1elts,"Inputs A, C, p, {z1=0}, {z2=0}: quaternion algebra A split at one real place, positive real number C, upper half plane point p, unit disc point z.\n Computes small norm 1 elements in the order of A, i.e. such that Q_{z_1,z_2}(g)<=C, where Q is defined on page 478 of Voight ``Computing fundamental domains''. The point p is the base for the mapping from the upper half plane model to the unit disc model, and z1, z2 are basepoints (the inverse radius is computed for h_2^(-1)*g*h_1, where h_i(0)=z_i; hence elements g with g(z_1) close to z_2 are found).");
	install("algswapab","G","algswapab","./libfdom.so");
	addhelp(algswapab,"Input A, a quaternion algebra=(a, b/F).\n Returns (b, a/F), i.e. swapping a and b.");
	
\\TEMPORARY
install("completebasisdet1","G","completebasisdet1","./libfdom.so");
install("algnormform","Gp","algnormform","./libfdom.so");
install("algfdom_test","GGD1,L,D0,G,D0,G,p","algfdom1","./libfdom.so");
install("algabsrednorm","GGD0,G,D0,G,p","algabsrednorm","./libfdom.so");
install("algenum","GGGGGLp","algenum","./libfdom.so");
install("balltester","GGGp","balltester","./libfdom.so");
install("bestAval","GGp","bestAval","./libfdom.so");
install("Ntries","GGGGGLLp","Ntries","./libfdom.so");

install("algsmallnorm1elts_condition","GGGD0,G,D0,G,D0,L,D0,L,p","algsmallnorm1elts1","./libfdom.so");
install("algsmallnorm1elts_condition2","GGGD0,G,D0,G,D0,L,D0,L,D0,L,p","algsmallnorm1elts2","./libfdom.so");


