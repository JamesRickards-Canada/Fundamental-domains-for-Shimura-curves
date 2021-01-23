print("\n\nType '?fdom' for help.\n\n");
addhelp(fdom, "This package can be used to compute fundamental domains for Shimura curves (the PARI code can be easily adapted to compute fundamental domains for any discrete subgroup of PSL(2, R)).\n For each subtopic ``P (p)'', call ?p to access a basic description and list of methods. Subtopics:\n Visualizing fundamental domains with Python (vfd)\n Quaternion methods (quat)");
default(help, "gphelp -detex");

\\Visualization
	addhelp(vfd,"These methods allow one to save fundamental domains and geodesics, and view them with a Python program. Available methods: python_plotviewer, python_printarcs, python_printfdom.");
	
	install("python_plotviewer","vr","python_plotviewer","./libfdom.so");
	addhelp(python_plotviewer,"Input S: string denoting the file names of data fundamental domains/geodesics.\n Launches the python file fdviewer.py to view the domain/geodesics. Enter the files separated by spaces (they must be stored in the sub-folder 'fdoms').");
	install("python_printarcs","vGrD0,L,Drp","python_printarcs","./libfdom.so");
	addhelp(python_printarcs,"Inputs arcs, filename, {view=0}, {extrainput=NULL}: a set of arcs arcs, string filename, view=0, 1, and extrainput=NULL or a string.\n Prints the arcs specified by arcs to the file fdoms/filename.dat, ready for plotviewer. If view=1, calls plotviewer with the additional input of extrainput if you want to include other arcs/fundamental domains.");
	install("python_printfdom","vGrp","python_printfdom","./libfdom.so");
	addhelp(python_printfdom,"Input U, filename: fundamental domain U, string filename.\n Prints the fundamental domain U to the file fdoms/filename.dat, ready for the plot viewer. The filename must start with 'fd' to work properly.")


\\Quaternion methods
	addhelp(quat,"These methods allow for the computation of fundamental domains for Eichler orders in quaternion algebras split at one real place. Available methods: algfdomarea, algfdom, algramifiedplacesf, algsmallnorm1elts.");

	install("algfdomarea","Gp","algfdomarea","./libfdom.so");
	addhelp(algfdomarea,"Input A, a quaternion algebra split at one real place.\n Returns the area of the fundamental domain associated to the group of units of norm 1 in the order stored in A. Requires the computation of the zeta_K(2) (Dedekind zeta function for the centre), which may require some calls to allocatemem() if K is ``large''.");
	install("algfdom","GGD1,L,D0,G,D0,G,p","algfdom","./libfdom.so");
	addhelp(algfdom,"Inputs A, p, {dispprogress=1}, {area=0}, {ANRdata=0}: quaternion algebra A split at one real place, upper half plane point p, dispprogress=0,1, ANRdata=0 or a length 5 vector, area=0 or the area of the fundamental domain.\n Computes and returns the fundamental domain for the group of units of norm 1 in A. We use the unit disc model, which the upper half plane is mapped to via p->0. p need to NOT be a fixed point of this group under the standard embedding into PSL(2, R) (p=I/2 is safe for quaternion algebras over Q). If area is passed in, this method will not re-compute it, which can save a little bit of time. ANRdata is a technical entry, and can be passed in to specify some/all constants used in the enumeration of Page (they greatly affect the running time, but the optimal choices are not totally clear).");
	install("algramifiedplacesf","G","algramifiedplacesf","./libfdom.so");
	addhelp(algramifiedplacesf,"Input A, an algebra.\n Returns the vector of finite places that ramify.");
	install("algsmallnorm1elts","GGGD0,G,p","algsmallnorm1elts","./libfdom.so");
	addhelp(algsmallnorm1elts,"Inputs A, C, p, {z=0}: quaternion algebra A split at one real place, positive real number C, upper half plane point p, unit disc point z.\n Computes small norm 1 elements in the order of A, i.e. such that absrednorm(g)<=C, where absrednorm is defined on page 478 of Voight ``Computing fundamental domains''. The point p is the base for the mapping from the upper half plane model to the unit disc model, and z is the basepoint in the unit disc model for absrednorm (the invrad part is computed with respect to z, hence elements with isometric circles close to z are found).");