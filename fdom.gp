default(help, "gphelp -detex");

install("algfdarea","Gp","algfdarea","./libfdom.so");
install("algfd","GGD1,L,D0,G,D0,G,p","algfd","./libfdom.so");
install("algramifiedplacesf","G","algramifiedplacesf","./libfdom.so");
install("algsmallnorm1elts","GGGD0,G,p","algsmallnorm1elts","./libfdom.so");


install("python_printarcs","vGrD0,L,Drp","python_printarcs","./libfdom.so");
addhelp(python_printarcs,"Inputs arcs, filename, {view=0}, {extrainput=NULL}: a set of arcs arcs, string filename, view=0, 1, and extrainput=NULL or a string.\n Prints the arcs specified by arcs to the file fdoms/filename.dat, ready for plotviewer. If view=1, calls plotviewer with the additional input of extrainput if you want to include other arcs/fundamental domains.");
install("python_plotviewer","vr","python_plotviewer","./libfdom.so");
addhelp(python_plotviewer,"Input S: string denoting the file names of data fundamental domains/geodesics.\n Launches the python file fdviewer.py to view the domain/geodesics. Enter the files separated by spaces (they must be stored in the sub-folder 'fdoms').");
install("python_printfdom","vGrp","python_printfdom","./libfdom.so");
addhelp(python_printfdom,"Input U, filename: fundamental domain U, string filename.\n Prints the fundamental domain U to the file fdoms/filename.dat, ready for the plot viewer. The filename must start with 'fd' to work properly.")