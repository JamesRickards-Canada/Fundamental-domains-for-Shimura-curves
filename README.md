# Fundamental domains for Shimura curves

This package is a PARI/GP implementation of fundamental domains for Shimura curves. You can compute fundamental domains and view them with a Python application (using Matplotlib).

This package is written for the development version of PARI/GP, 2.14. If you are using a stable release, you will have to call "make" in order to recompile the library (.so) file for your own system. In the Makefile, you will have to update the location of the PARI/GP installation to where it is on your own system. Alternatively, there is a .so file for the 2.13 version stored in the folder "2.13_sofile". This may not be updated as often as the code is, but will generally be fine.
Once you have the .so file compiled, you can call "gp fdom" to open gp and load the package.

For a longer introduction on quaternion algebras in PARI/GP & this package, see the file "QuaternionAlgebras_PARIGP.pdf".

Type ?fdom to access the help inside of gp.

Run the file "example.gp" to see an example of these methods in action.

The files paper_computations.gp, papertimetest.gp, and the folder data_in can be used to reproduce figures from my paper, "Improved computation of fundamental domains for arithmetic Fuchsian groups". They are not necessary for the operation of the package.
