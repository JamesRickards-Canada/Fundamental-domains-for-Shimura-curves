# Fundamental domains for congruence arithmetic Fuchsian groups.

## References
The code is originally based on the paper [Improved computation of fundamental domains for arithmetic Fuchsian groups](https://doi.org/10.1090/mcom/3777) ([Arxiv](https://arxiv.org/abs/2110.11503)), which builds off of the papers [Computing fundamental domains
for Fuchsian groups](https://math.dartmouth.edu/~jvoight/articles/funddom-jtnb-fixederrata.pdf) and [Computing arithmetic Kleinian groups](http://www.normalesup.org/~page/Recherche/Documents/articles/kln_gps.pdf). If you make use of the code, please cite this paper, as well as the GitHub repository.

## Installation Instructions
The code in the "paper" branch matches the code when the paper was written, and data in the paper can be recreated using this branch. This branch is no longer updated. The default branch, klein, is a significantly improved version of this code, and should be the branch of choice for most of users.

### Prerequisites
* PARI/GP, but _not_ the downloaded ready-to-go binary. The PARI/GP website has binaries for Windows and Mac avaliable, but these will not work with the package. See below for OS specific instructions.
* Version at least 2.15, though the more up-to-date the better.
* You should have a guess as to the location of the ```pari.cfg``` file for the version of PARI/GP you are running. Suggestions on how to do this can be found below.

### Operating systems
* **Linux** - No further requirements
* **Windows** - You need to use Windows Subsytem for Linux. Further instructions can be found [here](https://pari.math.u-bordeaux.fr/PDF/PARIwithWindows.pdf).
* **Mac** - You need to have [Homebrew](https://brew.sh/) installed. This is also an easy way to install PARI/GP: ```brew install pari```

### Where is pari.cfg?
* On Linux or WSL, if you build PARI/GP from source, it should be located in ```/usr/local/lib/pari/pari.cfg```, or at least somewhere in the ```/usr``` folder.
* On a Mac, if you install with Homebrew, it may be found in a folder like ```/opt/homebrew/Cellar/pari/VERSION/lib/pari```
* If you are obtaining it through SageMath, it might be found where the library files of SageMath are
* Assuming you open PARI/GP with the command ```gp```, try ```type -a gp```, which will display where this command lives. The corresponding file(s) are likely symbolic links, and you can call ```readlink -f LOCATION``` on each of them to see where it lives.
* In absolute doubt, the configuration method allows you to search the entire system for the file. This should only be done as a last resort, as the search could be quite slow!

### Configuring and building the package
* Call ```./configure``` to initialize the project. This helps you search for ```pari.cfg```, and stores the location to a file. It displays the corresponding versions of the found files, so if you have multiple versions, you can choose the correct one.
* As long as the location of the installation of PARI/GP does not change, you do not need to reconfigure ever.
* Call ```make``` to build the project, and ```make clean``` to remove all .o object files. If you update to a new version of PARI/GP, you must remake the project.
* Once this is done, a call to ```gp fdom``` starts gp with the package installed!

### Optional packages
* LaTeX compiler ```pdftex```, which has the ```standalone``` and ```tikz``` packages. This allows automatic compilation of fundamental domains in LaTeX.
* Python with ```matplotlib```, for an interactive viewer of fundamental domains and geodesics.

## How to use the methods

Full instructions can be found in the [User's Manual](Documentation/QuaternionAlgebras_PARIGP.pdf).
