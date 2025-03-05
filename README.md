# Fundamental domains for congruence arithmetic Fuchsian groups.
Installation: See below.

Using the project (and quaternion algebras in general): [User's Manual](Documentation/QuaternionAlgebras_PARIGP.pdf).

For a sample use of the main methods, see the file ```example.gp```.

## References
The code is originally based on the paper [Improved computation of fundamental domains for arithmetic Fuchsian groups](https://doi.org/10.1090/mcom/3777) ([Arxiv](https://arxiv.org/abs/2110.11503)), which builds off of the papers [Computing fundamental domains
for Fuchsian groups](https://math.dartmouth.edu/~jvoight/articles/funddom-jtnb-fixederrata.pdf) and [Computing arithmetic Kleinian groups](http://www.normalesup.org/~page/Recherche/Documents/articles/kln_gps.pdf). If you make use of the code, please cite this paper, as well as the GitHub repository. A suggested Bibtex entry for this repository is
```
@misc{FundamentalDomains,
  AUTHOR = {Rickards, James},
  TITLE = {Fundamental domains for {S}himura curves},
  YEAR = {2025},
  PUBLISHER = {GitHub},
  JOURNAL = {GitHub repository},
  HOWPUBLISHED = {\url{https://github.com/JamesRickards-Canada/Fundamental-domains-for-Shimura-curves}},
}
```

## Installation Instructions

### Downloading the code
Call ```git clone https://github.com/JamesRickards-Canada/Fundamental-domains-for-Shimura-curves.git```. If you are on Windows, be sure to git clone from WSL (see below), as Windows line endings (carriage returns) may be added to files, causing issues.

Ensure that you are on the branch ```klein```. The code in the ```paper``` branch matches the code when the paper was written, and data in the paper can be recreated using this branch. This branch is no longer updated. The default branch, ```klein```, is a significantly improved version of this code, and should be the branch of choice for most of users.

### Prerequisites
* PARI/GP, but _not_ the downloaded ready-to-go binary. The PARI/GP website has binaries for Windows and Mac avaliable, but these will not work with the package. See below for OS specific instructions.
* Version at least 2.15.1, though the more up-to-date the better. This may work on 2.13.x, but this is not suggested, as there are various quaternion algebra bugs that are unfixed.
* You should have a guess as to the location of the ```pari.cfg``` file for the version of PARI/GP you are running. Suggestions on how to do this can be found below.

### Operating systems
* **Linux** - No further requirements
* **Windows** - You need to use Windows Subsytem for Linux. See the [guide](https://pari.math.u-bordeaux.fr/PDF/PARIwithWindows.pdf) I wrote for additional instructions.
* **Mac** - You need to have [Homebrew](https://brew.sh/) installed. This is also an easy way to install PARI/GP: ```brew install pari```

### Where is pari.cfg?
* The configuration file will search for this, but it is preferrable to not search your entire hard drive (as this can be very slow). So, you should at least supply a guess as to the location of ```pari.cfg```. Often only the top-level folder (e.g. ```/usr``` or ```/opt```) suffices.
* On Linux or WSL, if you build PARI/GP from source, it should be located in ```/usr/local/lib/pari/pari.cfg```, or at least somewhere in the ```/usr``` folder.
* On a Mac, if you install PARI/GP with Homebrew, it may be found in a folder like ```/opt/homebrew/Cellar/pari/VERSION/lib/pari```. Searching ```/opt/homebrew``` should be fine.
* If you are obtaining it through SageMath, it might be found where the library files of SageMath are
* Assuming you open PARI/GP with the command ```gp```, try ```type -a gp```, which will display where this command lives. The corresponding file(s) are likely symbolic links, and you can call ```readlink -f LOCATION``` on each of them to see where it lives. This can provide a clue as to the place to search for ```pari.cfg```.
* Another clue comes from gp itself. Open gp, and type ```default()```. Look for the entries ```datadir``` and ```help```. It is _often_ the case that ```datadir``` is in ```X/share/pari```, ```help``` is in ```X/bin/gphelp```, and ```pari.cfg``` is in ```X/lib/pari/pari.cfg```.

### Configuring and building the package
* From inside the project folder, call ```./configure``` to initialize the project. This helps you search for ```pari.cfg```, and stores the location to a file. You should supply it with a folder to search in!
* The script displays the corresponding versions of the found files, so if you have multiple versions, you can choose the correct one. This can be useful if you keep multiple copies of PARI/GP around.
* If the location of the installation of PARI/GP does not change, you do not need to reconfigure. If when you update PARI/GP there is a new location (e.g. if the version number is in the file path of ```pari.cfg```), you should call ```./configure``` again.
* Call ```make``` to build the project, and ```make clean``` to remove all .o object files. If you update to a new version of PARI/GP, you must remake the project.
* Once this is done, a call to ```gp fdom``` starts gp with the package installed!
* Call ```?fdom``` or consult the [User's Manual](Documentation/QuaternionAlgebras_PARIGP.pdf) to access further help for this package.

### Optional packages
* LaTeX compiler ```pdftex```, which has the ```standalone``` and ```tikz``` packages. This allows automatic compilation of fundamental domains in LaTeX.
* Python with ```matplotlib``` and ```numpy```, for an interactive viewer of fundamental domains and geodesics.

### Troubleshooting
* If you are encoutering unexpected errors or warnings when building, ensure you found the correct ```pari.cfg``` file. You need to configure the project with the same version of PARI/GP you are using.
* If you upgrade your version of PARI/GP, you should ```make clean```, call ```configure``` again, and remake the project.
* If you still have trouble installing or using the package, please get in touch!

## How to use the methods

Full instructions can be found in the [User's Manual](Documentation/QuaternionAlgebras_PARIGP.pdf).
