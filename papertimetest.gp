\\This file contains methods used to generate data in the paper "IMPROVED COMPUTATION OF FUNDAMENTAL DOMAINS FOR ARITHMETIC FUCHSIAN GROUPS".

avgtime(fil, alg, case, N)={
  my(dom, T, st, extra);
  gettime();
  for(i=1,N,dom=algfdom(alg));
  T=gettime()/(1000*N);
  st=strprintf("Test case %d: %.3f s", case, T);
  filewrite(fil, st);
  print(st);
  extra=length(digits(case));
  st=strprintf("      Area: %.10f", dom[6]);
  for(i=1,extra,st=strprintf(" %s", st));\\Adjusting for the length of the case number
  filewrite(fil, st);
  print(st);
  st=strprintf("     Sides: %d", length(dom[1]));
  for(i=1,extra,st=strprintf(" %s", st));\\Adjusting for the length of the case number
  filewrite(fil, st);
  print(st);
}

default(parisize, "4096M");\\4 GB of memory.

\\Details
st1="Test cases for paper";
testspercase=5;
Ntests=8;

\\Initialization
F=vector(Ntests);
A=vector(Ntests);

\\The data
F[1]=nfinit(y);
A[1]=alginit(F[1], [33, -1]);\\area 20.943951023931954923084289221863352561;

F[2]=nfinit(y);\\Ramified at 13, 61
A[2]=alginit(F[2], [793, -7]);\\area 753.98223686155037723103441198708069220;

F[3]=nfinit(y^2-11);
A[3]=alginit(F[3], [7*y-12, -5]);\\area 571.76986295334236940020109575686952492;

F[4]=nfinit(y^2-33);
A[4]=alginit(F[4], [-6*y-26, -8*y-15]);\\area 226.19467105846511316931032359612420766

F[5]=nfinit(y^3-5*y+1);
A[5]=alginit(F[5], [5*y-7,-300]);\\area 418.87902047863909846168578443726705123;

F[6]=nfinit(y^4+2*y^3-4*y^2-6*y+1);
A[6]=alginit(F[6], [36*y^2-4*y-223, y^2-y-10]);\\area 469.14450293607579027708807856973909737;

\p50 \\We need more precision for this case
F[7]=nfinit(y^5-9*y^4-7*y^3+15*y^2+14*y+3);
A[7]=alginit(F[7], [-40*y^4 + 385*y^3 + 35*y^2 - 575*y - 240, 25*y^4 - 235*y^3 - 85*y^2 + 415*y + 151]);\\area 4490.3830995310111355092716091675027888;

\p80 \\We need even more precision for this case
F[8]=nfinit(y^7-y^6-6*y^5+4*y^4+10*y^3-4*y^2-4*y+1);
A[8]=alginit(F[8], [y^6-5*y^5+3*y^4+13*y^3-10*y^2-8*y-8, -264*y^6+136*y^5+1632*y^4-432*y^3-2444*y^2+268*y-11]);\\area 1507.9644737231007544620688239741613602;



\\Starting the file
fil=fileopen("papertimetest_results.txt", "a");
st2=strprintf("Testing the running time of algfdom. There are %d cases, and each case is tested %d times.\n The algorithm details are: %s", Ntests, testspercase, st1);
filewrite(fil, st2);
print(st2);

\\We put the tests on new lines (instead of a for loop) to reset the pari stack each time.
i=1;
\p38
avgtime(fil, A[i], i, testspercase);i++;\\1
avgtime(fil, A[i], i, testspercase);i++;\\2
avgtime(fil, A[i], i, testspercase);i++;\\3
avgtime(fil, A[i], i, testspercase);i++;\\4
avgtime(fil, A[i], i, testspercase);i++;\\5
avgtime(fil, A[i], i, testspercase);i++;\\6
\p50
avgtime(fil, A[i], i, testspercase);i++;\\7
\p80
avgtime(fil, A[i], i, testspercase);i++;\\8

\\End of tests
filewrite(fil, "End of tests\n\n");
fileclose(fil);
print("End of tests\n\n");