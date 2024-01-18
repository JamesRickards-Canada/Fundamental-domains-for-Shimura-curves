/*Runs the fundamental domain methods across a method suite to compare to Magma. We don't use orders, and compute for O^1 here.
Be aware that if gp is built with pthread, then these timings are not quite accurate! This gives the CPU time, which will not correlate with the timings on a single threaded gp installation.*/
avg_stdev(v)={
  my(avg, n, va);
  n=#v;
  avg=vecsum(v)/n;
  va=sum(i=1,n,(v[i]-avg)^2);
  return([avg, sqrt(va/n)]);
}


avgtime(fil, A, Or, case, N)={
  my(times, X, dat, dom, st, extra);
  times=vector(N);
  for(i=1,N,
    gettime();
    X=afuchinit(A, , 0, 1);
    times[i]=gettime()/1000;
  );
  dat=avg_stdev(times);
  st=strprintf("%4d   %2d   %P9.3f %8d %P12.3f %P14.3f", case, algcenter(A).r1, afucharea(X), #afuchelts(X), dat[1], dat[2]);
  filewrite(fil, st);
  print(st);
}

default(parisize, "4096M");/*4 GB of memory.*/

/*Details*/
testspercase=10;
Ntests=18;
tno=1;

/*Starting the file*/
#/*My timer is default on*/
fil=fileopen("guide_timings_gp.txt", "a");
rgen=getwalltime();
setrand(rgen);
st=strprintf("Testing the running time of algfdom, in seconds. There are %d cases, and each case is tested %d times. The initial random seed is %d.\nCase |  n |      Area | #Sides | Time (avg) | Time (stdev)", Ntests, testspercase, rgen);
filewrite(fil, st);
print(st);

/*The data. We put the tests on new lines (instead of a for loop) to reset the pari stack each time.*/
F=nfinit(y);

A=alginit(F, [-733, 775515]);/*Ramified at 2, 3, 5, 13*/
avgtime(fil, A, 0, tno, testspercase);tno++;/*Area 100.53096491487338363*/

A=alginit(F, [1022, 9680385]);/*Ramified at 3, 5, 11, 13*/
avgtime(fil, A, 0, tno, testspercase);tno++;/*Area 1005.3096491487338363*/

A=alginit(F, [-1, 2021]);/*Ramified at 43, 47*/
avgtime(fil, A, 0, tno, testspercase);tno++;/*Area 2023.1856689118268456*/

A=alginit(F, [2070, -45999]);/*Ramified at 19, 269*/
avgtime(fil, A, 0, tno, testspercase);tno++;/*Area 5051.680986972387528*/

F=nfinit(y^2-y-1);
A=alginit(F, [-16*y + 3, -6734*y + 5873]);/*Ram at primes over 2, 7, 19*/
avgtime(fil, A, 0, tno, testspercase);tno++;/*Area 542.8672105403162716*/

F=nfinit(y^2-11);
A=alginit(F, [-1, -7*y - 16]);/*Ram at prime over 283*/
avgtime(fil, A, 0, tno, testspercase);tno++;/*Area 2067.167966062083951*/

F=nfinit(y^3-3*y-1);
A=alginit(F, [-1, -4*y^2 - 80*y - 75]);/*Ram over 3, 503*/
avgtime(fil, A, 0, tno, testspercase);tno++;/*Area 350.4621138004613790*/

F=nfinit(y^3-y^2-4*y+3);
A=alginit(F, [580*y^2 + 656*y - 2851, -16*y^2 - 16*y + 7]);/*Ramified above 3, 47*/
avgtime(fil, A, 0, tno, testspercase);tno++;/*Area 770.7373976806959411*/

F=nfinit(y^4+2*y^3-4*y^2-6*y+1);
A=alginit(F, [-92*y^3 - 172*y^2 + 352*y + 221, 2*y^3 + 4*y^2 - 7*y - 18]);/*Ramified above 17*/
avgtime(fil, A, 0, tno, testspercase);tno++;/*Area 469.1445029360757902*/

F=nfinit(y^4 - 4*y^2 - y + 1);
A=alginit(F, [-1, -48*y^3 + 36*y^2 + 120*y - 111]);/*Ram above 3, 3, 31*/
avgtime(fil, A, 0, tno, testspercase);tno++;/*Area 1633.6281798666922729*/

\p57
F=nfinit(y^5-5*y^3-y^2+3*y+1);
A=alginit(F, [-104*y^4 + 116*y^3 + 672*y^2 - 660*y - 867, y^2 - y - 10]);/*Ramified above 5 and 83*/
avgtime(fil, A, 0, tno, testspercase);tno++;/*area 343.48079679248406073*/

F=nfinit(y^5 - 2*y^4 - 5*y^3 + 9*y^2 + 5*y - 7);
A=alginit(F, [8*y^4 - 104*y^3 + 164*y^2 + 188*y - 615, y - 11]);/*Ramified above 5 and 7*/
avgtime(fil, A, 0, tno, testspercase);tno++;/*area 829.38046054770541495*/

\p38
F=nfinit(y^6-y^5-5*y^4+4*y^3+6*y^2-3*y-1);
A=alginit(F, [4*y^4 - 12*y^3 + 12*y - 47, -1]);/*Ramified above 131*/
avgtime(fil, A, 0, tno, testspercase);tno++;/*area 198.96753472735357176*/

F=nfinit(y^6-y^5-7*y^4+2*y^3+7*y^2-2*y-1);
A=alginit(F, [-1, -8*y^5 + 4*y^4 + 48*y^3 + 32*y^2 - 24*y - 59]);/*Ramified above 491*/
avgtime(fil, A, 0, tno, testspercase);tno++;/*area 542.4483315198376323*/

\p57
F=nfinit(y^7 - 8*y^5 + 19*y^3 - y^2 - 13*y + 1);
A=alginit(F, [-4*y^5 + y^4 + 20*y^3 - 2*y^2 - 16*y - 3, -10*y^6 - 2*y^5 + 58*y^4 + 21*y^3 - 74*y^2 - 20*y - 15]);/*Unramified*/
avgtime(fil, A, 0, tno, testspercase);tno++;/*area 138.2300767579414814*/

F=nfinit(y^7 - 2*y^6 - 6*y^5 + 8*y^4 + 12*y^3 - 5*y^2 - 6*y - 1);
A=alginit(F, [-4*y^6 + 18*y^5 - 11*y^4 - 37*y^3 + 25*y^2 + 22*y - 1, -21*y^6 + 17*y^5 + 176*y^4 - 20*y^3 - 400*y^2 - 183*y - 33]);/*Unramified*/
avgtime(fil, A, 0, tno, testspercase);tno++;/*area 238.7610416728242861*/

F=nfinit(y^8 - 4*y^7 - y^6 + 17*y^5 - 5*y^4 - 23*y^3 + 6*y^2 + 9*y - 1);
A=alginit(F, [-1, -16*y^7 + 68*y^6 - 20*y^5 - 204*y^4 + 168*y^3 + 120*y^2 - 140*y - 3]);/*Ramified above 19*/
avgtime(fil, A, 0, tno, testspercase);tno++;/*area 422.23005264246821124937542704619344205*/

F=nfinit(y^8 - 2*y^7 - 7*y^6 + 11*y^5 + 14*y^4 - 18*y^3 - 8*y^2 + 9*y - 1);
A=alginit(F, [-1, -12*y^7 + 16*y^6 + 76*y^5 - 48*y^4 - 96*y^3 + 28*y^2 + 16*y - 35]);/*Ramified above 11*/
avgtime(fil, A, 0, tno, testspercase);tno++;/*area 423.06781068342548944634012072981623082*/


/*End of tests*/
filewrite(fil, "End of tests\n\n");
fileclose(fil);
print("End of tests\n\n");
#/*Turn the timer back on*/