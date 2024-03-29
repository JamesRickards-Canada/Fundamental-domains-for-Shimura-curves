/*Runs the fundamental domain methods across a testing suite to give a series of benchmarks.
Be aware that if gp is built with pthread, then these timings are not quite accurate! This gives the CPU time, which will not correlate with the timings on a single threaded gp installation.*/

avg_stdev(v) = {
  my(avg, n, va);
  n = #v;
  avg = vecsum(v)/n;
  va = sum(i = 1, n, (v[i] - avg)^2);
  return([avg, sqrt(va/n)]);
}


avgtime(fil, A, Or, case, N) = {
  my(times, X, dat, dom, st, extra);
  times = vector(N);
  for(i = 1, N,
    gettime();
    if(Or, X = afuchinit(A, Or, , 1), X = afuchinit(A, , , 1));
    times[i] = gettime()/1000;
  );
  dat = avg_stdev(times);
  st = strprintf("%4d   %2d   %P9.3f %8d %P12.3f %P14.3f", case, algcenter(A).r1, afucharea(X), #afuchelts(X), dat[1], dat[2]);
  filewrite(fil, st);
  print(st);
}

default(parisize, "4096M");/*4 GB of memory.*/

/*Details*/
testspercase = 3;
Ntests = 34;
tno = 1;

/*Starting the file*/
#/*My timer is default on*/
fil = fileopen("benchmarks.txt", "a");
st=strprintf("Testing the running time of algfdom, in seconds. There are %d cases, and each case is tested %d times.\nCase |  n |      Area | #Sides | Time (avg) | Time (stdev)", Ntests, testspercase);
filewrite(fil, st);
print(st);

/*The data. We put the tests on new lines (instead of a for loop) to reset the pari stack each time.*/
setrand(1);
F = nfinit(y);

setrand(1);
A = alginit(F, [793, -7]);/*Ramified at 13, 61*/
avgtime(fil, A, 0, tno, testspercase); tno++;/*area 753.98223686155037723103441198708069220;*/

setrand(1);
A = alginit(F, [14820, -343]);/*Ramified at 3, 5, 13, 19*/
avgtime(fil, A, 0, tno, testspercase); tno++;/*area 1809.5573684677209053544825887689936613*/

setrand(1);
A = alginit(F, [67, 61]);/*Ramified at 61, 67*/
avgtime(fil, A, 0, tno, testspercase); tno++;/*area 4146.9023027385270747706892659289438071*/

setrand(1);
A = alginit(F, [7975922955, -39733]);/*Ramified at 2, 3, 5, 7, 11, 13*/
avgtime(fil, A, 0, tno, testspercase); tno++;/*area 6031.8578948924030178482752958966455376*/

setrand(1);
A = alginit(F, [103, 101]);/*Ramified at 101, 103*/
avgtime(fil, A, 0, tno, testspercase); tno++;/*area 10681.415022205297010772987503150309806*/

setrand(1);
A = alginit(F, [3, 1]);/*Unramified*/
Or = algorderalgtoorder(A, [[1, 0]~, [900*x, 0]~, [468*x + 1/2, 1/2]~, [2375/6*x + 1/2, 1/6*x + 1/2]~]);/*Level 900*/
avgtime(fil, A, Or, tno, testspercase); tno++;/*area 2261.9467105846511316931032359612420766*/

setrand(1);
A = alginit(F, [11, -3]);/*Ramified at 3, 11*/
Or = algorderalgtoorder(A, [[1, 0]~, [115*x, 0]~, [197/2*x, 1/2]~, [63/2*x + 1/2, 1/2*x + 1/2]~]);/*Level 115*/
avgtime(fil, A, Or, tno, testspercase); tno++;/*area 3015.9289474462015089241376479483227688*/

setrand(1);
A = alginit(F, [31, 29]);/*Ramified at 29, 31*/
Or = algorderalgtoorder(A, [[1, 0]~, [x, 0]~, [5/2*x, 5/2]~, [2*x + 1/2, 1/2*x + 2]~]);/*Level 5*/
avgtime(fil, A, Or, tno, testspercase); tno++;/*area 5277.8756580308526406172408839095648454*/


setrand(1);
F = nfinit(y^2 - 2);
setrand(1);
A = alginit(F, [y, -17]);/*Ramified at primes above 2, 17, 17*/
Or = algorderalgtoorder(A, [[1, 0]~, [5*x, 0]~, [y, 0]~, [5*y*x, 0]~, [(3/2*y + 1)*x + (1/2*y + 1/2), 1/2]~, [(7/2*y + 3/2)*x, 1/2*x]~, [(y + 3)*x + 1/2*y, 1/2*y]~, [(3/2*y + 2)*x, 1/2*y*x]~]);/*Level 5*/
avgtime(fil, A, Or, tno, testspercase); tno++;/*area 3485.0734503822772992012257265180618662*/

setrand(1);
F = nfinit(y^2 - 11);
setrand(1);
A = alginit(F, [7*y - 12, -5]);/*ramified above 79*/
avgtime(fil, A, 0, tno, testspercase); tno++;/*area 571.76986295334236940020109575686952492;*/

setrand(1);
F = nfinit(y^2 - 33);
setrand(1);
A = alginit(F, [-6*y - 26, -8*y - 15]);/*Ramified at a prime above 37*/
Or = algorderalgtoorder(A, [[1, 0]~, [(-17/16*y + 119/16)*x, 0]~, [(-7/16*y + 57/16)*x, 0]~, [-1/2*y - 1/2, 0]~, [(-1/4*y + 7/4)*x + 1/2, 1/2]~, [(-29/32*y + 203/32)*x + 1/2, (-1/32*y + 7/32)*x + 1/2]~, [(-43/64*y + 309/64)*x, (1/64*y + 1/64)*x]~, [(-13/32*y + 107/32)*x + (-1/4*y + 3/4), (-23/1632*y + 59/544)*x + (-1/204*y + 41/68)]~]);/*Level a prime above 17*/
avgtime(fil, A, Or, tno, testspercase); tno++;/*area 4071.5040790523720370475858247302357379*/

setrand(1);
F = nfinit(y^2 - 122);
setrand(1);
A = alginit(F, [-y - 1, -3]);/*Ramified at 3.*/
avgtime(fil, A, 0, tno, testspercase); tno++;/*area 1935.2210746113126348929883241001737767*/

setrand(1);
F = nfinit(y^3 - 3*y - 1);
setrand(1);
A = alginit(F, [192*y^2 - 272*y - 487, -3*y^2 + 5*y - 6]);/*Ramified at primes above 3, 37*/
Or = algorderalgtoorder(A, [[1, 0]~, [-y^2 + y + 2, 0]~, [-y, 0]~, [17/2*x - 17/2, 0]~, [(-1/2*y + 2)*x + (-1/2*y - 2), 0]~, [(-1/2*y^2 + 8)*x + (1/2*y^2 - 8), 0]~, [6*x - 6, 1]~, [8*x - 8, -y^2 + y + 2]~, [3/2*x - 3/2, -y]~, [5/2*x - 5/2, 1/2*x - 1/2]~, [(-13/34*y^2 - 5/34*y + 80/17)*x + (-1/2*y^2 + 1/2*y - 2), (-1/34*y + 5/17)*x + (-9/17*y^2 + 13/34*y + 19/17)]~, [(-438199/810254*y^2 - 337709/810254*y + 5546293/810254)*x + (-1/2*y^2 + 1/2*y - 9/2), (-1/57528034*y^2 - 361491/57528034*y + 13975085/28764017)*x + (-141/142*y^2 + 101/142*y + 129/71)]~]);/*Level a prime above 17.*/
avgtime(fil, A, Or, tno, testspercase); tno++;/*area 452.38934211693022633862064719224841532*/

setrand(1);
F = nfinit(y^3 - y^2 - 4*y + 3);
setrand(1);
A = alginit(F, [580*y^2 + 656*y - 2851, -16*y^2 - 16*y + 7]);/*Ramified above 3, 47*/
avgtime(fil, A, 0, tno, testspercase); tno++;/*area 770.73739768069594116950184336457137425*/

setrand(1);
F = nfinit(y^3 - 5*y + 1);
setrand(1);
A = alginit(F, [5*y - 7, -300]);/*Ramified at primes above 3, 11*/
Or = algorderalgtoorder(A, [[1, 0]~, [y, 0]~, [y^2 - 3, 0]~, [5*x, 0]~, [(y + 1)*x, 0]~, [(y^2 + y)*x, 0]~, [(1/2*y^2 + y + 2)*x, 1/20]~, [(1/2*y + 3/2)*x, 1/20*y]~, [(1/2*y + 5)*x, 1/20*y^2 - 3/20]~, [9/2*x + (1/2*y^2 + 1/2*y - 3/2), 1/20*x + (1/20*y^2 + 1/20*y - 3/20)]~, [(1/2*y + 3/2)*x + 1/2*y, (1/20*y + 1/20)*x + 1/20*y]~, [(81/74*y^2 + 90/37*y + 16/37)*x + (1/2*y^2 + y - 1), (1/2220*y^2 + 31/2220*y + 17/555)*x + (1/30*y^2 + 1/30*y - 1/12)]~]);/*Level a prime above 5*/
avgtime(fil, A, Or, tno, testspercase); tno++;/*area 2513.2741228718345907701147066236023074*/

setrand(1);
F = nfinit(y^3 - y^2 - 2*y + 1);
setrand(1);
A = alginit(F, [2*y^2 - 7*y - 7, -39*y^2 + 67*y - 49]);/*Ramified at primes above 13, 757*/
avgtime(fil, A, 0, tno, testspercase); tno++;/*area 1357.1680263507906790158619415767452460*/


setrand(1);
F = nfinit(y^4 + 2*y^3 - 4*y^2 - 6*y + 1);
setrand(1);
A = alginit(F, [-92*y^3 - 172*y^2 + 352*y + 221, 2*y^3 + 4*y^2 - 7*y - 18]);/*Ramified above 17*/
avgtime(fil, A, 0, tno, testspercase); tno++;/*area 469.14450293607579027708807856973909737;*/

setrand(1);
F = nfinit(y^4 - y^3 - 3*y^2 + y + 1);
setrand(1);
A = alginit(F, [-8*y^3 + 2*y^2 + 27*y - 3, -15*y^3 - 27*y^2 + 88*y - 105]);/*Ramified at prime above 601*/
Or = algorderalgtoorder(A, [[1, 0]~, [-y^3 + y^2 + 2*y - 1, 0]~, [-y, 0]~, [-y^2 + y + 2, 0]~, [11/2*x + (-11/2*y^3 + 11*y^2 + 11*y - 33/2), 0]~, [(1/2*y^3 - 1/2*y^2 - y + 4)*x + (-7/2*y^3 + 13/2*y^2 + 15/2*y - 19/2), 0]~, [(1/2*y + 1)*x + (-1/2*y^3 + 3/2*y^2 + y - 5/2), 0]~, [(1/2*y^2 - 1/2*y + 5/2)*x + (-3*y^3 + 6*y^2 + 6*y - 19/2), 0]~, [3/2*x + (-3/2*y^3 + 3*y^2 + 3*y - 9/2), 1]~, [5*x + (-5*y^3 + 10*y^2 + 10*y - 15), -y^3 + y^2 + 2*y - 1]~, [3*x + (-3*y^3 + 6*y^2 + 6*y - 9), -y]~, [5*x + (-5*y^3 + 10*y^2 + 10*y - 15), -y^2 + y + 2]~, [3*x + (-3*y^3 + 6*y^2 + 6*y - 9), 1/2*x + (-1/2*y^3 + y^2 + y - 3/2)]~, [4*x + (-4*y^3 + 8*y^2 + 8*y - 12), (1/2*y^3 - 1/2*y^2 - y + 1/2)*x + (-1/2*y^2 + 1/2*y + 1)]~, [5*x + (-5*y^3 + 10*y^2 + 10*y - 15), 1/2*y*x + (1/2*y^3 - 1/2*y^2 - y + 1/2)]~, [(4822/38719*y^3 - 3040/38719*y^2 + 2369/77438*y + 223209/77438)*x + (-3*y^3 + 11/2*y^2 + 6*y - 8), (12289960120/24867626221*y^3 - 24579920239/49735252442*y^2 - 33604592547/49735252442*y + 47018067775/49735252442)*x + (-224150/642259*y^3 + 130210/642259*y^2 + 508827/642259*y + 467903/1284518)]~]);/*Level prime above 11*/
avgtime(fil, A, Or, tno, testspercase); tno++;/*area 1507.9644737231007544620688239741613844*/

\p57
setrand(1);
F = nfinit(y^4 - 2*y^3 - 7*y^2 + 8*y + 1);
setrand(1);
A = alginit(F, [-1932*y^3 + 8092*y^2 + 15120*y - 59591, -7*y^3 - 14*y^2 + 91*y - 427]);/*Ramified above 229*/
avgtime(fil, A, 0, tno, testspercase); tno++;/*area 573.0265000147782866180001803968480089*/
\p38

setrand(1);
F = nfinit(y^4 - 9*y^2 + 4);
setrand(1);
A = alginit(F, [788*y^3 - 632*y^2 - 6692*y + 1128, -y^3 + y^2 + 8*y - 24]);/*Ramified above 13*/
avgtime(fil, A, 0, tno, testspercase); tno++;/*area 562.97340352329094833250569428368691684*/


setrand(1);
F = nfinit(y^5 - y^4 - 4*y^3 + 3*y^2 + 3*y - 1);
setrand(1);
A = alginit(F, [4*y^4 + 8*y^3 + 8*y^2 + 4*y - 79, -1]);/*Ramified above 11 and 131*/
avgtime(fil, A, 0, tno, testspercase); tno++;/*area 618.79855297980775909112672700959905294*/

\p57
setrand(1);
F = nfinit(y^5 - 5*y^3 - y^2 + 3*y + 1);
setrand(1);
A = alginit(F, [-104*y^4 + 116*y^3 + 672*y^2 - 660*y - 867, y^2 - y - 10]);/*Ramified above 5 and 83*/
avgtime(fil, A, 0, tno, testspercase); tno++;/*area 343.48079679248406073858234323855898197*/
\p38

setrand(1);
F = nfinit(y^5 - 5*y^3 + 4*y - 1);
setrand(1);
A = alginit(F, [16*y^4 - 20*y^3 - 56*y^2 + 64*y - 31, -1]);/*Ramified above 7 and 59*/
avgtime(fil, A, 0, tno, testspercase); tno++;/*area 728.84949563283203132333326482460490553*/

\p57
setrand(1);
F = nfinit(y^5 - 2*y^4 - 5*y^3 + 9*y^2 + 5*y - 7);
setrand(1);
A = alginit(F, [8*y^4 - 104*y^3 + 164*y^2 + 188*y - 615, y - 11]);/*Ramified above 5 and 7*/
avgtime(fil, A, 0, tno, testspercase); tno++;/*area 829.38046054770541495413785318578876548*/


\p96
setrand(1);
F = nfinit(y^6 - y^5 - 7*y^4 + 2*y^3 + 7*y^2 - 2*y - 1);
setrand(1);
A = alginit(F, [-y^5 + 2*y^4 + 3*y^3 - y^2 - 12, 244*y^5 - 192*y^4 - 1396*y^3 - 596*y^2 + 696*y - 291]);/*Ramified above 281*/
avgtime(fil, A, 0, tno, testspercase); tno++;/*area 309.97047515419293286164748048357761790*/
\p38

setrand(1);
F = nfinit(y^6 - y^5 - 5*y^4 + 4*y^3 + 6*y^2 - 3*y - 1);
setrand(1);
A = alginit(F, [4*y^4 - 12*y^3 + 12*y - 47, -1]);/*Ramified above 131*/
avgtime(fil, A, 0, tno, testspercase); tno++;/*area 198.96753472735357176930074760770184933*/

\p96
setrand(1);
F = nfinit(y^6 - 6*y^4 - y^3 + 8*y^2 + y - 2);
setrand(1);
A = alginit(F, [-252*y^5 + 176*y^4 + 1048*y^3 - 196*y^2 - 834*y - 627, y^5 - 4*y^3 - y^2 + 2*y - 24]);/*Ramified above 2*/
avgtime(fil, A, 0, tno, testspercase); tno++;/*area 714.18872991607966287717426246554032234*/
\p38

\p57
setrand(1);
F = nfinit(y^6 - 6*y^4 + 8*y^2 - y - 1);
setrand(1);
A = alginit(F, [-28*y^5 + 64*y^4 + 52*y^3 - 268*y^2 + 84*y - 67, y^4 - y^3 - 4*y^2 + 2*y - 6]);/*Ramified above 13*/
avgtime(fil, A, 0, tno, testspercase); tno++;/*area 326.72563597333849680011491186106829978*/
\p38


\p115
setrand(1);
F = nfinit(y^7 - y^6 - 6*y^5 + 4*y^4 + 10*y^3 - 4*y^2 - 4*y + 1);
setrand(1);
A = alginit(F, [-627*y^6 - 4180*y^5 + 9614*y^4 + 18601*y^3 - 21109*y^2 - 18183*y + 627, 1045*y^6 - 2926*y^5 + 2299*y^4 + 2926*y^3 - 11077*y^2 - 836*y - 23198]);/*Unramified*/
avgtime(fil, A, 0, tno, testspercase); tno++;/*area 15.70796326794802006809187444602*/
\p38

\p57
setrand(1);
F = nfinit(y^7 - 8*y^5 + 19*y^3 - y^2 - 13*y + 1);
setrand(1);
A = alginit(F, [-4*y^5 + y^4 + 20*y^3 - 2*y^2 - 16*y - 3, -10*y^6 - 2*y^5 + 58*y^4 + 21*y^3 - 74*y^2 - 20*y - 15]);/*Unramified*/
avgtime(fil, A, 0, tno, testspercase); tno++;/*area 138.23007675794148146897200035709829535*/
\p38

setrand(1);
F = nfinit(y^7 - 7*y^5 + 13*y^3 - 5*y - 1);
setrand(1);
A = alginit(F, [4*y^6 - 20*y^4 - 8*y^3 + 20*y^2 + 12*y - 27, -1]);/*Ramified above 3, 11*/
avgtime(fil, A, 0, tno, testspercase); tno++;/*area 1361.3568165555770700005059216922265252*/

\p57
setrand(1);
F = nfinit(y^7 - 2*y^6 - 6*y^5 + 8*y^4 + 12*y^3 - 5*y^2 - 6*y - 1);
setrand(1);
A = alginit(F, [-4*y^6 + 18*y^5 - 11*y^4 - 37*y^3 + 25*y^2 + 22*y - 1, -21*y^6 + 17*y^5 + 176*y^4 - 20*y^3 - 400*y^2 - 183*y - 33]);
avgtime(fil, A, 0, tno, testspercase); tno++;/*area 238.76104167282428612316089712924221924*/
\p38


setrand(1);
F = nfinit(y^8 - 4*y^7 - y^6 + 17*y^5 - 5*y^4 - 23*y^3 + 6*y^2 + 9*y - 1);
setrand(1);
A = alginit(F, [-1, -16*y^7 + 68*y^6 - 20*y^5 - 204*y^4 + 168*y^3 + 120*y^2 - 140*y - 3]);/*Ramified above 19*/
avgtime(fil, A, 0, tno, testspercase); tno++;/*area 422.23005264246821124937542704619344205*/

setrand(1);
F = nfinit(y^8 - 2*y^7 - 7*y^6 + 11*y^5 + 14*y^4 - 18*y^3 - 8*y^2 + 9*y - 1);
setrand(1);
A = alginit(F, [-1, -12*y^7 + 16*y^6 + 76*y^5 - 48*y^4 - 96*y^3 + 28*y^2 + 16*y - 35]);/*Ramified above 11*/
avgtime(fil, A, 0, tno, testspercase); tno++;/*area 423.06781068342548944634012072981623082*/


/*End of tests*/
filewrite(fil, "End of tests\n\n");
fileclose(fil);
print("End of tests\n\n");
#/*Turn the timer back on*/