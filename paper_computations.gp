\p80

\\This file contains methods used to generate data in the paper "IMPROVED COMPUTATION OF FUNDAMENTAL DOMAINS FOR ARITHMETIC FUCHSIAN GROUPS".

\\FIGURES 7 AND 8
\\Make plots for the success rate heuristic.
successrate(n, compile, WSL)={
  my(F, A);
  if(n==1,
    F=nfinit(y);
	A=alginit(F, [589, -1]);
	return(1);
  );
  if(n==2,
    F=nfinit(y^2-21);\\disc(F)=21
	A=alginit(F, [-10*y + 18, 64*y - 647]);\\algnormdisc(A)=101
	return(enum_successrate_range(A, I/2, 2, 20, 400, 1000, 0, "d2_success", compile, WSL));
  );
  if(n==3,
    F=nfinit(y^3 - 4*y - 1);\\disc(F)=229
	A=alginit(F, [12*y^2 - 22*y - 26, 410*y^2 + 106*y - 2213]);\\algnormdisc(A)=404
	return(enum_successrate_range(A, I/2, 3, 25, 400, 1000, 0, "d3_success", compile, WSL));
  );
  if(n==4,
    F=nfinit(y^4 - y^3 - 7*y^2 + 3*y + 9);\\disc(F)=4525
	A=alginit(F, [-12*y^3 + 30*y^2 + 30*y - 144, 48*y^3 + 96*y^2 - 480*y - 1359]);\\algnormdisc(A)=149
	return(enum_successrate_range(A, I/2, 4, 20, 400, 1000, 0, "d4_success", compile, WSL));
  );
}

\\FIGURES 9 AND 10
\\Make plots for the time to compute small norm 1 elements.
Ctime(n, compile, WSL)={
  my(F, A);
  if(n==1,
    F=nfinit(y);\\disc(F)=1
    A=alginit(F, [2021, -1]);\\algnormdisc(A)=2021
	return(enum_time_range(A, I/2, 5, 5000, 1000, 200, "d1_Ctime", compile, WSL));
  );
  if(n==2,
    F=nfinit(y^2-163);\\disc(F)=652
    A=alginit(F, [-7*y + 10, -259*y - 7516]);\\algnormdisc(A)=114
	return(enum_time_range(A, I/2, 50, 500, 1000, 200, "d2_Ctime", compile, WSL));
  );
  if(n==3,
    F=nfinit(y^3 - 12*y - 3);\\disc(F)=6669
	A=alginit(F, [y^2 - y - 10, -4*y - 3]);\\algnormdisc(A)=51
	return(enum_time_range(A, I/2, 20, 100, 1000, 200, "d3_Ctime", compile, WSL));
  );
  if(n==4,
    F=nfinit(y^4 - y^3 - 6*y^2 + 6*y + 3);\\disc(F)=14013
    A=alginit(F, [-2*y^3 + 5*y - 7, -48*y^3 - 48*y^2 - 4*y - 403]);\\algnormdisc(A)=109
	return(enum_time_range(A, I/2, 15, 45, 1000, 200, "d4_Ctime", compile, WSL));
  );
}

\\FIGURES 11 AND 12
\\Make plots to vary the algebra over a fixed field of degree n.
Arange(n, compile, WSL)={
  my(F, ablist, Avec);
  if(n==1,
    F=nfinit(y);\\disc(F)=1
	ablist=readvec("data_in/Arangedata.txt")[1];\\length 400 list of ab's for F
	Avec=vector(400, i, alginit(F, ablist[i]));
	return(enum_bestC_range(Avec, I/2, 4, 20, 200, "d1_Arange", 1, compile, WSL));
  );
  if(n==2,
    F=nfinit(y^2-13);\\disc(F)=13
	ablist=readvec("data_in/Arangedata.txt")[2];\\length 400 list of ab's for F
	Avec=vector(400, i, alginit(F, ablist[i]));
	return(enum_bestC_range(Avec, I/2, 4, 20, 200, "d2_Arange", 1, compile, WSL));
  );
  if(n==3,
    F=nfinit(y^3 - y^2 - 4*y - 1);\\disc(F)=169
	ablist=readvec("data_in/Arangedata.txt")[3];\\length 400 list of ab's for F
	Avec=vector(400, i, alginit(F, ablist[i]));
	return(enum_bestC_range(Avec, I/2, 4, 20, 200, "d3_Arange", 1, compile, WSL));
  );
  if(n==4,
    F=nfinit(y^4 - 5*y^2 + 5);\\disc(F)=2000
	ablist=readvec("data_in/Arangedata.txt")[4];\\length 400 list of ab's for F
	Avec=vector(400, i, alginit(F, ablist[i]));
	return(enum_bestC_range(Avec, I/2, 4, 20, 200, "d4_Arange", 1, compile, WSL));
  );
}

\\Computes sample algebras of degree n. We take numfield fields, and numalgperfield algebras per field. The pairs [Fs, As] are written to "./data_in/filename.txt", and this produces output valid for Frange. A list of small fields must be found in the file pointed to by smallfields
getsamplealgebras(n, numfield, numalgperfield, filename, smallfieldsfile)={
  my(v, totN, Fs, As, F, u, k, outv);
  if(smallfieldsfile==0,smallfieldsfile=strprintf("../DATA/TotallyReal/d%d_3000.txt", n));\\Default small field file
  v=readvec(smallfieldsfile);
  totN=numfield*numalgperfield;
  Fs=vector(totN);\\Holds the fields
  for(i=1,numfield,
    for(j=numalgperfield*(i-1)+1,numalgperfield*i,Fs[j]=v[3][i]);\\Outputting the fields
  );
  As=vector(totN);\\Holds the algebras
  for(i=1,numfield,
    F=nfinit(Fs[numalgperfield*i]);
	u=smallalgebras(F, 1, 2);
	k=numalgperfield*(i-1)+1;
	As[k]=u[2][1];
    for(j=2,numalgperfield,
	  k++;
	  u=smallalgebras(F, 1, u[1][1]+1);
	  As[k]=u[2][1];
	);
	printf("Field number %d done.\n", i);
  );
  st=strprintf("./data_in/%s.txt", filename);
  outv=[Fs, As];
  write(st, outv);
  return(outv);
}

\\FIGURES 13 AND 14
\\Make plots to vary the field and algebra over a fixed degree n.
Frange(n, compile, WSL)={
  my(dat, len, Avec, F);
  dat=readvec("data_in/Frangedata.txt")[n];\\[Fields, algebras], length 400 list for n<=7, 165 for n=8.
  len=length(dat[1]);
  Avec=vector(len);
  for(i=1,len,
    F=nfinit(dat[1][i]);
    Avec[i]=alginit(F, dat[2][i]);
    if(i%50==0,printf("Algebra %d initialized.\n", i));
  );
  if(n==2,return(enum_bestC_range(Avec, I/2, 20, 60, 200, "d2_Frange", 0, compile, WSL)));
  if(n==3,return(enum_bestC_range(Avec, I/2, 20, 60, 200, "d3_Frange", 0, compile, WSL)));\\27mins gpu1
  if(n==4,return(enum_bestC_range(Avec, I/2, 20, 60, 200, "d4_Frange", 0, compile, WSL)));\\28mins gpu1 
  if(n==5,return(enum_bestC_range(Avec, I/2, 20, 60, 200, "d5_Frange", 0, compile, WSL)));\\29mins gpu1
  if(n==6,return(enum_bestC_range(Avec, I/2, 20, 60, 200, "d6_Frange", 0, compile, WSL)));\\56mins gpu1
  if(n==7,return(enum_bestC_range(Avec, I/2, 20, 60, 200, "d7_Frange", 0, compile, WSL)));\\78mins gpu1
  if(n==8,return(enum_bestC_range(Avec, I/2, 20, 60, 200, "d8_Frange", 0, compile, WSL)));\\
}

\\FIGURES 15, 16, 17, 18, 19, 20; ty=1 corresponds to FP, ty=2 corresponds to IFP
\\Computing the time taken to find N non-trivial elements.
timeforNelts(n, ty)={
  my(F, A);
  if(n==1,
    F=nfinit(y);
	A=alginit(F, [20, -119]);
	enum_timeforNelts_range(A, I/2, 5, 60, 1000, 1000, "d1_Nelttime", ty);
	A=alginit(F, [3, 2827]);
	enum_timeforNelts_range(A, I/2, 5, 120, 1000, 1000, "d1_Nelttime2", ty);
  );
  if(n==2,
    F=nfinit(y^2-5);
	A=alginit(F, [-2*y - 1, 2*y - 90]);
	enum_timeforNelts_range(A, I/2, 3, 13, 500, 1000, "d2_Nelttime", ty);
  );
  if(n==3,
    F=nfinit(y^3 - y^2 - 4*y + 2);
	A=alginit(F, [5*y^2 - 19*y - 1, -82*y^2 + 200*y - 863]);
	enum_timeforNelts_range(A, I/2, 5, 20, 500, 1000, "d3_Nelttime", ty);
  );
}

\\FIGURES 21 AND 22
\\fin is a file with a bunch of vectors. Each vector is [pol, ab1, ab2, ...] where F=nfinit(pol) and A=alginit(F, abi). This computes the fundamental domain, the number of sides, the area, and the number of generators. fin is in the folder "data_in", and fout is output in the folder "data".
sides_and_gens(fin, fout)={
  my(st, v, f, F, A, W);
  st=strprintf("data_in/%s.txt", fin);
  v=readvec(st);
  st=strprintf("data/%s.txt", fout);
  f=fileopen(st, "w");
  filewrite(f, "nelts Sides Area\n");
  printf("i j nelts Sides Area\n");
  for(i=1,length(v),
    F=nfinit(v[i][1]);
	for(j=2,length(v[i]),
	  A=alginit(F, v[i][j]);
	  W=algfdom_nelts(A, I/2);
	  st=strprintf("%d %d %.3f\n", W[1], W[2], W[3]);
	  filewrite1(f, st);
	  printf("%d %d: %s", i, j, st);
	);
  );
  fileclose(f);
}

\\FIGURES 23, 24, 25, 26
\\Set precision to >=60 first. This computes the data in the last 4 figures of my paper.
fdometime(n)={
  if(n==1,gtimedata("fdtime_d1_1to20000", "fdtime_d1_1to20000_results"));\\areas 1->20000
  if(n==2,gtimedata("fdtime_d2_1to4000", "fdtime_d2_1to4000_results"));\\areas 1->4000, 3 fields
  if(n==3,gtimedata("fdtime_d3_1to2500", "fdtime_d3_1to2500_results"));\\areas 1->2500, 4 fields
  if(n==4,gtimedata("fdtime_d4_1to1500", "fdtime_d4_1to1500_results"));\\areas 1->1500, 6 fields
}

\\fin is a file with a bunch of vectors. Each vector is [pol, ab1, ab2, ...] where F=nfinit(pol) and A=alginit(F, abi). This computes the time taken to compute the fundamental domain for each of these algebras. fin is in the folder "data_in", and fout is output in the folder "data".
gtimedata(fin, fout)={
  my(st, v, f, F, A, Fdisc, n, W, T Adisc);
  st=strprintf("data_in/%s.txt", fin);
  v=readvec(st);
  st=strprintf("data/%s.txt", fout);
  f=fileopen(st, "w");
  filewrite(f, "n disc(F) N(D) Sides Area Time\n");
  printf("i j n disc(F) N(D) Sides Area Time\n");
  for(i=1,length(v),
    F=nfinit(v[i][1]);
	Fdisc=F.disc;
	n=poldegree(F.pol);
	for(j=2,length(v[i]),
	  A=alginit(F, v[i][j]);
	  Adisc=algnormdisc(A);
	  gettime();
	  W=algfdom(A, , I/2, 0, 0, 0);
	  T=gettime()/1000.;
	  st=strprintf("%d %d %d %d %.3f %.3f\n", n, Fdisc, Adisc, length(W[1]), W[6], T);
	  filewrite1(f, st);
	  printf("%d %d: %s", i, j, st);
	);
  );
  fileclose(f);
}

