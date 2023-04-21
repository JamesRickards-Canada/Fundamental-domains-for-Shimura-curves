/*These methods are useful with the fdom.c package, but not essential.*/

/*INCLUSIONS*/

#ifndef PARILIB
#define PARILIB
#include <pari/pari.h>
#endif

#ifndef FDOMDECL
#define FDOMDECL
#include "fdom.h"
#endif

/*STATIC DECLARATIONS*/

/*SECTION 1: VISUALIZATION*/
/*SECTION 2: TUNING*/

/*2: BEST C*/
static GEN tune_bestC_givenAB(long n, GEN A, GEN B, long prec);
static void tune_bestC_plot(GEN reg, GEN Cmin, GEN Cmax, long n, char *fdata);
static GEN tune_time(GEN X, GEN Cset, long mintesttime, long prec);

/*2: REGRESSIONS AND PLOTS*/
static GEN OLS(GEN X, GEN y, int retrsqr);
static GEN OLS_nointercept(GEN X, GEN y, int retrsqr);
static GEN rsquared(GEN X, GEN y, GEN fit);
static void plot_compile(char *fname, int WSL);


/*MAIN BODY*/

/*SECTION 1: VISUALIZATION*/

/*Prints the fundamental domain to a latex file. model=0 means upper half plane, 1 is unit disc, and 2 is Klein.*/
void
afuchfdom_latex(GEN X, char *filename, int model, int boundcircle, int compile, int open)
{
  pari_sp av = avma;
  if (!model) pari_err(e_MISC, "Upper half plane not yet supported");
  GEN tol = gdat_get_tol(afuch_get_gdat(X));
  long prec = lg(tol);
  GEN U = afuch_get_fdom(X);
  GEN width = dbltor(6.0);/*Width of the picture in inches.*/
  GEN radius = shiftr(width, -1);/*Circle radius*/
  if (!pari_is_dir("plots/build")) {/*Checking the directory*/
      int s = system("mkdir -p plots/build");
      if (s == -1) pari_err(e_MISC, "ERROR CREATING DIRECTORY");
  }
  char *plotmake = stack_sprintf("plots/build/%s.tex", filename);
  FILE *f = fopen(plotmake, "w");
  pari_fprintf(f, "\\documentclass{standalone}\n\\usepackage{pgf}\n\\standaloneenv{pgfpicture}\n");/*Initial things*/
  pari_fprintf(f, "%%Constants: update these to your liking\n%%Colours\n\\definecolor{background}{rgb}{1,1,1}%%White\n");/*Some constants*/
  pari_fprintf(f, "\\definecolor{circle_edge}{rgb}{0,0,1}%%Blue\n\\definecolor{fdom_edge}{rgb}{0,0.5,0}%%Green\n");
  pari_fprintf(f, "\\definecolor{fdom_fill}{rgb}{0.85,0.85,0.85}%%Grey\n%%Lengths\n");
  pari_fprintf(f, "\\def\\circlewidth{0.02in}\n\\def\\fdomwidth{0.02in}\n\n\\begin{document}\n\\begin{pgfpicture}\n");
  /*Start with the background*/
  pari_fprintf(f, "%%Background\n\\pgfsetfillcolor{background}\n\\pgfpathrectangle{\\pgfpoint{-%P.8fin-0.5*\\circlewidth}{-%P.8fin-0.5*\\circlewidth}}{\\pgfpoint{%P.8fin+\\circlewidth}{%P.8fin+\\circlewidth}}\n\\pgfusepath{fill}\n\n", radius, radius, width, width);
  
  /*The circle*/
  if(!boundcircle) pari_fprintf(f,"\\iffalse\n");
  pari_fprintf(f, "%%Circle\n\\begin{pgfscope}\n\\pgfsetlinewidth{\\circlewidth}\n\\pgfsetstrokecolor{circle_edge}\n");/*Starting it*/
  pari_fprintf(f, "\\pgfpathmoveto{\\pgfpoint{%P.8fin}{0}}\n\\pgfpatharc{0}{360}{%P.8fin}\n", radius, radius);/*The circle*/
  pari_fprintf(f, "\\pgfpathclose\n\\pgfusepath{stroke}\n\\end{pgfscope}\n");
  if(!boundcircle) pari_fprintf(f,"\\fi\n");
  
  /*Fundamental domain*/
  pari_fprintf(f, "\n%%Fundamental domain\n\\begin{pgfscope}\n\\pgfsetroundjoin %%Rounded intersections\n");
  pari_fprintf(f, "\\pgfsetlinewidth{\\fdomwidth}\n\\pgfsetstrokecolor{fdom_edge}\n\\pgfsetfillcolor{fdom_fill}\n");
  GEN circles = normbound_get_sides(U);
  GEN vertices = normbound_get_vcors(U);
  long nsides = lg(circles) - 1, i;
  if (model == 1) {/*Unit ball model*/
    GEN first = klein_to_disc(gel(vertices, nsides), tol);/*The lower point on the first arc.*/
    pari_fprintf(f, "\\pgfpathmoveto{\\pgfqpoint{%P.8fin}{%P.8fin}}\n", gmul(real_i(first), radius), gmul(imag_i(first), radius));/*First point*/
    GEN ascale = divsr(180, mppi(prec));/*Radian to degrees*/
    GEN lastv = first;
    for (i = 1; i <= nsides; i++) {/*Each arc*/
      GEN arc = gel(circles, i);
	  GEN nextv = klein_to_disc(gel(vertices, i), tol);
	  GEN centre = mkcomplex(gel(arc, 1), gel(arc, 2));
	  GEN ang1 = garg(gsub(lastv, centre), prec);
	  GEN ang2 = garg(gsub(nextv, centre), prec);
	  if (gcmp(ang1, ang2) < 0) ang1 = gadd(ang1, Pi2n(1, prec));
      pari_fprintf(f, "\\pgfpatharc{%P.8f}{%P.8f}{%P.8fin}\n", gmul(ang1, ascale), gmul(ang2, ascale), gmul(gel(arc, 3), radius));
      lastv=nextv;/*Updating the previous vertex.*/
    }
    pari_fprintf(f, "\\pgfpathclose\n\\pgfusepath{stroke, fill}\n\\end{pgfscope}\n");
  }
  else if (model == 2) {/*Klein model*/
    GEN first = gel(vertices, nsides);/*The lower point on the first arc.*/
    pari_fprintf(f, "\\pgfpathmoveto{\\pgfqpoint{%P.8fin}{%P.8fin}}\n", gmul(real_i(first), radius), gmul(imag_i(first), radius));/*First point*/
    for (i = 1; i <= nsides; i++) {/*Each side*/
	  GEN ver = gel(vertices, i);
      pari_fprintf(f, "\\pgfpathlineto{\\pgfqpoint{%P.8fin}{%P.8fin}}\n", gmul(real_i(ver), radius), gmul(imag_i(ver), radius));
    }
    pari_fprintf(f, "\\pgfpathclose\n\\pgfusepath{stroke, fill}\n\\end{pgfscope}\n");
  }
  /*End*/
  pari_fprintf(f, "\\end{pgfpicture}\n\\end{document}");
  fclose(f);
  if (!compile) { set_avma(av); return; }
  
  /*Compile and open*/
  char *line = stack_sprintf("(cd ./plots/build && pdflatex --interaction=batchmode -shell-escape %s.tex)", filename);/*Build*/
  int s = system(line);
  if (s == -1) pari_err(e_MISC, "ERROR EXECUTING COMMAND");
  line = stack_sprintf("mv -f ./plots/build/%s.pdf ./plots/", filename);/*Move the file*/
  s = system(line);
  if (s == -1) pari_err(e_MISC, "ERROR EXECUTING COMMAND");
  if (open) {
    line = stack_sprintf("cmd.exe /C start plots/%s.pdf", filename);/*Open the file*/
    s = system(line);
    if (s == -1) pari_err(e_MISC, "ERROR EXECUTING COMMAND");
  }
  set_avma(av);
}


/*SECTION 2: TUNING*/

/*2: BEST C*/

/*Finds and returns the optimal C for a given algebra. Actually, returns [A, B, C, R^2]. We do ntrials trials in a range [Cmin, scale^(1/2n)*Cmin], where this is centred at the conjectural best C value.*/
GEN
tune_bestC(GEN X, GEN scale, long ntrials, long mintesttime)
{
  pari_sp av = avma;
  GEN Cbest = afuch_get_bestC(X);
  GEN A = afuch_get_alg(X);
  GEN nf = alg_get_center(A);
  long n = nf_get_degree(nf), twon = n << 1, fourn = twon << 1;
  long prec = lg(gdat_get_tol(afuch_get_gdat(X))), i;
  GEN scalefact = gpow(scale, gdivgs(gen_1, fourn), prec);
  GEN Cmin = gdiv(Cbest, scalefact);
  GEN Cmax = gmul(Cbest, scalefact);
  if (ntrials <= 1) ntrials = 2;/*Make sure it's at least 2.*/
  GEN blen = gdivgs(gsub(Cmax, Cmin), ntrials - 1);/*length.*/
  GEN Clist = cgetg(ntrials + 1, t_VEC);
  GEN C = Cmin;
  gel(Clist, 1) = C;
  for (i = 2; i <= ntrials; i++) {
    C = gadd(C, blen);
    gel(Clist, i) = C;
  }
  GEN times = tune_time(X, Clist, mintesttime, prec);/*Executing the computation*/
  GEN Xdata = cgetg(ntrials + 1, t_MAT);
  for (i = 1; i <= ntrials; i++) gel(Xdata, i) = mkcol2(gen_1, gpowgs(gel(Clist, i), twon));/*The X data*/
  GEN reg = OLS(Xdata, times, 1);
  GEN Aco = gmael(reg, 1, 1), Bco = gmael(reg, 1, 2);
  if(gsigne(Aco) == -1 || gsigne(Bco) == -1) pari_err(e_MISC, "One of A, B are negative! This is not correct.");
  C = tune_bestC_givenAB(n, Aco, Bco, prec);
  return gerepilecopy(av, mkvec4(Aco, Bco, C, gel(reg, 2)));/*[A, B, C, R^2].*/
}

/*Returns the unique real solution >n to (2n-1)C^2n-2n^2c^(2n-1)=A/B, which should be the optimal value of C to use for a given algebra.*/
static GEN
tune_bestC_givenAB(long n, GEN A, GEN B, long prec)
{
  pari_sp av = avma;
  long twon = n << 1;
  long var = fetch_var(), i;/*The temporary variable number*/
  GEN P = cgetg(twon + 3, t_POL);
  P[1] = evalvarn(var);
  gel(P, 2) = gdiv(gneg(A), B);/*Constant coefficient.*/
  for (i = 3; i <= twon;i++) gel(P, i)=gen_0;/*0's in the middle. P[i] corresp to x^(i-2)*/
  gel(P, twon + 1) = stoi(-twon*n);
  gel(P, twon + 2) = stoi(twon - 1);
  P = normalizepol(P);/*Normalize it in place.*/
  GEN rts = realroots(P, mkvec2(stoi(n), mkoo()), prec);/*Real roots in [n, oo)*/
  delete_var();/*Delete the variable.*/
  if (lg(rts) != 2) pari_err(e_MISC, "Could not find exactly one real root in [n, oo)!");
  return gerepilecopy(av, gel(rts, 1));
}

/*Prepares a basic latex plot of the data.*/
static void
tune_bestC_plot(GEN reg, GEN Cmin, GEN Cmax, long n, char *fdata)
{
  pari_sp av = avma;
  if (!pari_is_dir("plots/build")) {/*Checking the directory*/
      int s = system("mkdir -p plots/build");
      if (s == -1) pari_err(e_MISC, "ERROR CREATING DIRECTORY");
  }
  char *plotmake = stack_sprintf("plots/build/%s_plotter.tex", fdata);
  FILE *f = fopen(plotmake, "w");
  long invpower = n;
  pari_fprintf(f, "\\documentclass{article}\n\\usepackage{amsmath, amssymb, pgfplots}\n  \\usepgfplotslibrary{external}\n  \\tikzexternalize\n");
  pari_fprintf(f, "  \\pgfplotsset{compat=1.16}\n\\begin{document}\n\\tikzsetnextfilename{%s}\n\\begin{tikzpicture}\n  \\begin{axis}", fdata);
  pari_fprintf(f, "[xlabel=$\\text{disc}(F)$, ylabel=$\\text{C}/\\text{Nm}_{F/\\mathbb{Q}}(\\mathfrak{D})^{1/%d}$,\n", invpower << 1);
  pari_fprintf(f, "    xmin=0, ymin=0,\n");
  pari_fprintf(f, "    scatter/classes={a={mark=o}}, clip mode=individual,]\n");
  pari_fprintf(f, "    \\addplot[scatter, blue, only marks, mark size=0.9]\n      table[x=x,y=y,col sep=space]{%s.dat};\n", fdata);
  pari_fprintf(f, "    \\addplot[red, ultra thick, samples=1000, domain=%Pf:%Pf]{%Pf*(x)^(1/%d)}", Cmin, Cmax, gel(reg, 1), invpower);
  pari_fprintf(f, ";%%R^2=%Pf\n  \\end{axis}\n\\end{tikzpicture}\n\\end{document}", gel(reg, 2));
  fclose(f);
  set_avma(av);
}

/*Does tune_bestC for all algebras in Aset, where we assume n=deg(F) is fixed, and we regress along disc(F).*/
GEN
tune_bestC_range(GEN Aset, GEN scale, long ntrials, long mintesttime, char *fname, int compile, int WSL, long prec)
{
  pari_sp av = avma;
  if (!pari_is_dir("plots/build")) {/*Checking the directory*/
    int s = system("mkdir -p plots/build");
    if (s == -1) pari_err(e_MISC, "ERROR CREATING DIRECTORY");
  }
  char *fname_full = stack_sprintf("plots/build/%s.dat", fname);
  FILE *f = fopen(fname_full, "w");
  pari_fprintf(f, "x y rsqr\n");

  GEN nf = alg_get_center(gel(Aset, 1));
  long n = nf_get_degree(nf);
  GEN oneovertwon = mkfracss(1, n << 1), oneovern = (n == 1)? gen_1: mkfracss(1, n);
  long lgAset = lg(Aset);
  GEN Xdat = vectrunc_init(lgAset), Cdat = coltrunc_init(lgAset), lastX = gen_0, firstX = NULL;
  long i;
  for (i = 1; i < lgAset; i++) {
    GEN A = gel(Aset, i);
	GEN X = afuchinit(A, NULL, gen_0, NULL, 0, prec);
    GEN C;
    pari_CATCH (CATCH_ALL) {
      C = NULL;
    }
    pari_TRY {C = tune_bestC(X, scale, ntrials, mintesttime);} pari_ENDCATCH
    if (!C) {
      pari_printf("Algebra %d skipped due to error.\n", i);
      continue;
    }
    nf = alg_get_center(A);
    GEN nfdisc = nf_get_disc(nf);
    GEN Adisc = algdiscnorm(A);
    gel(C, 3) = gdiv(gel(C, 3), gpow(Adisc, oneovertwon, prec));
    pari_fprintf(f, "%Pd %Pf %Pf\n", nfdisc, gel(C, 3), gel(C, 4));
    if(!firstX) firstX = nfdisc;
    lastX = nfdisc;
    vectrunc_append(Xdat, gpow(nfdisc, oneovern, prec));
    vectrunc_append(Cdat, gel(C, 3));
    pari_printf("Algebra %d done.\n", i);
  }
  fclose(f);
  if (n == 1) {/*nfdisc=1, so no regression: just do average and variance.*/
	GEN S = vecsum(Cdat);
	long lgC = lg(Cdat);
	GEN Cavg = gdivgs(S, lgC - 1);
	GEN var = gen_0;
	for (i = 1; i < lgC; i++) var = gadd(var, gsqr(gsub(gel(Cdat, i), Cavg)));
	var = gdivgs(var, lgC - 1);
	return gerepilecopy(av, mkvec2(Cavg, var));
  }
  GEN reg = OLS_nointercept(Xdat, Cdat, 1);
  tune_bestC_plot(reg, firstX, lastX, n, fname);
  if(compile) plot_compile(fname, WSL);
  return gerepilecopy(av, reg);
}

/*Computes the average time to find algsmallnormelts(A, C, 0, z) for all C in Cset, and returns it as a column vector.*/
static GEN
tune_time(GEN X, GEN Cset, long mintesttime, long prec)
{
  pari_sp av = avma, av2;
  GEN R = afuch_get_R(X);
  long lgCset = lg(Cset), i;
  GEN avgtimes = cgetg(lgCset, t_COL);
  pari_timer T;
  timer_start(&T);
  for (i = 1; i < lgCset; i++){
    long t = 0, tries = 0;
    timer_delay(&T);
    while (t < mintesttime) {/*Make sure we do at least mintesttime per test*/
      av2 = avma;
      tries++;
      GEN z = hdiscrandom(R, prec);
	  GEN elts = afuchfindelts_i(X, z, gel(Cset, i), 0);
      if(!elts) pari_err_PREC("Precision too low");
      t = timer_get(&T);
      set_avma(av2);
    }
    gel(avgtimes, i) = rdivss(t, 1000*tries, prec);
    if (i%10 == 0) pari_printf("%d test cases done.\n", i);
  }
  return gerepilecopy(av, avgtimes);
}


/*We have pre-stored algebras of degree n in "data_in/fdom_Cn.dat". We compute the time to compute the fundamental domains of the algebras, where each test is repeated testsperalg times, with the range of C_n's.*/
GEN
tune_Cnrange(long n, GEN Cmin, GEN Cmax, long testsperalg, long tests, long prec)
{
  pari_sp av = avma, av1, av2, av3;
  if (n <= 0 || n >= 9) pari_err(e_MISC, "n must be between 1 and 8.");
  GEN dat = gel(gp_readvec_file("data_in/fdom_Cn.dat"), n);
  long ldat = lg(gel(dat, 1)), i, j, k;
  if (tests <= 1) tests = 2;
  GEN Cn = Cmin;
  GEN Cnadd = gdivgs(gsub(Cmax, Cmin), tests - 1);
  GEN times = cgetg(tests + 1, t_VECSMALL);
  pari_timer T;
  timer_start(&T);
  av1 = avma;
  for (i = 1; i <= tests; i++) {
	long t = 0;
	av2 = avma;
	for (j  = 1; j < ldat; j++) {
	  GEN F = nfinit(gmael(dat, 1, j), prec);
	  GEN A = alginit(F, gmael(dat, 2, j), 0, 1);
      GEN Adisc = algdiscnorm(A);/*Norm to Q of disc(A)*/
	  GEN discpart = gmul(nf_get_disc(F), gsqrt(Adisc, prec));/*disc(F)*sqrt(Adisc)*/
      GEN discpartroot = gpow(discpart, gdivgs(gen_1, n), prec);/*discpart^(1/n)=disc(F)^(1/n)*algdisc^(1/2n)*/
	  GEN C = gmul(Cn, discpartroot);
      if (gcmpgs(C, n + 2) <= 0) C = stoi(n + 2);
	  av3 = avma;
	  for (k = 1; k <= testsperalg; k++) {
		pari_printf("%P.8f %d %d %d\n", C, i, j, k);
	    timer_delay(&T);
	    GEN X = afuchinit(A, NULL, NULL, NULL, 0, prec);
	    gmael3(X, 7, 6, 2) = C;/*This isn't really safe but should be OK for now. If we change where C is stored, this must change.*/
	    afuchfdom(X);
	    t = t + timer_delay(&T);
		set_avma(av3);
	  }
	  set_avma(av2);
	}
	times[i] = t;
	pari_printf("The value C_n=%P.8f took %d time.\n", Cn, t);
	Cn = gerepileupto(av1, gadd(Cn, Cnadd));
  }
  Cn = Cmin;
  GEN Cs = cgetg(tests + 1, t_VEC);
  for (i = 1; i <= tests; i++) {
	gel(Cs, i) = Cn;
	Cn = gadd(Cn, Cnadd);
  }
  return gerepilecopy(av, mkvec2(Cs, times));
}


/*2: TIME FOR N ELTS*/

/*Returns the time taken to find nelts non-trivial elements*/
long
tune_Nelts(GEN X, GEN C, long nelts, long prec)
{
  pari_sp av = avma, av2;
  GEN R = afuch_get_R(X);
  long found = 0;
  pari_timer T;
  timer_start(&T);
  while (found < nelts) {
	av2 = avma;
	GEN z = hdiscrandom(R, prec);
	GEN v = afuchfindelts_i(X, z, C, 1);
	if (!v) pari_err_PREC("Precision is too low, please increase.");
	if (lg(v) > 1) found++;
	set_avma(av2);
  }
  long t = timer_get(&T);
  return gc_long(av, t);
}



/*2: REGRESSIONS AND PLOTS*/

/*Perform ordinary least squares regression. X is a matrix whose columns are the parameters, and y is a column vector of results. Must include linear term as first variable of X.
The formula is B=Bhat=(X*X^T)^(-1)Xy, for the ordinary least squares regression for y=X^T*B+error (formula differs to Wikipedia due to X here being the transpose of what they define there.
Returns [best fit, R^2]*/
static GEN
OLS(GEN X, GEN y, int retrsqr)
{
  pari_sp av = avma;
  if (typ(y) != t_COL) y = gtocol(y);
  if (lg(y) != lg(X)) pari_err_TYPE("The inputs must have the same length.", mkvec2(X, y));
  GEN Xy = RgM_RgC_mul(X, y);
  GEN XTX = RgM_multosym(X, shallowtrans(X));/*X*X^T, which is symmetric*/
  GEN fit = RgM_solve(XTX, Xy);/*Best fit.*/
  if (!fit) pari_err(e_MISC, "Could not compute matrix inverse. Error with given data or with precision?");
  if (!retrsqr) return gerepileupto(av, fit);
  GEN rsqr = rsquared(X, y, fit);
  return gerepilecopy(av, mkvec2(fit, rsqr));
}

GEN oln(GEN X, GEN y){return OLS_nointercept(X, y, 1);}

/*Performs OLS where we have one independant variable and assume the intercept is 0 (so y=ax). The formua is now sum(x*y)/sum(x^2).*/
static GEN
OLS_nointercept(GEN X, GEN y, int retrsqr)
{
  pari_sp av = avma;
  GEN xysum = gen_0, xsqrsum = gen_0;
  long lX = lg(X), i;
  if (lX != lg(y)) pari_err_TYPE("The inputs must have the same length.", mkvec2(X, y));
  for (i = 1; i < lX; i++) {
    xysum = gadd(xysum, gmul(gel(X, i), gel(y, i)));
    xsqrsum = gadd(xsqrsum, gsqr(gel(X, i)));
  }
  GEN fit = gdiv(xysum, xsqrsum);
  if (!retrsqr) return gerepileupto(av, fit);
  GEN M = cgetg(lX, t_MAT);
  for (i = 1; i < lX; i++) gel(M, i) = mkcol2(gen_1, gel(X, i));
  GEN rsqr = rsquared(M, y, mkcol2(gen_0, fit));
  return gerepilecopy(av, mkvec2(fit, rsqr));
}

/*Given inputs for OLS and the proposed linear fit, this returns the R^2 value of the regression.*/
static GEN
rsquared(GEN X, GEN y, GEN fit)
{
  pari_sp av = avma;
  long n = lg(y) - 1, i;/*Number of observations*/
  GEN predicted = RgV_RgM_mul(shallowtrans(fit), X);/*1xn matrix of the fitted values.*/
  GEN yavg = gen_0;
  for (i = 1; i <= n; i++) yavg = gadd(yavg, gel(y, i));
  yavg = gdivgs(yavg, n);/*Average value of y*/
  GEN sstot=gen_0, ssres=gen_0;
  for (i = 1; i <= n; i++) {
    sstot = gadd(sstot, gsqr(gsub(gel(y, i), yavg)));
    ssres = gadd(ssres, gsqr(gsub(gel(y, i), gel(predicted, i))));
  }
  return gerepileupto(av, gsubsg(1, gdiv(ssres, sstot)));
}

/*Assumes plots/build/fname_plotter.tex points to a file making a plot. This compiles the plot, and views it if WSL=1. If there isn't enough memory, something like lualatex -shell-escape NAME_plotter.tex works.*/
static void
plot_compile(char *fname, int WSL)
{
  pari_sp av = avma;
  char *line = stack_sprintf("(cd ./plots/build && pdflatex --interaction=batchmode -shell-escape %s_plotter.tex)", fname);/*Build*/
  int s = system(line);
  if (s == -1) pari_err(e_MISC, "ERROR EXECUTING COMMAND");
  line = stack_sprintf("mv -f ./plots/build/%s.pdf ./plots/", fname);/*Move the file*/
  s = system(line);
  if (s == -1) pari_err(e_MISC, "ERROR EXECUTING COMMAND");
  if (WSL) {
    line = stack_sprintf("cmd.exe /C start plots/%s.pdf", fname);/*Open the file*/
    s = system(line);
    if (s == -1) pari_err(e_MISC, "ERROR EXECUTING COMMAND");
  }
  set_avma(av);
}

