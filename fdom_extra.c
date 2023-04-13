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

/*SECTION 3: FINCKE POHST PRUNING TESTING*/

/*3: SUPPORTING METHODS TO FINCKE POHST*/
static int check_bound(GEN B, GEN xk, GEN yk, GEN zk, GEN vk);
static GEN clonefill(GEN S, long s, long t);
static int mplessthan(GEN x, GEN y);
static int mpgreaterthan(GEN x, GEN y);
static GEN norm_aux(GEN xk, GEN yk, GEN zk, GEN vk);

/*3: MAIN METHODS*/
static GEN smallvectors_prune(GEN q, GEN C);


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
      t=timer_get(&T);
      set_avma(av2);
    }
    gel(avgtimes, i) = rdivss(t, 1000*tries, prec);
    if (i%10 == 0) pari_printf("%d test cases done.\n", i);
  }
  return gerepilecopy(av, avgtimes);
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





/*SECTION 3: FINCKE POHST PRUNING TESTING*/

/*3: SUPPORTING METHODS TO FINCKE POHST*/

/* increment ZV x, by incrementing cell of index k. Initial value x0[k] was
 * chosen to minimize qf(x) for given x0[1..k-1] and x0[k+1,..] = 0
 * The last nonzero entry must be positive and goes through x0[k]+1,2,3,...
 * Others entries go through: x0[k]+1,-1,2,-2,...*/
INLINE void
step(GEN x, GEN y, GEN inc, long k)
{
  if (!signe(gel(y, k))) /* x[k+1..] = 0 */
    gel(x,k) = addiu(gel(x, k), 1); /* leading coeff > 0 */
  else
  {
    long i = inc[k];
    gel(x,k) = addis(gel(x,k), i),
    inc[k] = (i > 0)? -1-i: 1-i;
  }
}

/* x a t_INT, y  t_INT or t_REAL */
INLINE GEN
mulimp(GEN x, GEN y)
{
  if (typ(y) == t_INT) return mulii(x, y);
  return signe(x) ? mulir(x, y): gen_0;
}

/* x + y*z, x,z two mp's, y a t_INT */
INLINE GEN
addmulimp(GEN x, GEN y, GEN z)
{
  if (!signe(y)) return x;
  if (typ(z) == t_INT) return mpadd(x, mulii(y, z));
  return mpadd(x, mulir(y, z));
}

/* yk + vk * (xk + zk)^2 < B + epsilon */
static int
check_bound(GEN B, GEN xk, GEN yk, GEN zk, GEN vk)
{
  pari_sp av = avma;
  int f = mpgreaterthan(norm_aux(xk,yk,zk,vk), B);
  return gc_bool(av, !f);
}

static GEN
clonefill(GEN S, long s, long t)
{ /* initialize to dummy values */
  GEN T = S, dummy = cgetg(1, t_STR);
  long i;
  for (i = s+1; i <= t; i++) gel(S,i) = dummy;
  S = gclone(S); if (isclone(T)) gunclone(T);
  return S;
}

/* 1 if we are "sure" that x < y, up to few rounding errors, i.e.
 * x < y - epsilon. More precisely :
 * - sign(x - y) < 0
 * - lgprec(x-y) > 3 || expo(x - y) - expo(x) > -24 */
static int
mplessthan(GEN x, GEN y)
{
  pari_sp av = avma;
  GEN z = mpsub(x, y);
  set_avma(av);
  if (typ(z) == t_INT) return (signe(z) < 0);
  if (signe(z) >= 0) return 0;
  if (realprec(z) > LOWDEFAULTPREC) return 1;
  return ( expo(z) - mpexpo(x) > -24 );
}

/* 1 if we are "sure" that x > y, up to few rounding errors, i.e.
 * x > y + epsilon */
static int
mpgreaterthan(GEN x, GEN y)
{
  pari_sp av = avma;
  GEN z = mpsub(x, y);
  set_avma(av);
  if (typ(z) == t_INT) return (signe(z) > 0);
  if (signe(z) <= 0) return 0;
  if (realprec(z) > LOWDEFAULTPREC) return 1;
  return ( expo(z) - mpexpo(x) > -24 );
}

/* yk + vk * (xk + zk)^2 */
static GEN
norm_aux(GEN xk, GEN yk, GEN zk, GEN vk)
{
  GEN t = mpadd(xk, zk);
  if (typ(t) == t_INT) { /* probably gen_0, avoid loss of accuracy */
    yk = addmulimp(yk, sqri(t), vk);
  } else {
    yk = mpadd(yk, mpmul(sqrr(t), vk));
  }
  return yk;
}


/*3: MAIN METHODS*/



/* Increment the ZV x at step k. If x[k+1]=...=x[n]=0, then we start at x[k]=0 and add 1 every time. Otherwise, we start at the value of x[k] that minimizes (x[k]+partnorms[k])^2 (closest integer to -partnorms[k]), and then do x[k]+1, x[k]-1, x[k]+2, x[k]-2, .... inc is a Vecsmall that tracks the variation from x[k] that we have gone. If we are at x[k]+u, then inc=-2*u (u>0), and if we are at x[k]-u, then inc=2*u+1 (u<=0).*/
INLINE void
sv_step(GEN x, GEN partnorms, GEN inc, long k)
{
  if (!signe(gel(partnorms, k + 1))) gel(x, k) = addiu(gel(x, k), 1); /*partnorms[k]=0 <==> x[k+1]=...=x[n]=0*/
  else {
    long i = inc[k];
    gel(x, k) = addis(gel(x, k), i),
    inc[k] = (i > 0)? -1 - i: 1 - i;
  }
}


/*
The quadratic form is given by sum_{i=1}^n q_{i,i}(x_i+sum_{j=i+1}^n q_{i,j}*xj)^2
prune should be a vector of length n, entries real/integer in (0, 1], used to prune the search space.
*/
static GEN
sv_prune(GEN q, GEN C, GEN prune)
{
  long N, i;
  GEN bounds = cgetg_copy(prune, &N);/*N=lg(q) as well.*/
  for (i = 1; i < N; i++) gel(bounds, i) = gmul(C, gel(prune, i));/*The pruning bounds.*/
  pari_sp av = avma;
  long maxv = 2000, found = 0;/*Maximum found vectors, may grow.*/
  GEN v = cgetg(maxv + 1, t_VEC);/*To store them*/
  long n = N - 1;
  /*Page 22 of https://iacr.org/archive/eurocrypt2010/66320257/66320257.pdf claims to give a 40% speedup by carefull tracking the partial sums of q_{i,j}x_j. The matrix partsums accomplishes this. If j>i, then partsums[i, j]=sum_{l>=j}q_{i, l}*x_l IF it has been updated. The key is we don't update everything at every step. Thus the ith term in q(x) is q_ii(x_i+partsums[i, i+1])^2*/
  GEN partsums = zeromatcopy(n, N);/*Only need n-1 rows and n columsn, but use n/N for convenience: the last row and column is always zero.*/
  GEN r = cgetg(N, t_VECSMALL);/*For a fixed i, partsums[i, j] is correct from j=r[i] to j=n. If j>n, it is incorrect everywhere, and clearly r[i]>=i+1 for all i. Initially, since x=0, everything is correct.*/
  for (i = 1; i < N; i++) r[i] = i + 1;
  GEN x = zerovec(n);/*The vector we are making.*/
  GEN partnorms = zerovec(N);/*partnorms[i]=sum_{j=i}^n q_{j, j}(x_j+partsums[j, j+1])^2*/
  GEN inc = const_vecsmall(n, 1);/*Tracks variation from the starting point for x[k].*/
  long k = n, j;
  int down = 1;/*down = 1 means the last step was down, 0 means up or k=1 and the same.*/
  for (;;) {/*Work down the tree until we hit k=1.*/
    if (gc_needed(av, 2)) {/*Garbage day!*/
        v = clonefill(v, found, maxv);/*Clone to the heap the already found vectors.*/
        gerepileall(av, 5, &partsums, &r, &x, &partnorms, &inc);
    }
	if (down) {/*We have just gone down, so must initialize inc[k], x[k], partsums[k,], and r[k]*/
	  inc[k] = 1;
	  for (j = r[k] - 1; j > k; j--) {/*Update the row from right to left*/
	    gcoeff(partsums, k, j) = addmulimp(gcoeff(partsums, k, j + 1), gel(x, j), gcoeff(q, k, j));
	  }
	  r[k] = k + 1;/*Whole row is correct.*/
	  gel(x, k) = mpround(mpneg(gcoeff(partsums, k, k + 1)));/*Smallest possible value for x_k+sum_{j=k+1}^n q_{i, j}x_j*/
	}
	else sv_step(x, partnorms, inc, k);/*We went up, so increment x[k].*/
	GEN nextnorm = norm_aux(gel(x, k), gel(partnorms, k + 1), gcoeff(partsums, k, k + 1), gcoeff(q, k, k));/*Next norm*/
	if (mpgreaterthan(nextnorm, gel(bounds, k))) {
	  sv_step(x, partnorms, inc, k);/*If we are too big, then increasing x once MIGHT still give a valid x. Increasing it twice 100% won't.*/
	  nextnorm = norm_aux(gel(x, k), gel(partnorms, k + 1), gcoeff(partsums, k, k + 1), gcoeff(q, k, k));/*Next norm*/
	  if (mpgreaterthan(nextnorm, gel(bounds, k))) {/*We have run out of steam at this level, go on to the next one.*/
	    long next = k + 1;
	    if (next > n) break;/*Done!*/
	    for (j = k; j > 0; j--) {/*Update r[j]*/
	      if (r[j] <= next) r[j] = next + 1;
	      else break;/*As soon as r[j]>k, it is true for the rest of the j's to 1.*/
	    }
	    k = next;
	    down = 0;/*Going up!*/
	    continue;
	  }
	}
	gel(partnorms, k) = nextnorm;/*Still valid!*/
	if (k > 1) { k--; down = 1; continue; }/*Go down the tree.*/
	/*Now, there is a solution to add.*/
	down = 0;/*Want to stay here and just increment.*/
	if (!signe(gel(partnorms, 1))) continue;/*Exclude 0.*/
	found++;
	gel(v, found) = leafcopy(x);
	if (found != maxv) continue;/*Still room.*/
	maxv <<= 1;/*Double the size*/
    GEN vnew = clonefill(vec_lengthen(v, maxv), found, maxv);
    if (isclone(v)) gunclone(v);
    v = vnew;
  }
  setlg(v, found + 1);
  settyp(v, t_MAT);
  if (isclone(v)) {GEN p1 = v; v = gcopy(v); gunclone(p1); }
  return v;
}





/* q is the Cholesky decomposition of a quadratic form
 * Enumerate vectors whose norm is less than BORNE (Algo 2.5.7),
 * minimal vectors if BORNE = NULL (implies check = NULL).
 * If (check != NULL) consider only vectors passing the check, and assumes
 *   we only want the smallest possible vectors */
static GEN
smallvectors_prune(GEN q, GEN C)
{
  long N = lg(q), n = N - 1, i, j, k, s, stockmax;
  pari_sp av, av1;
  GEN inc, S, x, y, z, v, p1, alpha;
  GEN norme1, normax1, borne1, borne2;

  alpha = dbltor(0.95);
  normax1 = gen_0;

  v = cgetg(N, t_VEC);
  inc = const_vecsmall(n, 1);

  av = avma;
  stockmax = 2000;
  S = cgetg(stockmax + 1, t_VEC);
  x = cgetg(N, t_COL);
  y = cgetg(N, t_COL);
  z = cgetg(N, t_COL);
  for (i = 1; i< N; i++) {
    gel(v, i) = gcoeff(q, i, i);
    gel(x, i) = gel(y, i) = gel(z, i) = gen_0;
  }
  borne1 = C;
  if (gsigne(borne1) <= 0) retmkvec3(gen_0, gen_0, cgetg(1, t_MAT));
  if (typ(borne1) != t_REAL) {
    long prec = nbits2prec(gexpo(borne1) + 10);
    borne1 = gtofp(borne1, maxss(prec, DEFAULTPREC));
  }
  borne2 = mulrr(borne1, alpha);
  s = 0; k = n;
  for(;; step(x, y, inc, k)) {/* main */
	 /* x (supposedly) small vector, ZV.
     * For all t >= k, we have
     *   z[t] = sum_{j > t} q[t,j] * x[j]
     *   y[t] = sum_{i > t} q[i,i] * (x[i] + z[i])^2
     *        = 0 <=> x[i]=0 for all i>t */
    do {
      int skip = 0;
      if (k > 1) {
        long l = k-1;
        av1 = avma;
        p1 = mulimp(gel(x, k), gcoeff(q, l, k));
        for (j = k + 1; j < N; j++) p1 = addmulimp(p1, gel(x, j), gcoeff(q, l, j));
        gel(z, l) = gerepileuptoleaf(av1, p1);

        av1 = avma;
        p1 = norm_aux(gel(x, k), gel(y, k), gel(z, k), gel(v, k));
        gel(y, l) = gerepileuptoleaf(av1, p1);
        /* skip the [x_1,...,x_skipfirst,0,...,0] */
        if (mplessthan(borne1, gel(y, l))) skip = 1;
        else /* initial value, minimizing (x[l] + z[l])^2, hence qf(x) for
                the given x[1..l-1] */
          gel(x, l) = mpround(mpneg(gel(z, l)));
        k = l;
      }
      for(;; step(x,y,inc,k)) { /* at most 2n loops */
        if (!skip) {
          if (check_bound(borne1, gel(x, k), gel(y, k), gel(z, k), gel(v, k))) break;
          step(x, y, inc, k);
          if (check_bound(borne1, gel(x, k), gel(y, k), gel(z, k), gel(v, k))) break;
        }
        skip = 0; inc[k] = 1;
        if (++k > n) goto END;
      }

      if (gc_needed(av,2)) {
        S = clonefill(S, s, stockmax);
        gerepileall(av, 6, &x, &y, &z, &normax1, &borne1, &borne2);
      }
    }
    while (k > 1);
	output(x);
    if (!signe(gel(x, 1)) && !signe(gel(y, 1))) continue; /* exclude 0 */

    av1 = avma;
    norme1 = norm_aux(gel(x, 1), gel(y, 1), gel(z, 1), gel(v, 1));
    if (mpgreaterthan(norme1, borne1)) { set_avma(av1); continue; /* main */ }

    norme1 = gerepileuptoleaf(av1, norme1);
    if (mpcmp(norme1, normax1) > 0) normax1 = norme1;
    if (++s > stockmax) continue; /* too many vectors: no longer remember */
    gel(S, s) = leafcopy(x);
    if (s != stockmax) continue; /* still room, get next vector */
    stockmax <<= 1;
    GEN Snew = clonefill(vec_lengthen(S,stockmax), s, stockmax);
    if (isclone(S)) gunclone(S);
    S = Snew;
  }
END:
  if (s < stockmax) stockmax = s;
  setlg(S, stockmax + 1);
  settyp(S, t_MAT);
  if (isclone(S)) { p1 = S; S = gcopy(S); gunclone(p1); }
  return mkvec3(utoi(s << 1), borne1, S);
}

GEN
qfminim_prune(GEN M, GEN C, int prunetype, long prec)
{
  GEN res = fp_prune(M, C, prunetype, prec);
  if (!res) pari_err_PREC("qfminim");
  return res;
}


/*Solve q(x)=x~*M*x <= C, M is positive definite with real entries. We only require M to be upper triangular. We use pruning. This is mostly a copy of fincke_pohst from bibli1.c. For some reason, this is 2-3 times as fast WITHOUT changing anything else, i.e. the smallvector method.*/
GEN
fp_prune(GEN M, GEN C, int prunetype, long PREC)
{
  pari_sp av = avma;
  long prec = PREC;
  long lM = lg(M);
  if (lM == 1) retmkvec3(gen_0, gen_0, cgetg(1, t_MAT));
  GEN U = lllfp(M, 0.75, LLL_GRAM | LLL_IM);/*LLL reduce our input matrix*/
  if (lg(U) != lM) return gc_NULL(av);
  GEN R = qf_apply_RgM(M, U);/*U~*M*U*/
  long i = gprecision(R), j;
  if (i) prec = i;
  else {
    prec = DEFAULTPREC + nbits2extraprec(gexpo(R));
    if (prec < PREC) prec = PREC;
  }
  R = qfgaussred_positive(R);
  if (!R) return gc_NULL(av);
  for (i = 1; i < lM; i++) {
    GEN s = gsqrt(gcoeff(R, i, i), prec);
    gcoeff(R, i, i) = s;
    for (j = i + 1; j < lM; j++) gcoeff(R, i, j) = gmul(s, gcoeff(R, i, j));
  }
  /* now R~ * R = a in LLL basis */
  GEN Rinv = RgM_inv_upper(R);
  if (!Rinv) return gc_NULL(av);
  GEN Rinvtrans = shallowtrans(Rinv);
  GEN V = lll(Rinvtrans);
  if (lg(V) != lM) return gc_NULL(av);

  Rinvtrans = RgM_mul(Rinvtrans, V);
  V = ZM_inv(shallowtrans(V), NULL);
  R = RgM_mul(R, V);
  U = ZM_mul(U, V);

  GEN Vnorm = cgetg(lM, t_VEC);
  for (j = 1; j < lM; j++) gel(Vnorm,j) = gnorml2(gel(Rinvtrans, j));
  GEN Rperm = cgetg(lM, t_MAT);
  GEN Uperm = cgetg(lM, t_MAT), perm = indexsort(Vnorm);
  for (i = 1; i < lM; i++) { Uperm[lM - i] = U[perm[i]]; Rperm[lM - i] = R[perm[i]]; }
  U = Uperm;
  R = Rperm;
  GEN res = NULL;
  /*
  pari_CATCH(e_PREC) { }
  pari_TRY {
    GEN q = gaussred_from_QR(R, gprecision(Vnorm));
    if (q) res = smallvectors_prune(q, C);
  } pari_ENDCATCH;
  if (!res) return gc_NULL(av);
  
  GEN z = cgetg(4,t_VEC);
  gel(z, 1) = gcopy(gel(res, 1));
  gel(z, 2) = gcopy(gel(res, 2));
  gel(z, 3) = ZM_mul(U, gel(res,3));
  return gerepileupto(av, z);
  */
  GEN prune;
  if (!prunetype) prune = const_vec(lM - 1, gen_1);/*No funny business, just normal Fincke-Pohst.*/
  else {/*Linear pruning*/
	GEN con = dbltor(1.05);
	prune = cgetg(lM, t_VEC);
	long n = lM - 1;
	for (i = 1; i < lM; i++) gel(prune, i) = gmin_shallow(gen_1, divrs(mulrs(con, lM - i), n));
  }
  pari_CATCH(e_PREC) { }
  pari_TRY {
    GEN q = gaussred_from_QR(R, gprecision(Vnorm));
    if (q) res = sv_prune(q, C, prune);
  } pari_ENDCATCH;
  if (!res) return gc_NULL(av);
  return gerepileupto(av, ZM_mul(U, res));
}












