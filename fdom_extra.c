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
static GEN tune_time(GEN X, GEN Cset, long mintesttime, long prec);

/*2: REGRESSIONS*/
static GEN OLS(GEN X, GEN y, int retrsqr);
static GEN rsquared(GEN X, GEN y, GEN fit);

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

/*2: REGRESSIONS*/

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

