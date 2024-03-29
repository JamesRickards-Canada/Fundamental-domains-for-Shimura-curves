/*These methods are useful with the fdom.c package, but not essential.*/

/*INCLUSIONS*/

#include <pari.h>
#include "fdom.h"

/*STATIC DECLARATIONS*/

/*SECTION 1: VISUALIZATION*/
/*SECTION 2: TESTING AND TUNING*/
static int alg_in_centre(GEN A, GEN g);
static GEN afuchmulvec(GEN X, GEN G, GEN L);
/*SECTION 3: EXTRA METHODS OVER Q*/

/*MAIN BODY*/

/*SECTION 1: VISUALIZATION*/

/*Prints the fundamental domain to a latex file. model=0 means Klein, 1 is unit disc, and 2 is upper half plane.*/
void
afuchfdom_latex(GEN X, char *filename, int model, int boundcircle, int compile, int open)
{
  pari_sp av = avma;
  if (model == 2) pari_err(e_MISC, "Upper half plane not yet supported");
  GEN tol = gdat_get_tol(afuch_get_gdat(X));
  long prec = realprec(tol);
  GEN U = afuch_get_fdom(X);
  if (gequal0(U)) pari_err(e_MISC, "Please initialize the fundamental domain first with X = afuchmakefdom(X).");
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
  else if (model == 1) {/*Klein model*/
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

/*Writes the fundamental domain corresponding to U to fdoms/filename/dat, to be used with the Python program fdomviewer. We output to the unit disc model.*/
void
afuchfdom_python(GEN X, char *filename)
{
  pari_sp av = avma;
  if (!pari_is_dir("fdoms")) {/*Checking the directory*/
    int s = system("mkdir -p fdoms");
    if (s == -1) pari_err(e_MISC, "ERROR CREATING DIRECTORY fdoms");
  }
  char *fullfile = stack_sprintf("fdoms/%s.dat", filename);
  FILE *f = fopen(fullfile, "w");/*Now we have created the output file f.*/
  GEN U = afuch_get_fdom(X);
  if (gequal0(U)) pari_err(e_MISC, "Please initialize the fundamental domain first with X = afuchmakefdom(X).");
  GEN pair = normbound_get_spair(U);
  pari_fprintf(f, "%d", pair[1]);
  long i, lp = lg(pair);
  for (i = 2; i < lp; i++) pari_fprintf(f, " %d", pair[i]);/*Print side pairing.*/
  pari_fprintf(f, "\n");
  GEN arcs = normbound_get_sides(U);
  GEN verts = normbound_get_vcors(U);
  GEN tol = gdat_get_tol(afuch_get_gdat(X));
  long prec = realprec(tol);
  GEN radtodeg = divsr(180, mppi(prec));
  for (i = 1; i < lp; i++) {
    GEN arc = gel(arcs, i), v1;
    if (i == 1) v1 = gel(verts, lp - 1);
    else v1 = gel(verts, i - 1);
    GEN v2 = gel(verts, i);/*The two vertices in the Klein model.*/
    v1 = klein_to_disc(v1, tol);
    v2 = klein_to_disc(v2, tol);
    GEN xc = gel(arc, 1);
    GEN yc = gel(arc, 2);/*The coords of the centre of the isometric circle.*/
    GEN centre = mkcomplex(xc, yc);
    GEN r = gel(arc, 3);
    GEN ang1 = mulrr(garg(gsub(v1, centre), prec), radtodeg);
    GEN ang2 = mulrr(garg(gsub(v2, centre), prec), radtodeg);
    pari_fprintf(f, "%P.20f %P.20f %P.20f %P.20f %P.20f\n", xc, yc, r, ang1, ang2);
  }
  fclose(f);
  set_avma(av);
}

/*Writes the geodesic given by g (or the output of afuchgeodesic(X, g)) to fdoms/filename.dat, to be used with the Python program fdomviewer.*/
void
afuchgeodesic_python(GEN X, GEN g, char *filename)
{
  pari_sp av = avma;
  if (!pari_is_dir("fdoms")) {/*Checking the directory*/
    int s = system("mkdir -p fdoms");
    if (s == -1) pari_err(e_MISC, "ERROR CREATING DIRECTORY fdoms");
  }
  char *fullfile = stack_sprintf("fdoms/%s.dat", filename);
  FILE *f = fopen(fullfile, "w");/*Now we have created the output file f.*/
  GEN tol = gdat_get_tol(afuch_get_gdat(X));
  long prec = realprec(tol);
  if (typ(g) == t_COL) g = afuchgeodesic(X, g);
  GEN radtodeg = divsr(180, mppi(prec));
  long i, lgeo = lg(g);
  for (i = 1; i < lgeo; i++) {
    GEN dat = gel(g, i);/*The geodesic segment.*/
    GEN v1 = klein_to_disc(gel(dat, 4), tol);
    GEN v2 = klein_to_disc(gel(dat, 5), tol);
    GEN eqn = gel(dat, 6);
    if (gequal1(gel(eqn, 3))) {/*Arc*/
      GEN xc = gel(eqn, 1);
      GEN yc = gel(eqn, 2);/*The coords of the centre of the circle.*/
      GEN centre = mkcomplex(xc, yc);
      GEN r = sqrtr(subrs(addrr(sqrr(xc), sqrr(yc)), 1));/*sqrt(xc^2+yc^2-1)=r.*/
      GEN ang1 = mulrr(garg(gsub(v1, centre), prec), radtodeg);
      GEN ang2 = mulrr(garg(gsub(v2, centre), prec), radtodeg);
      GEN diff = gmodgs(subrr(ang2, ang1), 360);
      long dir;
      if (gcmpgs(diff, 180) > 0) { dir = -1; GEN temp = ang2; ang2 = ang1; ang1 = temp; }
      else dir = 1;
      ang1 = gmodgs(ang1, 360);
      ang2 = gmodgs(ang2, 360);
      if (cmprr(ang2, ang1) < 0) ang2 = addrs(ang2, 360);
      pari_fprintf(f, "0 %P.20f %P.20f %P.20f %P.20f %P.20f %d\n", xc, yc, r, ang1, ang2, dir);
    }
    else {/*Segment through the origin.*/
      pari_fprintf(f, "1 %P.20f %P.20f %P.20f %P.20f\n", real_i(v1), imag_i(v1), real_i(v2), imag_i(v2));
    }
  }
  fclose(f);
  set_avma(av);
}

/*Launches the Python viewer. Only works with WSL where Python is installed in Windows.*/
void
fdomviewer(char *input)
{
  pari_sp av = avma;
  char *command;
  command = stack_sprintf("cmd.exe /C start py fdomviewer.py %s", input);
  int s = system(command);
  if (s == -1) pari_err(e_MISC, "ERROR EXECUTING COMMAND");
  set_avma(av);
}


/*SECTION 2: TESTING AND TUNING*/

/*Runs a series of checks on X with the fundamental domain and presentation initialized. Returns 0 if all passed, and something non-zero else. These return codes are:
  1: signature area formula does not match computed area;
  2: presentation has wrong number of generators;
  3: presentation has wrong number of relations;
  4: one of the relations fails;
  5: one of the side pairing element relations fails;
  6: afuchfdomword fails on a random element (15 tested).
  7: there are too many / few elliptic elements, or some of their orders are wrong.
*/
long
afuchcheck(GEN X)
{
  pari_sp av = avma;
  GEN pres = afuch_get_presentation(X);
  if (gequal0(pres)) pari_err(e_MISC, "Please initialize the fundamental domain and presentation first with X = afuchmakefdom(X).");
  long prec = afuch_get_prec(X);
  GEN sig = afuch_get_signature(X);
  long genus = itos(gel(sig, 1)), npar = itos(gel(sig, 3));
  GEN asig = stoi(((genus - 1) << 1) + npar);/*2g-2 + e_oo*/
  long i, nell = lg(gel(sig, 2)) - 1;
  for (i = 1; i <= nell; i++) asig = gadd(asig, mkfracss(gel(sig, 2)[i] - 1, gel(sig, 2)[i]));/*Add 1-1/e_i*/
  asig = gmul(asig, Pi2n(1, prec));
  GEN aX = normbound_get_area(afuch_get_fdom(X));
  GEN diff = gabs(gsub(asig, aX), prec);
  if (gcmp(diff, real2n(-16, prec)) > 0) return gc_long(av, 1);/*Area not within 2^-16 of signature formula.*/
  
  long ngens = (genus << 1) + nell + npar;/*Generators corresponding to signature*/
  if (nell || npar) ngens--;/*Can eliminate one if there is at least one elliptic or parabolic cycle.*/
  GEN gens = gel(pres, 1), rels = gel(pres, 2);/*Minimal set of generators and the relations.*/
  if (ngens + 1 != lg(gens)) return gc_long(av, 2);/*Wrong number of generators.*/
  
  long nrel;
  if (nell) nrel = nell;/*If there are elliptic cycles, then their count is the number of relations.*/
  else {
    if (npar) nrel = 0;/*No elliptic cycles, but parabolic cycles. No relations.*/
    else nrel = 1;/*No elliptic cycles or parabolic cycles. One relation.*/
  }
  if (nrel + 1 != lg(rels)) return gc_long(av, 3);/*Wrong number of relations.*/
  
  GEN A = afuch_get_alg(X), O = afuch_get_O(X);
  for (i = 1; i <= nrel; i++) {/*Check the relations actually work.*/
    GEN g = afuchmulvec(X, gens, gel(rels, i));
    if (!gequal1(O)) g = QM_QC_mul(O, g);/*Move to base field if required.*/
    if (!alg_in_centre(A, g)) return gc_long(av, 4);/*Relation fails.*/
  }
  
  GEN fdomelts = afuchelts(X);/*Get the side pairing elements in terms of A, not the order.*/
  long nfdomelts = lg(fdomelts) - 1;
  GEN fdomwords = gel(pres, 3);/*fdomelts as words in gens*/
  for (i = 1; i <= nfdomelts; i++) {/*Check that the claimed expressions for these elements are valid.*/
    GEN g = afuchmulvec(X, gens, gel(fdomwords, i));
    if (!gequal1(O)) g = QM_QC_mul(O, g);/*Move to base field if required.*/
    g = algmul(A, g, alginv(A, gel(fdomelts, i)));
    if (!alg_in_centre(A, g)) return gc_long(av, 5);/*Expression for side pairing element fails.*/
  }
  
  long wordtries;
  for (wordtries = 1; wordtries <= 15; wordtries++) {/*Try a random element.*/
    long lword = 6 + random_Fl(15);
    GEN tomul = cgetg(lword, t_VECSMALL);
    for (i = 1; i < lword; i++) {
      tomul[i] = random_Fl(nfdomelts) + 1;
      if (random_Fl(2)) tomul[i] = -tomul[i];/*Invert some elements.*/
    }
    GEN g = algmulvec(A, fdomelts, tomul);
    GEN word = afuchword(X, g);
    GEN gword = afuchmulvec(X, gens, word);
    if (!gequal1(O)) gword = QM_QC_mul(O, gword);/*Move to base field if required.*/
    GEN diff = algmul(A, g, alginv(A, gword));
    if (!alg_in_centre(A, diff)) return gc_long(av, 6);/*Expressing a random element as a word failed.*/
  }
  
  GEN elliptic = afuchelliptic(X);/*Make sure we get the base field version of our elements, not in terms of a stored order.*/
  if (lg(elliptic) != (nell + 1)) return gc_long(av, 7);/*Not the same number of orders as elliptic elements.*/
  for (i = 1; i <= nell; i++) {
    GEN g, gbase;
    g = gbase = gel(elliptic, i);
    long j;
    for (j = 1; j < gel(sig, 2)[i]; j++) {
      if (afuchistriv(X, g)) return gc_long(av, 7);/*Elliptic element has smaller order than claimed.*/
      g = algmul(A, g, gbase);
    }
    if (!afuchistriv(X, g)) return gc_long(av, 7);/*Elliptic element has larger order than claimed.*/
  }
  
  return gc_long(av, 0);
}

/*Returns 1 iff g is in the centre of A, i.e. the base number field.*/
static int
alg_in_centre(GEN A, GEN g)
{
  pari_sp av = avma;
  GEN x;
  if (typ(g) == t_COL) x = algbasisto1ijk(A, g);
  else x = algalgto1ijk(A, g);
  long i;
  for (i = 2; i <= 4; i++) if (!gequal0(gel(x, i))) return gc_int(av, 0);/*i,j,k component not 0 menas not in centre.*/
  return gc_int(av, 1);
}

/*Returns G[L[1]]*G[L[2]]*...*G[L[n]], where L is a vecsmall. We use the faster multiplication afforded by X, and assume the G elements are written in terms of the basis O in X. Rather than inverses, we use the conjugates.*/
static GEN
afuchmulvec(GEN X, GEN G, GEN L)
{
  pari_sp av = avma;
  GEN g;
  if (L[1] > 0) g = gel(G, L[1]);
  else g = afuchconj(X, gel(G, -L[1]));
  long lL = lg(L), i;
  for (i = 2; i < lL; i++) {
    GEN g1;
    if (L[i] > 0) g1 = gel(G, L[i]);
    else g1 = afuchconj(X, gel(G, -L[i]));
    g = afuchmul(X, g, g1);
  }
  return gerepileupto(av, g);
}

/*We have pre-stored algebras of degree n in "testing/fdom_Cn.dat". We compute the time to compute the fundamental domains of the algebras, where each test is repeated testsperalg times, with the range of C_n's. We only have tests for 1<=n<=9, as algebras with n>=10 are too big to compute very quickly.*/
GEN
tune_Cn(long n, GEN Cmin, GEN Cmax, long testsperalg, long tests, long prec)
{
  pari_sp av = avma, av1, av2, av3;
  if (n <= 0 || n >= 10) pari_err(e_MISC, "n must be between 1 and 9.");
  GEN dat = gel(gp_readvec_file("testing/fdom_Cn.dat"), n);
  long ldat = lg(gel(dat, 1)), i, j, k;
  if (tests <= 1) tests = 2;
  GEN Cn = Cmin;
  GEN Cnadd = gdivgs(gsub(Cmax, Cmin), tests - 1);
  GEN times = cgetg(tests + 1, t_VECSMALL);
  FILE *f = fopen("testing/Cn_timings.dat", "a");
  pari_fprintf(f, "Testing n=%d, %d trials per algebra, %d values of C between %P.8f and %P.8f\n", n, testsperalg, tests, Cmin, Cmax);
  pari_printf("Testing n=%d, %d trials per algebra, %d values of C between %P.8f and %P.8f\n", n, testsperalg, tests, Cmin, Cmax);
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
        timer_delay(&T);
        GEN X = afuchinit(A, NULL, gen_0, 0, prec);
        gmael(X, afuch_FDOMDAT, fdomdat_BESTC) = C;
        X = afuchmakefdom(X);
        t = t + timer_delay(&T);
        set_avma(av3);
      }
      set_avma(av2);
    }
    times[i] = t;
    pari_printf("The value C_n=%P.8f took %d time.\n", Cn, t);
    pari_fprintf(f, "%P.8f %d\n", Cn, t);
    Cn = gerepileupto(av1, gadd(Cn, Cnadd));
  }
  fclose(f);
  Cn = Cmin;
  GEN Cs = cgetg(tests + 1, t_VEC);
  for (i = 1; i <= tests; i++) {
    gel(Cs, i) = Cn;
    Cn = gadd(Cn, Cnadd);
  }
  return gerepilecopy(av, mkvec2(Cs, times));
}


/*SECTION 3: EXTRA METHODS OVER Q*/

/*Initializes the quaternion algebra over Q with discriminant being the product of the prime factors of D. Can also pass in D as a vector of the prime factors*/
GEN
alginit_Qdisc(GEN D, long prec)
{
  pari_sp av = avma;
  GEN y = pol_x(1);/*Variable y*/
  GEN F = nfinit(y, prec);
  GEN ps;
  if (typ(D) == t_INT) {
    if (signe(D) == -1) D = negi(D);
    ps = gel(Z_factor(D), 1);/*First column.*/
    settyp(ps, t_VEC);
  }
  else if (typ(D) == t_VEC) ps = D;
  else { pari_err_TYPE("D must be an integer or a vector of prime numbers", D); ps = gen_0; }/*Set ps to avoid error.*/
  long lps, i;
  GEN pdecs = cgetg_copy(ps, &lps);
  for (i = 1; i < lps; i++) gel(pdecs, i) = gel(idealprimedec(F, gel(ps, i)), 1);
  GEN finram = const_vec(lps - 1, gen_1);
  GEN infram;
  if (lps % 2) infram = gen_0;/*oo needs to be ramified, even number of primes*/
  else infram = gen_1;
  return gerepileupto(av, alginit(F, mkvec3(gen_2, mkvec2(pdecs, finram), mkvec(infram)), 0, 3));
}

