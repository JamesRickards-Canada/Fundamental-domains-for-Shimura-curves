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

/*We have pre-stored algebras of degree n in "data_in/fdom_Cn.dat". We compute the time to compute the fundamental domains of the algebras, where each test is repeated testsperalg times, with the range of C_n's. We only have tests for 1<=n<=9, as algebras with n>=10 are too big to compute very quickly.*/
GEN
tune_Cn(long n, GEN Cmin, GEN Cmax, long testsperalg, long tests, long prec)
{
  pari_sp av = avma, av1, av2, av3;
  if (n <= 0 || n >= 10) pari_err(e_MISC, "n must be between 1 and 9.");
  GEN dat = gel(gp_readvec_file("data_in/fdom_Cn.dat"), n);
  long ldat = lg(gel(dat, 1)), i, j, k;
  if (tests <= 1) tests = 2;
  GEN Cn = Cmin;
  GEN Cnadd = gdivgs(gsub(Cmax, Cmin), tests - 1);
  GEN times = cgetg(tests + 1, t_VECSMALL);
  FILE *f = fopen("data_in/Cn_timings.dat", "a");
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
	    GEN X = afuchinit(A, NULL, NULL, NULL, 0, prec);
	    gmael3(X, 7, 6, 2) = C;/*This isn't really safe but should be OK for now. If we change where C is stored, this must change.*/
	    afuchfdom(X);
	    t = t + timer_delay(&T);
		obj_free(X);/*Kill it off.*/
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




