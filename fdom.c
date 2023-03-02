/*
POSSIBLE FUTURE ADDITIONS:
1) Parallelization of element enumeration along with partial domain computations.
2) Methods of Imbert (see Voight's original paper) to find a minimal presentation in canonical form.
3) Computation of cohomology groups.
*/

/*CHANGE THE LONG DECLARATIONS OUT OF FOR LOOPS*/

/*INCLUSIONS*/

#ifndef PARILIB
#define PARILIB
#include <pari/pari.h>
#endif

#ifndef PARIPRIVLIB
#define PARIPRIVLIB
#include <pari/paripriv.h>
#endif

#ifndef FDOMDECL
#define FDOMDECL
#include "fdom.h"
#endif


/*STATIC DECLARATIONS*/

/*1: SHORT VECTORS IN LATTICES*/
//static GEN quadraticintegernf(GEN nf, GEN A, GEN B, GEN C, long prec);
//static GEN smallvectors_cholesky(GEN Q, GEN C, long maxelts, GEN condition, long prec);
//static GEN smallvectors_nfcondition(GEN A, GEN C, long maxelts, GEN condition, long prec);

/*2: INTERSECTION OF LINES/CIRCLES*/
//static GEN arc_int(GEN c1, GEN c2, GEN tol, long prec);
//static GEN arcseg_int(GEN c, GEN l, GEN tol, long prec);
//static GEN circle_int(GEN c1, GEN c2, GEN tol, long prec);
//static GEN circleline_int(GEN c, GEN l, GEN tol, long prec);
//static GEN line_int(GEN l1, GEN l2, GEN tol, long prec);
//static int onarc(GEN c, GEN p, GEN tol, long prec);
//static int onseg(GEN l, GEN p, GEN tol, long prec);

/*2: DISTANCES*/
//static GEN hdist_ud(GEN z1, GEN z2, long prec);
//static GEN hpolygon_area(GEN circles, GEN vertices, GEN tol, long prec);

/*2: FUNDAMENTAL DOMAIN COMPUTATION*/
//static GEN edgepairing(GEN U, GEN tol, int rboth, long prec);
//static GEN normalizedbasis_shiftpoint(GEN c, GEN r, int initial, long prec);
//static GEN normalizedboundary_append(GEN Ubase, GEN G, GEN mats, GEN id, GEN tol, long prec);
//static GEN normalizedboundary_givenU(GEN Ubase, GEN G, GEN mats, GEN id, GEN data, GEN (*gamtopsl)(GEN, GEN, long), GEN tol, long prec);
//static GEN normalizedboundary_givencircles(GEN G, GEN mats, GEN id, GEN tol, long prec);
//static long normalizedboundary_outside(GEN U, GEN z, GEN tol, long prec);
//static GEN normalizedboundary_sideint(GEN U, GEN c, int start, GEN tol, long prec);
//static GEN psl_roots(GEN M, GEN tol, long prec);

/*2: FUNDAMENTAL DOMAIN OTHER COMPUTATIONS*/
//static GEN minimalcycles(GEN pair);
//static void presentation_updatewords(GEN words, long ind, GEN repl);
//static GEN word_collapse(GEN word);
//static GEN word_substitute(GEN word, long ind, GEN repl, GEN invrepl);

/*2: GEOMETRIC HELPER METHODS*/
//static GEN anglediff(GEN ang, GEN bot, GEN tol, long prec);
//static GEN atanoo(GEN x, long prec);
//static int gcmp_strict(void *data, GEN x, GEN y);
//static int geom_check(GEN c);
//static GEN shiftangle(GEN ang, GEN bot, GEN tol, long prec);
static int tolcmp(GEN x, GEN y, GEN tol);
static int tolcmp_sort(void *data, GEN x, GEN y);
static int toleq(GEN x, GEN y, GEN tol);
static int toleq0(GEN x, GEN tol);

/*3: QUATERNION ALGEBRA NON-FDOM METHODS*/
//static GEN algd(GEN A, GEN a);
//static GEN voidalgmul(void *A, GEN x, GEN y);

/*3: FUNDAMENTAL DOMAIN COMPUTATION*/
//static GEN qalg_fdom(GEN Q, GEN p, int dispprogress, GEN C, GEN R, GEN passes, int type, long prec);

/*3: HELPER METHODS*/
//static long algsplitoo(GEN A);
//static GEN qalg_normform_givenbasis(GEN Q, GEN basis);
//static GEN qalg_basis_conj(GEN Q, GEN x);


/*EXTRA*/
//static GEN elementabsmultable(GEN mt, GEN x);
//static GEN elementabsmultable_Fp(GEN mt, GEN x, GEN p);
//static GEN algbasismultable(GEN al, GEN x);
//static GEN algtracebasis(GEN al);
//static GEN elementabsmultable_Z(GEN mt, GEN x);
//static GEN FpM_trace(GEN x, GEN p);
//static GEN ZM_trace(GEN x);

/*SECTION 1: BASE METHODS*/



/*INFINITY */

/*LISTS*/

/*TARGET REMOVAL OF THIS SECTION? MAYBE*/

/*Appends x to v, returning v, and updating vind to vind++. If vind++>vlen, then we double the length of v as well. If this happens, the resulting vector is not suitable for gerepileupto; this must be done at the end (necessary anyway since it's likely we have to call vec_shorten at some point).*/
GEN
veclist_append(GEN v, long *vind, long *vlen, GEN x)
{
  if(*vind==*vlen){/*Need to lengthen!*/
    *vlen=2**vlen;
    v=vec_lengthen(v, *vlen);
  }
  *vind=*vind+1;
  gel(v, *vind)=x;
  return v;
}

/*Appends x to v, returning v, and updating vind to vind++. If vind++>vlen, then we double the length of v as well. Don't forget to call vec_shorten at the end, since some positions are uninitialized.*/
GEN
vecsmalllist_append(GEN v, long *vind, long *vlen, long x)
{
  if(*vind==*vlen){/*Need to lengthen!*/
    *vlen=2**vlen;
    v=vecsmall_lengthen(v, *vlen);
  }
  *vind=*vind+1;
  v[*vind]=x;
  return v;
}



/*SHORT VECTORS IN LATTICES*/


/*SECTION 2: GEOMETRY METHODS*/



/*A line is stored as [a, b, c], representing ax+by=c. We will normalize so that c=1 or 0. It is assumed that at least one of a, b is non-zero
A segment is stored as [a, b, c, x0, x1]. [a, b, c] gives the line, which has start point x0 and ends at x1, which are complex. We do not allow segments going through oo
GEN tol -> The tolerance, which MUST be of type t_REAL.*/


/*BASIC LINE, CIRCLE, AND POINT OPERATIONS*/

/*INTERSECTION OF LINES/CIRCLES*/


/*DISTANCES/AREAS*/

/*FUNDAMENTAL DOMAIN COMPUTATION*/

/*FUNDAMENTAL DOMAIN OTHER COMPUTATIONS*/

/*GEOMETRIC HELPER METHODS*/

/*Returns the default tolerance given the precision.*/
GEN
deftol(long prec)
{
  pari_sp av=avma;
  return gerepileupto(av, shiftr(gtofp(gen_1, prec), BITS_IN_LONG/2*(2-prec)));
}

/*Returns -1 if x<y, 0 if x==y, 1 if x>y (x, y are t_REAL). Accounts for the tolerance, so will deem x==y if they are equal up to tol AND at least one is inexact*/
static int
tolcmp(GEN x, GEN y, GEN tol)
{
  pari_sp av=avma;
  GEN d=gsub(x, y);
  switch(typ(d)){
    case t_FRAC:/*t_FRAC cannot be 0*/
	  return gc_int(av, signe(gel(d, 1)));
	case t_INT:/*Given exactly*/
	  return gc_int(av, signe(d));
	case t_REAL:
	  if(abscmprr(d, tol)<0) return 0;/*|d|<tol*/
	  return gc_int(av, signe(d));
  }
  pari_err_TYPE("Tolerance comparison only valid for type t_INT, t_FRAC, t_REAL", d);
  return 0;/*So that there is no warning*/
}

/*Data points to tol. Used to sort/search a list with tolerance.*/
static int
tolcmp_sort(void *data, GEN x, GEN y){return tolcmp(x, y, *(GEN*)data);}

/*Returns 1 if x==y up to tolerance tol. If x and y are both t_INT/t_FRAC, will only return 1 if they are exactly equal.*/
static int
toleq(GEN x, GEN y, GEN tol)
{
  pari_sp av=avma;
  GEN d=gsub(x, y);
  return gc_int(av, toleq0(d, tol));/*Just compare d with 0.*/
}

/*If x is of type t_INT or t_FRAC, returns 1 iff x==0. Otherwise, x must be of type t_REAL or t_COMPLEX, and returns 1 iff x=0 up to tolerance tol.*/
static int
toleq0(GEN x, GEN tol)
{
  switch(typ(x)){
    case t_FRAC:/*t_FRAC cannot be 0*/
	  return 0;
	case t_INT:/*Given exactly*/
	  return !signe(x);
	case t_REAL:
	  if(abscmprr(x, tol)<0) return 1;/*|x|<tol*/
	  return 0;
	case t_COMPLEX:;
	  long i;
	  for(i=1;i<=2;i++){
		switch(typ(gel(x, i))){
		  case t_FRAC:/*Fraction component, cannot be 0*/
		    return 0;
		  case t_INT:
		    if(signe(gel(x, i))) return 0;
			break;
		  case t_REAL:
		    if(abscmprr(gel(x, i), tol)>=0) return 0;/*Too large*/
			break;
		  default:/*Illegal input*/
		    pari_err_TYPE("Tolerance equality only valid for type t_INT, t_FRAC, t_REAL, t_COMPLEX", x);
		}
	  }
	  return 1;/*We passed*/
  }
  pari_err_TYPE("Tolerance equality only valid for type t_INT, t_FRAC, t_REAL, t_COMPLEX", x);
  return 0;/*So that there is no warning*/
}





