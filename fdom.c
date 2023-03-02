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

/*1: INFINITY */
static GEN divoo(GEN a, GEN b);

/*1: LISTS*/
//static long gen_search_old(GEN T, GEN x, long flag, void *data, int (*cmp)(void*, GEN, GEN));

/*1: SHORT VECTORS IN LATTICES*/
//static GEN quadraticintegernf(GEN nf, GEN A, GEN B, GEN C, long prec);
//static GEN smallvectors_cholesky(GEN Q, GEN C, long maxelts, GEN condition, long prec);
//static GEN smallvectors_nfcondition(GEN A, GEN C, long maxelts, GEN condition, long prec);

/*2: BASIC LINE, CIRCLE, AND POINT OPERATIONS*/
//static GEN arc_init(GEN c, GEN p1, GEN p2, int dir, long prec);
//static GEN arc_midpoint(GEN c, GEN p1, GEN p2, GEN tol, long prec);
//static GEN circle_angle(GEN c1, GEN c2, GEN p, GEN tol, long prec);
//static GEN circle_fromcp(GEN cent, GEN p, long prec);
//static GEN circle_fromppp(GEN p1, GEN p2, GEN p3, GEN tol, long prec);
//static GEN circle_tangentslope(GEN c, GEN p, long prec);
//static GEN line_fromsp(GEN s, GEN p);
//static GEN line_frompp(GEN p1, GEN p2, GEN tol, long prec);
//static GEN midpoint(GEN p1, GEN p2);
//static GEN mobius(GEN M, GEN c, GEN tol, long prec);
//static GEN mobius_arcseg(GEN M, GEN c, int isarc, GEN tol, long prec);
//static GEN mobius_circle(GEN M, GEN c, GEN tol, long prec);
//static GEN mobius_line(GEN M, GEN l, GEN tol, long prec);
//static GEN perpbis(GEN p1, GEN p2, GEN tol, long prec);
//static GEN radialangle(GEN c, GEN p, GEN tol, long prec);
//static GEN slope(GEN p1, GEN p2, GEN tol, long prec);

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
//static int tolcmp(GEN x, GEN y, GEN tol, long prec);
//static int tolcmp_sort(void *data, GEN x, GEN y);
//static int toleq(GEN x, GEN y, GEN tol, long prec);
//static int toleq0(GEN x, GEN tol, long prec);

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



/*Divides a and b, and allows for oo and division by 0. Returns oo for 0/0.*/
static GEN 
divoo(GEN a, GEN b)
{
  if(gequal0(b)){/*b=0*/
    if(gcmpgs(a, 0)>=0) return mkoo();
    return mkmoo();
  }
  if(typ(a)==t_INFINITY){/*a=+/-oo*/
    if(gsigne(a)==gsigne(b)) return mkoo();
    return mkmoo();
  }
  if(typ(b)==t_INFINITY) return gen_0;
  return gdiv(a, b);
}



/*LISTS*/

/*TARGET REMOVAL OF THIS SECTION? MAYBE*/

/*Appends x to v, returning v, and updating vind to vind++. If vind++>vlen, then we double the length of v as well. If this happens, the resulting vector is not suitable for gerepileupto; this must be done at the end (necessary anyway since it's likely we have to call vec_shorten at some point).*/
GEN
veclist_append(GEN v, long *vind, long *vlen, GEN x){
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
vecsmalllist_append(GEN v, long *vind, long *vlen, long x){
  if(*vind==*vlen){/*Need to lengthen!*/
    *vlen=2**vlen;
    v=vecsmall_lengthen(v, *vlen);
  }
  *vind=*vind+1;
  v[*vind]=x;
  return v;
}

