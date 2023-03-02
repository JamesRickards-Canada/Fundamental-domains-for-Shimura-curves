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

/*DEFINITIONS*/

/*The length (lg, so technically length+1) of a circle/line and arc/segment, and a normalized boundary*/
#define LINE_LG 4
#define SEG_LG 6

/*STATIC DECLARATIONS*/

/*1: SHORT VECTORS IN LATTICES*/
//static GEN quadraticintegernf(GEN nf, GEN A, GEN B, GEN C, long prec);
//static GEN smallvectors_cholesky(GEN Q, GEN C, long maxelts, GEN condition, long prec);
//static GEN smallvectors_nfcondition(GEN A, GEN C, long maxelts, GEN condition, long prec);

/*2: INTERSECTION OF LINES/CIRCLES*/
static GEN line_int(GEN l1, GEN l2, GEN tol);
static GEN line_line_detNULL(GEN l1, GEN l2, GEN tol);
static int onseg(GEN l, GEN p, GEN tol);
static GEN seg_int(GEN l1, GEN l2, GEN tol);

/*2: DISTANCES*/
static GEN hdist_ud(GEN z1, GEN z2, long prec);
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
  if(*vind==*vlen)
  {/*Need to lengthen!*/
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
  if(*vind==*vlen)
  {/*Need to lengthen!*/
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

/*Given ab=[a, b], gives the action on x. We take x to be an element of unit disc in the Klein model, and the action is given by the corresponding action of [a, b;conj(b), conj(a)] on the disc model (with |a|^2-|b|^2=1). The action is via K(x)=(a^2x+b^2xconj+2ab)/(|a|^2+|b|^2+a*x*bconj+aconj*xconj*b.*/
GEN klein_act(GEN ab, GEN x){
  pari_sp av=avma;
  GEN a=gel(ab, 1), b=gel(ab, 2);
  GEN ac=conj_i(a), bc=conj_i(b), xc=conj_i(x);
  GEN num=gadd(gadd(gmul(gsqr(a), x), gmul(gsqr(b), xc)), gmulsg(2, gmul(a, b)));/*a^2x+b^2xc+2ab*/
  GEN axbbar_real=real_i(gmul(a, gmul(x, bc)));/*real part of a*x*bc */
  GEN anorm=real_i(gmul(a, ac));/*real part of a*ac=|a|^2*/
  GEN den=gsubgs(gmulsg(2, gadd(anorm, axbbar_real)), 1);/*|a|^2+|b|^2+a*x*bc+ac*xc*b=2|a|^2+2*real(a*q*bc)-1*/
  return gerepileupto(av, gdiv(num, den));
}



/*INTERSECTION OF LINES/CIRCLES*/


/*Returns the intersection of two lines. If parallel/concurrent, returns NULL*/
static GEN
line_int(GEN l1, GEN l2, GEN tol)
{
  pari_sp av=avma;
  GEN c1=seg_get_c(l1), c2=seg_get_c(l2);/*Get c values*/
  GEN det=line_line_detNULL(l1, l2, tol);/*ad-bc*/
  if(!det) return gc_NULL(av);/*Lines are parallel*/
  if(!signe(c1))
  {/*ax+by=0*/
    if(gequal0(c2)) return gc_const(av, gen_0);/*Both must pass through 0*/
	return gerepilecopy(av, mkcomplex(gdiv(gneg(seg_get_b(l1)), det), gdiv(seg_get_a(l1), det)));/*-b/det, a/det*/
  }
  /*Next, cx+dy=0*/
  if(!signe(c2)) return gerepilecopy(av, mkcomplex(gdiv(seg_get_b(l2), det), gdiv(gneg(seg_get_a(l2)), det)));/*d/det, -c/det*/
  /*Now ax+by=cx+dy=1*/
  GEN x=gdiv(gsub(seg_get_b(l2), seg_get_b(l1)), det);/*(d-b)/det*/
  GEN y=gdiv(gsub(seg_get_a(l1), seg_get_a(l2)), det);/*(a-c)/det*/
  return gerepilecopy(av, mkcomplex(x, y));
}

/*Given two lines given by a1x+b1y=c1 and a2x+b2y=c2, returns a1b2-a2b1, unless this is within tolerance of 0, when we return NULL. Not stack clean.*/
static GEN
line_line_detNULL(GEN l1, GEN l2, GEN tol){
  pari_sp av=avma;
  GEN d=gsub(gmul(seg_get_a(l1), seg_get_b(l2)), gmul(seg_get_b(l1), seg_get_a(l2)));
  if(toleq0(d, tol)) return gc_NULL(av);/*d=0 up to tolerance*/
  return d;
}

/*p is assumed to be on the line defined by l; this checks if it is actually on the segment l. Returns 0 if not, 1 if p is in the interior, 2 if p is the start point, and 3 if p is the end point (all up to tolerance tol).*/
static int
onseg(GEN l, GEN p, GEN tol)
{
  if(lg(l)==LINE_LG) return 1;/*Lines contain all points*/
  pari_sp av=avma;
  GEN pstart=seg_get_start(l), px=real_i(p);
  int xcmp1=tolcmp(real_i(pstart), px, tol);/*Compare x coeff of starting point and x*/
  if(!xcmp1)
  {/*Same x-value*/
	GEN py=imag_i(p);
	int ycmp1=tolcmp(imag_i(pstart), py, tol);
	if(!ycmp1) return gc_int(av, 2);/*We must be the starting point.*/
	GEN pend=seg_get_end(l);/*At this point, the slope must be oo*/
	int ycmp2=tolcmp(py, imag_i(pend), tol);
	if(!ycmp2) return gc_int(av, 3);/*Since we assume that p is on the line, we must be at endpoint 2 now.*/
	if(ycmp1==ycmp2) return gc_int(av, 1);/*Must be between*/
	return gc_int(av, 0);/*Above or below*/
  }/*Now we have distinct x-values*/
  int xcmp2=tolcmp(px, real_i(seg_get_end(l)), tol);
  if(xcmp2==0) return gc_int(av, 3);/*End point, as the slope cannot be oo as pstartx!=pendx*/
  if(xcmp1==xcmp2) return gc_int(av, 1);
  return gc_int(av, 0);
}

/*Returns the intersection of two segments. If parallel/concurrent/do not intersect, returns NULL*/
static GEN 
seg_int(GEN l1, GEN l2, GEN tol){
  pari_sp av=avma;
  GEN ipt=line_int(l1, l2, tol);/*Intersect the lines*/
  if(!onseg(l1, ipt, tol)) return gc_NULL(av);
  if(!onseg(l2, ipt, tol)) return gc_NULL(av);
  return ipt;
}


/*DISTANCES/AREAS*/


/*Given the area of a hyperbolic disc, this returns the radius. The formula is area=4*Pi*sinh(R/2)^2, or R=2arcsinh(sqrt(area/4Pi))*/
GEN
hdiscradius(GEN area, long prec)
{
  pari_sp av=avma;
  return gerepileupto(av, gtofp(gmulgs(gasinh(gsqrt(gdiv(area, Pi2n(2, prec)), prec), prec), 2), prec));
}

/*Returns a random point z in the unit disc, uniform inside the ball of radius R. See page 19 of Page (before section 2.5).*/
GEN
hdiscrandom(GEN R, long prec)
{
  pari_sp av=avma;
  GEN arg=gmul(randomr(prec), Pi2n(1, prec));/*Random angle*/
  GEN zbound=expIr(arg);/*The boundary point. Now we need to scale by a random hyperbolic distance in [0, R]*/
  /*a(r)=Area of hyperbolic disc radius r=4*Pi*sinh^2(r/2).*/
  GEN r=gmulsg(2, gasinh(gmul(gsinh(gdivgs(R, 2), prec), gsqrt(randomr(prec), prec)), prec));
  GEN e2r=gexp(r, prec);
  return gerepileupto(av, gmul(zbound, gdiv(gsubgs(e2r, 1), gaddgs(e2r, 1))));
}

/*Returns a random point z in the unit disc, uniform inside the ball of radius R, with argument uniform in [ang1, ang2]. See page 19 of Page (before section 2.5).*/
GEN
hdiscrandom_arc(GEN R, GEN ang1, GEN ang2, long prec)
{
  pari_sp av=avma;
  GEN arg=gadd(ang1, gmul(randomr(prec), gsub(ang2, ang1)));/*Random angle in [ang1, ang2]*/
  GEN zbound=expIr(arg);/*The boundary point. Now we need to scale by a random hyperbolic distance in [0, R]*/
  /*a(r)=Area of hyperbolic disc radius r=4*Pi*sinh^2(r/2).
  dist=gmul(gsqr(gsinh(gdivgs(R, 2), prec)), randomr(prec)); A random element in [0, a(R)/4Pi].
  r=gmulsg(2, gasinh(gsqrt(dist, prec), prec)); The radius
  */
  GEN r=gmulsg(2, gasinh(gmul(gsinh(gdivgs(R, 2), prec), gsqrt(randomr(prec), prec)), prec));
  GEN e2r=gexp(r, prec);
  return gerepileupto(av, gmul(zbound, gdiv(gsubgs(e2r, 1), gaddgs(e2r, 1))));
}

/*z1 and z2 are complex numbers, this computes the hyperbolic distance between them. If flag=0, assumes upper half plane model; if flag=1, assumes unit disc model.*/
GEN
hdist(GEN z1, GEN z2, long flag, long prec)
{
  pari_sp av=avma;
  if(flag) return hdist_ud(z1, z2, prec);
  GEN x1=real_i(z1), y1=imag_i(z1);
  GEN x2=real_i(z2), y2=imag_i(z2);
  GEN x=gaddsg(1, gdiv(gadd(gsqr(gsub(x2, x1)), gsqr(gsub(y2, y1))), gmul(gmulsg(2, y1), y2)));
  GEN expd=gadd(x, gsqrt(gsubgs(gsqr(x), 1), prec));
  return gerepileupto(av, glog(expd, prec));  
}

/*The hyperbolic distance between z1 and z2 in the unit disc model*/
static GEN
hdist_ud(GEN z1, GEN z2, long prec)
{
  pari_sp av=avma;
  GEN a=gabs(gsubsg(1, gmul(z1, conj_i(z2))), prec);/*|1-z1*conj(z2)|*/
  GEN b=gabs(gsub(z1, z2), prec);/*|z1-z2|*/
  GEN num=gadd(a, b);
  GEN denom=gsub(a, b);
  if(gequal0(denom))
  {
    pari_warn(warner, "You may not have enough precision to compute the hyperbolic distance");
    set_avma(av);
    return mkoo();
  }
  return gerepileupto(av, glog(gdiv(num, denom), prec));/*log((a+b)/(a-b))*/
}


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
  switch(typ(d))
  {
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
  switch(typ(x))
  {
    case t_FRAC:/*t_FRAC cannot be 0*/
	  return 0;
	case t_INT:/*Given exactly*/
	  return !signe(x);
	case t_REAL:
	  if(abscmprr(x, tol)<0) return 1;/*|x|<tol*/
	  return 0;
	case t_COMPLEX:;
	  long i;
	  for(i=1;i<=2;i++)
	  {
		switch(typ(gel(x, i)))
		{
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





