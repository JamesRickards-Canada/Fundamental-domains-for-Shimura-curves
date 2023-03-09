/*TO DO
1. Is it more efficient to convert to unit ball, act there, and convert back? OR ALSO SHOULD I STORE MUCH MORE DATA WITH IT?
2. Do I remove lists section? Should I use hash tables?
3. CHANGE THE LONG DECLARATIONS OUT OF FOR LOOPS
*/

/*
POSSIBLE FUTURE ADDITIONS:
1) Parallelization of element enumeration along with partial domain computations.
2) Methods of Imbert (see Voight's original paper) to find a minimal presentation in canonical form.
3) Computation of cohomology groups.
*/


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

/*SECTION 0: MISCELLANEOUS METHODS*/

/*SECTION 1: GEOMETRIC METHODS*/

/*1: INTERSECTION OF LINES/CIRCLES*/
static GEN line_int(GEN l1, GEN l2, GEN tol);
static GEN line_int11(GEN l1, GEN l2, GEN tol);
static GEN line_line_detNULL(GEN l1, GEN l2, GEN tol);
static int onseg(GEN l, GEN p, GEN tol);
static GEN seg_int(GEN l1, GEN l2, GEN tol);

/*1: MATRIX ACTION ON GEOMETRY*/
static GEN psl_to_klein(GEN M, GEN gdat);

/*1: DISTANCES/AREA*/
static GEN hdist_ud(GEN z1, GEN z2, long prec);

/*1: TOLERANCE*/
static int tolcmp(GEN x, GEN y, GEN tol);
static int tolcmp_sort(void *data, GEN x, GEN y);
static int toleq(GEN x, GEN y, GEN tol);
static int toleq0(GEN x, GEN tol);

/*SECTION 2: FUNDAMENTAL DOMAIN GEOMETRY*/

/*2: ISOMETRIC CIRCLES*/
static GEN icirc_psu(GEN M, GEN gdat);
static GEN icirc_elt(GEN X, GEN g, GEN (*Xtopsl)(GEN, GEN, long), GEN gdat);


/*1: SHORT VECTORS IN LATTICES*/
//static GEN quadraticintegernf(GEN nf, GEN A, GEN B, GEN C, long prec);
//static GEN smallvectors_cholesky(GEN Q, GEN C, long maxelts, GEN condition, long prec);
//static GEN smallvectors_nfcondition(GEN A, GEN C, long maxelts, GEN condition, long prec);

//static GEN hpolygon_area(GEN circles, GEN vertices, GEN tol, long prec);

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



/*SECTION 0: MISCELLANEOUS METHODS*/

/*1: INFINITY */

/*1: LISTS*/

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



/*0: SHORT VECTORS IN LATTICES*/




/*SECTION 1: GEOMETRIC METHODS*/


/*LINES, SEGMENTS, TOLERANCE
Line   ->	[a, b, c]			Representing ax+by=c. We will normalize so that c=1 or 0. It is assumed that at least one of a, b is non-zero.
Segment->	[a, b, c, x0, x1]	[a, b, c] gives the line, which has start point x0 and ends at x1, which are complex. We do not allow segments going
								through oo.
GEN tol->	The tolerance, which MUST be of type t_REAL. The default choice is tol=deftol(prec). Two points are declared as equal if they are equal up
			to tolerance.
*/

/* GEOMETRIC DATA
We will need to deal with passing from the upper half plane model -> unit ball model -> Klein model, as well as doing computations with tolerance. For this, we will fix a "geometric data" argument:
gdat = 	[prec, tol, p, pc, pscale]
prec  ->	The precision we are working with everywhere, stored as a GEN. We normally want it to be a long, which can be done with prec[2] (as we 
			assume prec>0 and it does not overflow.
tol   ->	The tolerance, which should be initially set with tol=deftol(prec).
p     ->	The point in the upper half plane that is mapped to 0 in the unit ball model. This should be chosen to have trivial stabilizer in Gamma, 
			otherwise issues may arise. We convert it to type t_REAL.
pc    ->	The conjugate of p.
pscale->	1/(p-pc), which is required when we convert from the upper half plane to the Klein model.
*/

/*1: INTERSECTION OF LINES/CIRCLES*/


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

/*Line intersection where c1=c2=1, the typical case. l1 and l2 can have length 2 in this case.*/
static GEN
line_int11(GEN l1, GEN l2, GEN tol)
{
  pari_sp av=avma;
  GEN det=line_line_detNULL(l1, l2, tol);/*ad-bc*/
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



/*1: MATRIX ACTION ON GEOMETRY*/


/*EQUATIONS FOR ACTION
UPPER HALF PLANE->	M=[a, b;c, d] in (P)GL(2, R)^+ acts on z via (az+b)/(cz+d). In this program, we normally normalize so that det(M)=1.
UNIT DISC		->	M=[A, B;conj(B), conj(A)] in (P)SU(1, 1) acts on z in the same way as PGL. Note that |A|^2-|B|^2=1 as det(M)=1
KLEIN			->	M=[A, B] corresponding to the same (A, B) as for the unit disc action. The corresponding equation on a point in the Klein model is
					via Mz=(A^2z+B^2conj(z)+2AB)/(|A|^2+|B|^2+A*z*conj(B)+conj(A)*conj(z)*B.
					=(A(Az+B)+B(B*conj(z)+A))/(conj(B)(Az+B)+conj(A)(B*conj(z)+A)).
*/

/*This gives the action in the Klein model, as described above.*/
GEN
klein_act(GEN M, GEN z)
{
  pari_sp av=avma;
  GEN A=gel(M, 1), B=gel(M, 2);
  GEN AzpB=gadd(gmul(A, z), B);/*Az+B*/
  GEN BzcpA=gadd(gmul(B, conj_i(z)), A);/*B*conj(z)+A*/
  GEN num=gadd(gmul(A, AzpB), gmul(B, BzcpA));/*A(Az+B)+B(B*conj(z)+A)*/
  GEN denom=gadd(gmul(conj_i(B), AzpB), gmul(conj_i(A), BzcpA));/*conj(B)(Az+B)+conj(A)(B*conj(z)+A)*/
  return gerepileupto(av, gdiv(num, denom));
}

/*Gives the action of a matrix in the upper half plane/unit disc model. We assume that the input/output are not infinity, which could happen with the upper half plane model.*/
GEN
pgl_act(GEN M, GEN z){
  pari_sp av=avma;
  GEN numer=gadd(gmul(gcoeff(M, 1, 1), z), gcoeff(M, 1, 2));
  GEN denom=gadd(gmul(gcoeff(M, 2, 1), z), gcoeff(M, 2, 2));
  return gerepileupto(av, gdiv(numer, denom));
}

/*Coverts M in PSL(2, R) to [A, B] which acts on the Klein model. If M1=1/(p-conj(p))[1,-p;1,-conj(p)] and M2=[conj(p), -p;1, -1], then this is via M1*M*M2=[A, B;conj(B), conj(A)]. If M=[a, b;c, d], the explicit formula are A=(a*conj(p)-|p|^2*c+b-pd)/(p-conj(p)), and B=(-ap+p^2c-b+pd)/(p-conj(p)).*/
static GEN
psl_to_klein(GEN M, GEN gdat)
{
  pari_sp av=avma;
  GEN p=gdat_get_p(gdat), pc=gdat_get_pc(gdat), pscale=gdat_get_pscale(gdat);
  GEN ampc=gsub(gcoeff(M, 1, 1), gmul(p, gcoeff(M, 2, 1)));/*a-pc*/
  GEN bmpd=gsub(gcoeff(M, 1, 2), gmul(p, gcoeff(M, 2, 2)));/*b-pd*/
  GEN ampc_pconj=gmul(ampc, pc);/*(a-pc)*conj(p)*/
  GEN Apre=gadd(ampc_pconj, bmpd);/*(a-pc)*conj(p)+b-pd*/
  GEN ampc_p=gmul(ampc, p);/*(a-pc)*p*/
  GEN Bpre=gneg(gadd(ampc_p, bmpd));/*(-a+pc)p-b+pd*/
  GEN AB=cgetg(3, t_VEC);
  gel(AB, 1)=gmul(Apre, pscale);
  gel(AB, 2)=gmul(Bpre, pscale);
  return gerepileupto(av, AB);
}



/*1: DISTANCES/AREAS*/


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



/*1: TOLERANCE*/


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
	case t_REAL:
	  if(abscmprr(d, tol)<0) return 0;/*|d|<tol*/
	  return gc_int(av, signe(d));
	case t_INT:/*Given exactly*/
	  return gc_int(av, signe(d));
	case t_FRAC:/*t_FRAC cannot be 0*/
	  return gc_int(av, signe(gel(d, 1)));
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
	case t_INT:/*Given exactly*/
	  return !signe(x);
	case t_FRAC:/*t_FRAC cannot be 0*/
	  return 0;
  }
  pari_err_TYPE("Tolerance equality only valid for type t_INT, t_FRAC, t_REAL, t_COMPLEX", x);
  return 0;/*So that there is no warning*/
}



/*SECTION 2: FUNDAMENTAL DOMAIN GEOMETRY*/



/* DEALING WITH GENERAL STRUCTURES
Let X be a structure, from which Gamma, a discrete subgroup of PSL(2, R), can be extracted. Given a vector of elements, we want to be able to compute the normalized boundary, and the normalized basis. In order to achieve this, we need to pass in various methods to deal with operations in Gamma, namely:
	i) Multiply elements of Gamma: format as GEN Xmul(GEN X, GEN g1, GEN g2).
	ii) Invert elements of Gamma: format as GEN Xinv(GEN X, GEN g).
	iii) Embed elements of Gamma in PSL(2, R): format as GEN Xtopsl(GEN X, GEN g, long prec). We REQUIRE t_REAL entries in various places, so you should ensure that the coefficients of the output are converted to t_REAL.
	iv) Identify if an element is trivial in Gamma: format as int Xistriv(GEN X, GEN g). Since we are working in PSL, we need to be careful that -g==g, since most representations of elements in X would be for SL.
	v) Pass in the identity element of Gamma and find the area of the fundamental domain. These methods are not passed in; just the values.
We do all of our computations in the Klein model.
*/

/*2: ISOMETRIC CIRCLES*/


/* ISOMETRIC CIRCLE FORMATTING
icirc ->	[[a, b], r, ang]
[a, b]->	ax+by=1 is the Klein model equation of the boundary
			(x-a)^2+(y-b)^2=r^2 is the unit disc model equation of the boundary
r     ->	r^2+1=a^2+b^2 as it is orthogonal to the unit disc
ang   ->	icirc intersects the unit disc in P1 and P2, and assume that the angle from P1 to P2 is <pi (which uniquely determines them). Then ang is
			the argument of P1, in the range [0, 2*pi).
*/


/*Given an element M=[a, b;c, d] of PSU(1, 1), this returns the isometric circle associated to it. This has centre -d/c, radius 1/|c|. In the Klein model, the centre being u+iv -> xu+yv=1, and this intersects the unit disc at (a+/-br/(a^2+b^2), (b+/-ar)/(a^2+b^2)). If we pass in an element giving everything (i.e. c=0), we return 0.*/
static GEN
icirc_psu(GEN M, GEN gdat)
{
  GEN tol=gdat_get_tol(gdat);
  long prec=gdat_get_prec(gdat);
  if(toleq0(gcoeff(M, 2, 1), tol)) return gen_0;/*Isometric circle is everything, don't want to call it here.*/
  pari_sp av=avma;
  GEN centre=gneg(gdiv(gcoeff(M, 2, 2), gcoeff(M, 2, 1)));/*-d/c, the centre. If this is u+iv, then ux+vy=1 is the Klein version.*/
  GEN r=invr(gabs(gcoeff(M, 2, 1), prec));/*1/|c| is the centre.*/
  GEN a=real_i(centre), b=imag_i(centre);/*The coords of the centre.*/
  GEN ar=mulrr(a, r), br=mulrr(b, r);
  GEN apbr=addrr(a, br), ambr=subrr(a, br);/*a+/-br*/
  GEN bpar=addrr(b, ar), bmar=subrr(b, ar);/*b+/-ar*/
  GEN pi=mppi(prec);
  GEN theta1, theta2;/*The intersection point angles are tan^-1((b+ar)/(a-br)) and tan^-1((b-ar)/(a+br)). We normalize so they are between 0 and 2Pi*/
  if(toleq0(apbr, tol))
  {
	theta1=Pi2n(-1, prec);/*Pi/2*/
	if(signe(bmar)==-1) theta1=addrr(theta1, pi);/*3Pi/2*/
  }
  else
  {
	theta1=gatan(divrr(bmar, apbr), prec);
	if(signe(apbr)==-1) theta1=addrr(theta1, pi);/*atan lies in (-Pi/2, Pi/2), so must add by Pi to get it correct, as x=(a+br)/(a^2+b^2).*/
	else if(signe(bmar)==-1) theta1=addrr(theta1, Pi2n(1, prec));/*Add 2Pi to get in the right interval*/
  }
  if(toleq0(ambr, tol))
  {
	theta2=Pi2n(-1, prec);/*Pi/2*/
	if(signe(bpar)==-1) theta2=addrr(theta2, pi);/*3Pi/2*/
  }
  else
  {
	theta2=gatan(divrr(bpar, ambr), prec);
	if(signe(ambr)==-1) theta2=addrr(theta2, pi);/*atan lies in (-Pi/2, Pi/2), so must add by Pi to get it correct, as x=(a-br)/(a^2+b^2).*/
	else if(signe(bpar)==-1) theta2=addrr(theta2, Pi2n(1, prec));/*Add 2Pi to get in the right interval*/
  }
  GEN thetadiff=subrr(theta2, theta1);
  if(signe(thetadiff)==1)
  {
	if(cmprr(thetadiff, pi)>0) theta1=theta2;/*theta2 is actually the "first" angle, so we swap it in.*/
  }
  else
  {
	if(cmprr(thetadiff, negr(pi))>0) theta1=theta2;/*theta2 is actually the "first" angle, so we swap it in.*/  
  }
  return gerepilecopy(av, mkvec3(mkvec2(a, b), r, theta1));/*Return the output!*/
}

/*Computes the isometric circle for g in Gamma, returning [g, M, icirc], where M gives the action of g on the Klein model.*/
static GEN
icirc_elt(GEN X, GEN g, GEN (*Xtopsl)(GEN, GEN, long), GEN gdat)
{
  pari_sp av=avma;
  GEN psl=Xtopsl(X, g, gdat_get_prec(gdat));
  GEN ret=cgetg(4, t_VEC);
  gel(ret, 1)=gcopy(g);
  gel(ret, 2)=psl_to_klein(psl, gdat);
  gel(ret, 3)=icirc_psu(gel(ret, 2), gdat);
  return gerepileupto(av, ret);
}





/*FUNDAMENTAL DOMAIN OTHER COMPUTATIONS*/

/*GEOMETRIC HELPER METHODS*/







