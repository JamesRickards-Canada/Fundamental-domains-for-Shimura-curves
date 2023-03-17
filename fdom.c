/*TO DO
2. Do I remove lists section? Should I use hash tables?
3. CHANGE THE LONG DECLARATIONS OUT OF FOR LOOPS
4. Make the insertion methods in normbound inline? Not sure if this does anything or not.
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
static int angle_onarc(GEN a1, GEN a2, GEN a, GEN tol);

/*1: MATRIX ACTION ON GEOMETRY*/
static GEN gdat_initialize(GEN p, long prec);

/*1: TRANSFER BETWEEN MODELS*/
static GEN disc_to_klein(GEN z);
static GEN disc_to_plane(GEN z, GEN p);
static GEN klein_to_disc(GEN z, GEN tol, long prec);
static GEN klein_to_plane(GEN z, GEN p, GEN tol, long prec);
static GEN plane_to_disc(GEN z, GEN p);
static GEN plane_to_klein(GEN z, GEN p);
static GEN psl_to_klein(GEN M, GEN gdat);

/*1: DISTANCES/AREA*/
static GEN hdist_ud(GEN z1, GEN z2, long prec);

/*1: TOLERANCE*/
static int tolcmp(GEN x, GEN y, GEN tol);
static int tolcmp_sort(void *data, GEN x, GEN y);
static int toleq(GEN x, GEN y, GEN tol);
static int toleq0(GEN x, GEN tol);
static int tolsigne(GEN x, GEN tol);

/*SECTION 2: FUNDAMENTAL DOMAIN GEOMETRY*/

/*2: ISOMETRIC CIRCLES*/
static GEN icirc_angle(GEN c1, GEN c2, long prec);
static GEN icirc_klein(GEN M, GEN gdat);
static GEN icirc_elt(GEN X, GEN g, GEN (*Xtopsl)(GEN, GEN, long), GEN gdat);
static GEN argmod(GEN x, GEN y, GEN tol, long prec);

/*2: NORMALIZED BOUNDARY*/
static GEN normbound_icircs(GEN C, GEN gdat);
static int cmp_icircangle(void *nul, GEN c1, GEN c2);
static long normbound_icircs_bigr(GEN C, GEN order);
static void normbound_icircs_insinfinite(GEN elts, GEN vcors, GEN vargs, GEN curcirc, long *found);
static void normbound_icircs_insclean(GEN elts, GEN vcors, GEN vargs, GEN curcirc, long toins, long *found);
static void normbound_icircs_phase2(GEN elts, GEN vcors, GEN vargs, GEN curcirc, GEN firstab, GEN tol, long prec, long toins, long *found);
static GEN normbound_area(GEN C, long prec);

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
  if (*vind == *vlen) {/*Need to lengthen!*/
    *vlen = 2**vlen;
    v = vec_lengthen(v, *vlen);
  }
  *vind = *vind+1;
  gel(v, *vind) = x;
  return v;
}

/*Appends x to v, returning v, and updating vind to vind++. If vind++>vlen, then we double the length of v as well. Don't forget to call vec_shorten at the end, since some positions are uninitialized.*/
GEN
vecsmalllist_append(GEN v, long *vind, long *vlen, long x)
{
  if (*vind == *vlen) {/*Need to lengthen!*/
    *vlen = 2**vlen;
    v = vecsmall_lengthen(v, *vlen);
  }
  *vind = *vind+1;
  v[*vind] = x;
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
  pari_sp av = avma;
  GEN c1 = seg_get_c(l1), c2 = seg_get_c(l2);/*Get c values*/
  GEN det = line_line_detNULL(l1, l2, tol);/*ad-bc*/
  if (!det) return gc_NULL(av);/*Lines are parallel*/
  if (!signe(c1)) {/*ax+by = 0*/
    if (gequal0(c2)) return gc_const(av, gen_0);/*Both must pass through 0*/
	return gerepilecopy(av, mkcomplex(gdiv(gneg(seg_get_b(l1)), det), gdiv(seg_get_a(l1), det)));/*-b/det, a/det*/
  }
  /*Next, cx+dy = 0*/
  if (!signe(c2)) return gerepilecopy(av, mkcomplex(gdiv(seg_get_b(l2), det), gdiv(gneg(seg_get_a(l2)), det)));/*d/det, -c/det*/
  /*Now ax+by = cx+dy = 1*/
  GEN x = gdiv(gsub(seg_get_b(l2), seg_get_b(l1)), det);/*(d-b)/det*/
  GEN y = gdiv(gsub(seg_get_a(l1), seg_get_a(l2)), det);/*(a-c)/det*/
  return gerepilecopy(av, mkcomplex(x, y));
}

/*Line intersection where c1=c2=1, the typical case. l1 and l2 can have length 2 in this case.*/
static GEN
line_int11(GEN l1, GEN l2, GEN tol)
{
  pari_sp av = avma;
  GEN det = line_line_detNULL(l1, l2, tol);/*ad-bc*/
  if(!det) return gc_NULL(av);/*Lines are concurrent*/
  GEN x = gdiv(gsub(seg_get_b(l2), seg_get_b(l1)), det);/*(d-b)/det*/
  GEN y = gdiv(gsub(seg_get_a(l1), seg_get_a(l2)), det);/*(a-c)/det*/
  return gerepilecopy(av, mkcomplex(x, y));
}

/*Given two lines given by a1x+b1y=c1 and a2x+b2y=c2, returns a1b2-a2b1, unless this is within tolerance of 0, when we return NULL. Not stack clean.*/
static GEN
line_line_detNULL(GEN l1, GEN l2, GEN tol){
  pari_sp av = avma;
  GEN d = gsub(gmul(seg_get_a(l1), seg_get_b(l2)), gmul(seg_get_b(l1), seg_get_a(l2)));
  if (toleq0(d, tol)) return gc_NULL(av);/*d = 0 up to tolerance*/
  return d;
}

/*p is assumed to be on the line defined by l; this checks if it is actually on the segment l. Returns 0 if not, 1 if p is in the interior, 2 if p is the start point, and 3 if p is the end point (all up to tolerance tol).*/
static int
onseg(GEN l, GEN p, GEN tol)
{
  if (lg(l) == LINE_LG) return 1;/*Lines contain all points*/
  pari_sp av = avma;
  GEN pstart = seg_get_start(l), px = real_i(p);
  int xcmp1 = tolcmp(real_i(pstart), px, tol);/*Compare x coeff of starting point and x*/
  if (!xcmp1) {/*Same x-value*/
	GEN py = imag_i(p);
	int ycmp1 = tolcmp(imag_i(pstart), py, tol);
	if (!ycmp1) return gc_int(av, 2);/*We must be the starting point.*/
	GEN pend = seg_get_end(l);/*At this point, the slope must be oo*/
	int ycmp2 = tolcmp(py, imag_i(pend), tol);
	if (!ycmp2) return gc_int(av, 3);/*Since we assume that p is on the line, we must be at endpoint 2 now.*/
	if (ycmp1 == ycmp2) return gc_int(av, 1);/*Must be between*/
	return gc_int(av, 0);/*Above or below*/
  }/*Now we have distinct x-values*/
  int xcmp2 = tolcmp(px, real_i(seg_get_end(l)), tol);
  if (xcmp2 == 0) return gc_int(av, 3);/*End point, as the slope cannot be oo as pstartx! = pendx*/
  if (xcmp1 == xcmp2) return gc_int(av, 1);
  return gc_int(av, 0);
}

/*Returns the intersection of two segments. If parallel/concurrent/do not intersect, returns NULL*/
static GEN 
seg_int(GEN l1, GEN l2, GEN tol)
{
  pari_sp av = avma;
  GEN ipt = line_int(l1, l2, tol);/*Intersect the lines*/
  if (!onseg(l1, ipt, tol)) return gc_NULL(av);
  if (!onseg(l2, ipt, tol)) return gc_NULL(av);
  return ipt;
}

/*A is the arc on the unit disc from angle a1 to a2, where a1 and a2 are in [0, 2Pi). For another ange a in [0, 2Pi), this returns 1 if a==a1, 2 if a==a2, 3 if a in in the interior of the counterclockwise arc from a1 to a2, and 0 if it is outside of this arc. All angles must be t_REAL, even if they are 0.*/
static int
angle_onarc(GEN a1, GEN a2, GEN a, GEN tol)
{
  int a1cmp = tolcmp(a1, a, tol);
  if (!a1cmp) return 1;/*Equal to a1*/
  int cross = cmprr(a1, a2);/*No need for tolerance here*/
  if (cross > 0){/*The arc crosses the x-axis.*/
	if (a1cmp < 0) return 3;/*Between a1 and the x-axis*/
	int a2cmp = tolcmp(a, a2, tol);
	if (!a2cmp) return 2;/*Equal to a2*/
	if (a2cmp < 0) return 3;/*Inside*/
	return 0;/*Outside*/
  }
  if (a1cmp > 0) return 0;/*Bewteen x-axis and a1*/
  int a2cmp = tolcmp(a, a2, tol);
  if (!a2cmp) return 2;/*Equal to a2*/
  if (a2cmp > 0) return 0;/*Outside*/
  return 3;/*Inside*/
}



/*1: MATRIX ACTION ON GEOMETRY*/


/*EQUATIONS FOR ACTION
UPPER HALF PLANE->	M=[a, b;c, d] in (P)GL(2, R)^+ acts on z via (az+b)/(cz+d). In this program, we normally normalize so that det(M)=1.
UNIT DISC		->	M=[A, B;conj(B), conj(A)] in (P)SU(1, 1) acts on z in the same way as PGL. Note that |A|^2-|B|^2=1 as det(M)=1
KLEIN			->	M=[A, B] corresponding to the same (A, B) as for the unit disc action. The corresponding equation on a point in the Klein model is
					via Mz=(A^2z+B^2conj(z)+2AB)/(|A|^2+|B|^2+A*z*conj(B)+conj(A)*conj(z)*B.
					=(A(Az+B)+B(B*conj(z)+A))/(conj(B)(Az+B)+conj(A)(B*conj(z)+A)).
*/

/*Initializes gdat for a given p and precision.*/
static GEN
gdat_initialize(GEN p, long prec)
{
  pari_sp av = avma;
  GEN tol = deftol(prec);
  GEN pc = conj_i(p);
  GEN pmpc = gsub(p, pc);
  GEN pscale = ginv(pmpc);
  return gerepilecopy(av, mkvec5(stoi(prec), tol, p, pc, pscale));
}

/*This gives the action in the Klein model, as described above.*/
GEN
klein_act(GEN M, GEN z)
{
  pari_sp av = avma;
  GEN A = gel(M, 1), B = gel(M, 2);
  GEN AzpB = gadd(gmul(A, z), B);/*Az+B*/
  GEN BzcpA = gadd(gmul(B, conj_i(z)), A);/*B*conj(z)+A*/
  GEN num = gadd(gmul(A, AzpB), gmul(B, BzcpA));/*A(Az+B)+B(B*conj(z)+A)*/
  GEN denom = gadd(gmul(conj_i(B), AzpB), gmul(conj_i(A), BzcpA));/*conj(B)(Az+B)+conj(A)(B*conj(z)+A)*/
  return gerepileupto(av, gdiv(num, denom));
}

/*Gives the action of a matrix in the upper half plane/unit disc model. We assume that the input/output are not infinity, which could happen with the upper half plane model.*/
GEN
pgl_act(GEN M, GEN z){
  pari_sp av = avma;
  GEN numer = gadd(gmul(gcoeff(M, 1, 1), z), gcoeff(M, 1, 2));
  GEN denom = gadd(gmul(gcoeff(M, 2, 1), z), gcoeff(M, 2, 2));
  return gerepileupto(av, gdiv(numer, denom));
}



/*1: TRANSFER BETWEEN MODELS*/


/*Given a point z in the unit disc model, this transfers it to the Klein model.*/
static GEN
disc_to_klein(GEN z)
{
  pari_sp av = avma;
  GEN znorm = gnorm(z);
  GEN scale = gdivsg(2, gaddsg(1, znorm));//2/(1+|z|^2)
  return gerepileupto(av, gmul(scale, z));/*2z/(1+|z|^2)*/
}

/*Given a point z in the unit disc model, this transfers it to the upper half plane model.*/
static GEN
disc_to_plane(GEN z, GEN p)
{
  pari_sp av = avma;
  GEN num = gsub(gmul(conj_i(p), z), p);/*conj(p)*z-p*/
  GEN denom = gsubgs(z, 1);/*z-1*/
  return gerepileupto(av, gdiv(num, denom));
}

/*Given a point z in the Klein model, this transfers it to the unit disc model.*/
static GEN
klein_to_disc(GEN z, GEN tol, long prec)
{
  pari_sp av = avma;
  if (!tol) tol = deftol(prec);
  GEN znm1 = gsubsg(1, gnorm(z)), rt;/*1-|z|^2*/
  if (toleq0(znm1, tol)) rt = gen_0;/*sqrt(0) can cause great precision loss, so we check for equality with 0 before square rooting.*/
  else rt = gsqrt(znm1, prec);/*sqrt(1-|z|^2)*/
  GEN scale = invr(gaddsg(1, rt));//1/(1+sqrt(1-|z|^2))
  return gerepileupto(av, gmul(scale, z));/*z/(1+sqrt(1-|z|^2))*/
}

/*Given a point z in the Klein model, this transfers it to the upper half plane model.*/
static GEN
klein_to_plane(GEN z, GEN p, GEN tol, long prec)
{
  pari_sp av = avma;
  GEN zdisc = klein_to_disc(z, tol, prec);
  return gerepileupto(av, disc_to_plane(zdisc, p));//Klein -> disc -> plane
}

/*Given a point z in the upper half plane model, this transfers it to the unit disc model.*/
static GEN
plane_to_disc(GEN z, GEN p)
{
  pari_sp av = avma;
  GEN num = gsub(z, p);/*z-p*/
  GEN denom = gsub(z, conj_i(p));/*z-conj(p)*/
  return gerepileupto(av, gdiv(num, denom));
}

/*Given a point z in the upper half plane model, this transfers it to the Klein model.*/
static GEN
plane_to_klein(GEN z, GEN p)
{
  pari_sp av = avma;
  GEN zdisc = plane_to_disc(z, p);
  return gerepileupto(av, disc_to_klein(zdisc));/*Plane -> disc -> Klein*/
}

/*Coverts M in PSL(2, R) to [A, B] which acts on the Klein model. If M1=1/(p-conj(p))[1,-p;1,-conj(p)] and M2=[conj(p), -p;1, -1], then this is via M1*M*M2=[A, B;conj(B), conj(A)]. If M=[a, b;c, d], the explicit formula are A=(a*conj(p)-|p|^2*c+b-pd)/(p-conj(p)), and B=(-ap+p^2c-b+pd)/(p-conj(p)).*/
static GEN
psl_to_klein(GEN M, GEN gdat)
{
  pari_sp av = avma;
  GEN p = gdat_get_p(gdat), pc = gdat_get_pc(gdat), pscale = gdat_get_pscale(gdat);
  GEN ampc = gsub(gcoeff(M, 1, 1), gmul(p, gcoeff(M, 2, 1)));/*a-pc*/
  GEN bmpd = gsub(gcoeff(M, 1, 2), gmul(p, gcoeff(M, 2, 2)));/*b-pd*/
  GEN ampc_pconj = gmul(ampc, pc);/*(a-pc)*conj(p)*/
  GEN Apre = gadd(ampc_pconj, bmpd);/*(a-pc)*conj(p)+b-pd*/
  GEN ampc_p = gmul(ampc, p);/*(a-pc)*p*/
  GEN Bpre = gneg(gadd(ampc_p, bmpd));/*(-a+pc)p-b+pd*/
  GEN AB = cgetg(3, t_VEC);
  gel(AB, 1) = gmul(Apre, pscale);
  gel(AB, 2) = gmul(Bpre, pscale);
  return gerepileupto(av, AB);
}


/*1: DISTANCES/AREAS*/


/*Given the area of a hyperbolic disc, this returns the radius. The formula is area=4*Pi*sinh(R/2)^2, or R=2arcsinh(sqrt(area/4Pi))*/
GEN
hdiscradius(GEN area, long prec)
{
  pari_sp av = avma;
  return gerepileupto(av, gtofp(gmulgs(gasinh(gsqrt(gdiv(area, Pi2n(2, prec)), prec), prec), 2), prec));
}

/*Returns a random point z in the unit disc, uniform inside the ball of radius R. See page 19 of Page (before section 2.5).*/
GEN
hdiscrandom(GEN R, long prec)
{
  pari_sp av = avma;
  GEN arg = gmul(randomr(prec), Pi2n(1, prec));/*Random angle*/
  GEN zbound = expIr(arg);/*The boundary point. Now we need to scale by a random hyperbolic distance in [0, R]*/
  /*a(r) = Area of hyperbolic disc radius r = 4*Pi*sinh^2(r/2).*/
  GEN r = gmulsg(2, gasinh(gmul(gsinh(gdivgs(R, 2), prec), gsqrt(randomr(prec), prec)), prec));
  GEN e2r = gexp(r, prec);
  return gerepileupto(av, gmul(zbound, gdiv(gsubgs(e2r, 1), gaddgs(e2r, 1))));
}

/*Returns a random point z in the unit disc, uniform inside the ball of radius R, with argument uniform in [ang1, ang2]. See page 19 of Page (before section 2.5).*/
GEN
hdiscrandom_arc(GEN R, GEN ang1, GEN ang2, long prec)
{
  pari_sp av = avma;
  GEN arg = gadd(ang1, gmul(randomr(prec), gsub(ang2, ang1)));/*Random angle in [ang1, ang2]*/
  GEN zbound = expIr(arg);/*The boundary point. Now we need to scale by a random hyperbolic distance in [0, R]*/
  /*a(r) = Area of hyperbolic disc radius r = 4*Pi*sinh^2(r/2).
  dist = gmul(gsqr(gsinh(gdivgs(R, 2), prec)), randomr(prec)); A random element in [0, a(R)/4Pi].
  r = gmulsg(2, gasinh(gsqrt(dist, prec), prec)); The radius
  */
  GEN r = gmulsg(2, gasinh(gmul(gsinh(gdivgs(R, 2), prec), gsqrt(randomr(prec), prec)), prec));
  GEN e2r = gexp(r, prec);
  return gerepileupto(av, gmul(zbound, gdiv(gsubgs(e2r, 1), gaddgs(e2r, 1))));
}

/*z1 and z2 are complex numbers, this computes the hyperbolic distance between them. If flag=0, assumes upper half plane model; if flag=1, assumes unit disc model.*/
GEN
hdist(GEN z1, GEN z2, long flag, long prec)
{
  pari_sp av = avma;
  if (flag) return hdist_ud(z1, z2, prec);
  GEN x1 = real_i(z1), y1 = imag_i(z1);
  GEN x2 = real_i(z2), y2 = imag_i(z2);
  GEN x = gaddsg(1, gdiv(gadd(gsqr(gsub(x2, x1)), gsqr(gsub(y2, y1))), gmul(gmulsg(2, y1), y2)));
  GEN expd = gadd(x, gsqrt(gsubgs(gsqr(x), 1), prec));
  return gerepileupto(av, glog(expd, prec));  
}

/*The hyperbolic distance between z1 and z2 in the unit disc model*/
static GEN
hdist_ud(GEN z1, GEN z2, long prec)
{
  pari_sp av = avma;
  GEN a = gabs(gsubsg(1, gmul(z1, conj_i(z2))), prec);/*|1-z1*conj(z2)|*/
  GEN b = gabs(gsub(z1, z2), prec);/*|z1-z2|*/
  GEN num = gadd(a, b);
  GEN denom = gsub(a, b);
  if (gequal0(denom)) {
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
  pari_sp av = avma;
  return gerepileupto(av, shiftr(gtofp(gen_1, prec), BITS_IN_LONG/2*(2-prec)));
}

/*Returns -1 if x<y, 0 if x==y, 1 if x>y (x, y are t_REAL/t_INT/t_FRAC). Accounts for the tolerance, so will deem x==y if they are equal up to tol AND at least one is inexact*/
static int
tolcmp(GEN x, GEN y, GEN tol)
{
  pari_sp av = avma;
  GEN d = gsub(x, y);
  return gc_int(av, tolsigne(d, tol));/*Return sign(x-y)*/
}

/*Data points to tol. Used to sort/search a list with tolerance.*/
static int
tolcmp_sort(void *data, GEN x, GEN y){return tolcmp(x, y, *(GEN*)data);}

/*Returns 1 if x==y up to tolerance tol. If x and y are both t_INT/t_FRAC, will only return 1 if they are exactly equal.*/
static int
toleq(GEN x, GEN y, GEN tol)
{
  pari_sp av = avma;
  GEN d = gsub(x, y);
  return gc_int(av, toleq0(d, tol));/*Just compare d with 0.*/
}

/*If x is of type t_INT or t_FRAC, returns 1 iff x==0. Otherwise, x must be of type t_REAL or t_COMPLEX, and returns 1 iff x=0 up to tolerance tol.*/
static int
toleq0(GEN x, GEN tol)
{
  switch (typ(x)) {
	case t_REAL:
	  if (abscmprr(x, tol) < 0) return 1;/*|x|<tol*/
	  return 0;
	case t_COMPLEX:;
	  long i;
	  for (i = 1; i <= 2; i++) {
		switch (typ(gel(x, i))) {
		  case t_FRAC:/*Fraction component, cannot be 0*/
		    return 0;
		  case t_INT:
		    if (signe(gel(x, i))) return 0;
			break;
		  case t_REAL:
		    if (abscmprr(gel(x, i), tol) >= 0) return 0;/*Too large*/
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

/*Returns the sign of x, where it is 0 if is equal to 0 up to tolerance.*/
static int
tolsigne(GEN x, GEN tol)
{
  switch (typ(x)) {
	case t_REAL:
	  if (abscmprr(x, tol) < 0) return 0;/*|x|<tol*/
	case t_INT:/*Given exactly or real and not 0 up to tolerance.*/
	  return signe(x);
	case t_FRAC:/*t_FRAC cannot be 0*/
	  return signe(gel(x, 1));
  }
  pari_err_TYPE("Tolerance comparison only valid for type t_INT, t_FRAC, t_REAL", x);
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
icirc  ->	[[a, b], r, p1, p2, ang1, ang2]
[a, b] ->	ax+by=1 is the Klein model equation of the boundary
			(x-a)^2+(y-b)^2=r^2 is the unit disc model equation of the boundary
r      ->	r^2+1=a^2+b^2 as it is orthogonal to the unit disc
p1/p2  ->	Intersection points with the unit disc; see below for the ordering.  
ang1   ->	Assume that the angle from p1 to p2 is <pi (which uniquely determines them). Then ang1 is the argument of p1, in the range [0, 2*pi).
ang2   ->	The argument of p2, in the range [0, 2*pi).
a, b, r, ang1, ang2 must be t_REAL.
p1 and p2 must be t_COMPLEX with t_REAL components.
*/


/*Given two isometric circles c1 and c2, assumed to intersect, this computes the angle they form (counterclockwise from the tangent to c1 to the tangent to c2 at the intersection point).
If c1 is given by (x-a)^2+(y-b)^2=r^2 and c2 by (x-c)^2+(y-d)^2=s^2 (in the unit ball model) and theta is the angled formed, then a long computation shows that cos(theta)=(1-ac-bd)/rs.
This method is not stack clean.
*/
static GEN
icirc_angle(GEN c1, GEN c2, long prec)
{
  GEN ac = mulrr(gmael(c1, 1, 1), gmael(c2, 1, 1));
  GEN bd = mulrr(gmael(c1, 1, 2), gmael(c2, 1, 2));
  GEN omacmbd = subsr(1, addrr(ac, bd));/*1-ac-bd*/
  GEN cost = divrr(omacmbd, mulrr(gel(c1, 2), gel(c2, 2)));/*cos(theta)*/
  return gacos(cost, prec);/*acos is in the interval [0, Pi], which is what we want.*/
}

/*Given M=[A, B] acting on the Klein model, this returns the isometric circle associated to it. This has centre -conj(A/B), radius 1/|B|. In the Klein model, the centre being a+bi -> xa+yb=1, and this intersects the unit disc at (a+/-br/(a^2+b^2), (b-/+ar)/(a^2+b^2)). Assumptions:
	If we pass in an element giving everything (i.e. B=0), we return 0.
	A and B are t_REAL OR t_COMPLEX with t_REAL components
*/
static GEN
icirc_klein(GEN M, GEN gdat)
{
  GEN tol = gdat_get_tol(gdat);
  long prec = gdat_get_prec(gdat);
  if (toleq0(gel(M, 2), tol)) return gen_0;/*Isometric circle is everything, don't want to call it here.*/
  pari_sp av = avma;
  GEN centre_pre = gdiv(gel(M, 1), gel(M, 2));/*The centre is -conj(centre_pre)*/
  GEN r = invr(gabs(gel(M, 2), prec));/*1/|B| is the radius.*/
  GEN a = negr(real_i(centre_pre)), b = imag_i(centre_pre);/*The coords of the centre.*/
  GEN ar = mulrr(a, r), br = mulrr(b, r);
  GEN apbr = addrr(a, br), ambr = subrr(a, br);/*a+/-br*/
  GEN bpar = addrr(b, ar), bmar = subrr(b, ar);/*b+/-ar*/
  GEN theta1 = argmod(apbr, bmar, tol, prec);/*First point angle in [0, 2pi)*/
  GEN theta2 = argmod(ambr, bpar, tol, prec);/*Second point angle in [0, 2pi)*/
  GEN asqbsq = addrs(sqrr(r), 1);/*a^2+b^2 = r^2+1*/
  GEN p1 = mkcomplex(divrr(apbr, asqbsq), divrr(bmar, asqbsq));/*First intersection point.*/
  GEN p2 = mkcomplex(divrr(ambr, asqbsq), divrr(bpar, asqbsq));/*Second intersection point.*/
  GEN thetadiff = subrr(theta2, theta1);
  if (signe(thetadiff) == 1) {/*If the next if block executes, theta2 is actually the "first" angle, so we swap them.*/
	if (cmprr(thetadiff, mppi(prec)) > 0) gerepilecopy(av, mkvecn(6, mkvec2(a, b), r, p2, p1, theta2, theta1));
  }
  else {/*Same as above*/
	if (cmprr(thetadiff, negr(mppi(prec))) > 0) gerepilecopy(av, mkvecn(6, mkvec2(a, b), r, p2, p1, theta2, theta1));
  }
  return gerepilecopy(av, mkvecn(6, mkvec2(a, b), r, p1, p2, theta1, theta2));/*Return the output!*/
}

/*Computes the isometric circle for g in Gamma, returning [g, M, icirc], where M gives the action of g on the Klein model.*/
static GEN
icirc_elt(GEN X, GEN g, GEN (*Xtopsl)(GEN, GEN, long), GEN gdat)
{
  pari_sp av = avma;
  GEN psl = Xtopsl(X, g, gdat_get_prec(gdat));
  GEN ret = cgetg(4, t_VEC);
  gel(ret, 1) = gcopy(g);
  gel(ret, 2) = psl_to_klein(psl, gdat);
  gel(ret, 3) = icirc_klein(gel(ret, 2), gdat);
  return gerepileupto(av, ret);
}

/*Returns the argument of x+iy in the range [0, 2*pi). Assumes x and y are not both 0. Not stack clean.*/
static GEN
argmod(GEN x, GEN y, GEN tol, long prec)
{
  int xsign = tolsigne(x, tol);/*Sign of x up to tolerance.*/
  if (xsign == 0) {/*Fixing theta when x == 0, so the line is vertical.*/
	GEN theta = Pi2n(-1, prec);/*Pi/2*/
	if (signe(y) == -1) theta = addrr(theta, mppi(prec));/*3Pi/2*/
	return theta;
  }
  GEN theta = gatan(gdiv(y, x), prec);
  if (xsign == -1) return addrr(theta, mppi(prec));/*atan lies in (-Pi/2, Pi/2), so must add by Pi to get it correct.*/
  int ysign = tolsigne(y, tol);
  if (ysign == 0) return gtofp(gen_0, prec);/*Convert to real number.*/
  if (ysign == -1) theta = addrr(theta, Pi2n(1, prec));/*Add 2Pi to get in the right interval*/
  return theta;
}


/*2: NORMALIZED BOUNDARY*/


/*NORMALIZED BOUNDARY FORMATTING
A normalized boundary is represented by U, where
U=[elts, sides, vcors, vargs, kact, area, spair, gdat]
elts     ->	Elements of Gamma whose isometric circles give the sides of the normalied boundary. An infinite side corresponds to the element 0.
sides    ->	The ith entry is the isometric circle corresponding to elts[i], stored as [[a, b], r] representing ax+by=1 and a^2+b^2=r^2+1. Infinite
			side -> 0
vcors    ->	Vertex coordinates, stored as a t_COMPLEX. The side sides[i] has vertices vcor[i-1] and vcor[i], appearing in this order going 
			counterclockwise about the origin.
vargs    ->	The argument of the corresponding vertex, normalized to lie in [0, 2*Pi). We will also shift the elements so that this set is sorted, i.e. 
			vangs[1] is the vertex with minimal argument.
kact     ->	The action of the corresponding element on the Klein model, for use in klein_act. Infinite side -> 0
area     ->	The hyperbolic area of U, which will be oo unless we are a finite index subgroup of Gamma.
spair    ->	Stores the side pairing of U, if it exists/has been computed. When computing the normalized boundary, this will be stored as 0.
gdat	 ->	Stores the geometric data associated to the computations.
*/


/*Given C, the output of icirc_elt for each of our elements, this computes and returns the normalized boundary. Assumptions:
	C[i]=[g, M, icirc], where g=elt, M=Kleinian action, and icirc=[[a, b], r, p1, p2, ang1, ang2].
	None of the entries should be infinite sides, sort this out in the method that calls this.
	C has at least one element, so the boundary is non-trivial.
	Not stack clean/gerepileupto safe.
*/
static GEN
normbound_icircs(GEN C, GEN gdat)
{
  long prec = gdat_get_prec(gdat);
  GEN tol = gdat_get_tol(gdat);
  GEN order = gen_indexsort(C, NULL, &cmp_icircangle);/*Order the sides by initial angle.*/
  long lc = lg(C), maxsides = 2*lc, istart = normbound_icircs_bigr(C, order);/*Largest r value index (wrt order).*/
  GEN elts = cgetg(maxsides, t_VECSMALL);/*Stores indices in C of the sides. 0 represents an infinite side.*/
  GEN vcors = cgetg(maxsides, t_VEC);/*Stores coordinates of the vertices.*/
  GEN vargs = cgetg(maxsides, t_VEC);/*Stores arguments of the vertices.*/
  elts[1] = order[istart];
  GEN firstcirc = gmael(C, elts[1], 3);/*The equation for the fist line, used in phase 2 and for detecting phase 2 starting.*/
  gel(vcors, 1) = gel(firstcirc, 4);/*The terminal vertex of the side is the first vertex.*/
  gel(vargs, 1) = gel(firstcirc, 6);
  long found = 1, lenc = lc-1, absind;/*found=how many sides we have found up to now. This can increase and decrease.*/
  int phase2 = 0, infinitesides = 0;/*Which phase we are in, and if there are infinite sides or not.*/
  /*PHASE 1: inserting sides, where the end vertex does not intersect the first side. All we need to update / keep track of are:
	elts, vcors, vargs, found, absind
	PHASE 2: The same, but we have looped back around, and need to insert the last edge into the first (which will never disappear).
  */
  for (absind = 1; absind < lenc; absind++) {/*We are trying to insert the next side.*/
	long toins = order[1+(istart+absind-1)%lenc];/*The index of the side to insert.*/
	GEN curcirc = gmael(C, toins, 3);/*The isometric circle we are trying to insert*/
	GEN lastcirc = gmael(C, elts[found], 3);/*The isometric circle of the last side inserted.*/
	int termloc = angle_onarc(gel(lastcirc, 5), gel(lastcirc, 6), gel(curcirc, 6), tol);/*Whether the terminal angle lies inside the last arc.*/
	switch (termloc) {
	  case 3:
	    if (angle_onarc(gel(lastcirc, 5), gel(lastcirc, 6), gel(curcirc, 5), tol)) continue;/*The initial angle also lies here, so we are totally enveloped, and we move on.*/
	  case 1:
	    phase2 = 1;/*We have looped back around and are intersecting from the right.*/
		infinitesides = 1;
		normbound_icircs_insinfinite(elts, vcors, vargs, curcirc, &found);/*Insert oo side*/
		normbound_icircs_phase2(elts, vcors, vargs, curcirc, gel(firstcirc, 1), tol, prec, toins, &found);/*Phase 2 insertion.*/
		continue;
	  case 2:
	    continue;/*The terminal angle is the same as the previous, so we are enveloped. It is not possible for our new side to envelop the old side.*/
	}
	/*If we make it here, then the terminal point lies outside of the last range.*/
	int initloc = angle_onarc(gel(lastcirc, 5), gel(lastcirc, 6), gel(curcirc, 5), tol);/*Whether the initial angle lies inside the last arc.*/
	switch (initloc) {
	  case 3:/*We have an intersection to consider, will do so outside the switch.*/
	    break;
	  case 1:
	    absind--;/*We completely envelop the previous side. We need to delete it and redo this case.*/
		found--;
		continue;
	  case 0:
	    infinitesides = 1;
	    normbound_icircs_insinfinite(elts, vcors, vargs, curcirc, &found);/*Insert oo side!*/
	  case 2:
	    if (phase2 || angle_onarc(gel(lastcirc, 5), gel(lastcirc, 6), gel(firstcirc, 5), tol)) {/*Phase2 has started*/
		  phase2 = 1;
		  normbound_icircs_phase2(elts, vcors, vargs, curcirc, gel(firstcirc, 1), tol, prec, toins, &found);/*Phase 2 insertion.*/
		  continue;
		}
	    normbound_icircs_insclean(elts, vcors, vargs, curcirc, toins, &found);/*New(initial)=Old(terminal), so there is no infinite side coming first*/
		continue;
	}
	/*If we make it here, our current circle intersects the last one, so we need to see if it is "better" than the previous intersection.*/
	GEN ipt = line_int11(gel(curcirc, 1), gel(lastcirc, 1), tol);/*Find the intersection point.*/
	GEN iptarg = argmod(real_i(ipt), imag_i(ipt), tol, prec);/*Argument*/
	if (found == 1) {/*Straight up insert it*/
	   gel(vcors, found) = ipt;
	   gel(vargs, found) = iptarg;/*Fix the last vertex*/
	   normbound_icircs_insclean(elts, vcors, vargs, curcirc, toins, &found);/*Clean insert*/
	   continue;
	}
	int newptdat = angle_onarc(gel(vargs, found-1), gel(vargs, found), iptarg, tol);/*The new point wrt the previous side.*/
	switch (newptdat) {/*newptdat=2 is impossible.*/
	  case 3:/*Insert it!*/
	    gel(vcors, found) = ipt;
	    gel(vargs, found) = iptarg;/*Fix the last vertex*/
		if (phase2 || angle_onarc(gel(lastcirc, 5), gel(lastcirc, 6), gel(firstcirc, 5), tol)) {/*Phase2 has started*/
		  phase2 = 1;
		  normbound_icircs_phase2(elts, vcors, vargs, curcirc, gel(firstcirc, 1), tol, prec, toins, &found);/*Phase 2 insertion.*/
		  continue;
		}
	    normbound_icircs_insclean(elts, vcors, vargs, curcirc, toins, &found);/*Clean insert*/
	    continue;
	  case 0:
	    if (!phase2) {/*We supercede the last side, so delete and try again.*/
	      absind--;
	      found--;
	      continue;
		}
		int supercede = angle_onarc(gel(curcirc, 5), gel(curcirc, 6), gel(lastcirc, 5), tol);/*We are in phase 2, and we either miss the side to the left (so ignore and move on), or to the right (supercede previous edge).*/
		if (supercede) {
		  absind--;
		  found--;
		}
		continue;
	  case 1:/*We just replace the last side, we have the same vertex and vertex angle.*/
	    elts[found]=toins;
		if(phase2 || angle_onarc(gel(lastcirc, 5), gel(lastcirc, 6), gel(firstcirc, 5), tol)) {/*Phase2 has started*/
		  gel(vcors, found) = line_int11(gel(curcirc, 1), gel(firstcirc, 1), tol);/*Intersect with initial side*/
          gel(vargs, found) = argmod(real_i(gel(vcors, found)), imag_i(gel(vcors, found)), tol, prec);/*Argument*/
		}
	}
  }
  if (!phase2) {/*We never hit phase 2, so there is one final infinite edge to worry about.*/
    infinitesides = 1;
	normbound_icircs_insinfinite(elts, vcors, vargs, firstcirc, &found);
  }
  /*Now we can compile everything into the return vector.*/
  long i;
  GEN rv = cgetg(9, t_VEC);
  GEN rv_elts = cgetg(found+1, t_VEC);
  GEN rv_sides = cgetg(found+1, t_VEC);
  GEN rv_kact = cgetg(found+1, t_VEC);
  for (i=1; i<=found; i++) {/*Make sure we treat infinite sides correctly!*/
	gel(rv_elts, i) = elts[i] ? gmael(C, elts[i], 1) : gen_0;
	gel(rv_sides, i) = elts[i] ? vec_shorten(gmael(C, elts[i], 3), 2) : gen_0;
	gel(rv_kact, i) = elts[i] ? gmael(C, elts[i], 2) : gen_0;
  }
  gel(rv, 1) = rv_elts;/*The elements*/
  gel(rv, 2) = rv_sides;/*The sides*/
  gel(rv, 3) = vec_shorten(vcors, found);/*Vertex coords*/
  gel(rv, 4) = vec_shorten(vargs, found);/*Vertex arguments*/
  gel(rv, 5) = rv_kact;/*Kleinian action*/
  if(infinitesides) gel(rv, 6)=mkoo();/*Infinite side means infinite area.*/
  else gel(rv, 6) = normbound_area(rv_sides, prec);
  gel(rv, 7) = gen_0;/*Side pairing*/
  gel(rv, 8) = gdat;/*Geometric data*/
  return rv;
}


/*Used for sorting C by initial angles.*/
static int
cmp_icircangle(void *nul, GEN c1, GEN c2){return cmprr(gmael(c1, 3, 5), gmael(c2, 3, 5));}

/*Finds the index i such that C[order[i]] has the largest r value, which is guaranteed to be on the normalized boundary.*/
static long
normbound_icircs_bigr(GEN C, GEN order)
{
  pari_sp av = avma;
  long i, best = 1;
  GEN maxr = gmael3(C, order[1], 3, 2);
  for (i = 2; i < lg(C); i++) {
	if (cmprr(gmael3(C, order[i], 3, 2), maxr) <= 0) continue;/*Not as big. No need for tolerance here, it won't affect anything.*/
	best = i;
	maxr = gmael3(C, order[i], 3, 2);
  }
  return gc_long(av, best);
}

/*curcirc gives a new side that does not loop back around, but there is an infinite side first. This inserts the infinite side.*/
static void
normbound_icircs_insinfinite(GEN elts, GEN vcors, GEN vargs, GEN curcirc, long *found)
{
  (*found)++;
  elts[*found] = 0;/*Infinite side*/
  gel(vcors, *found) = gel(curcirc, 3);/*Initial vertex of curcirc*/
  gel(vargs, *found) = gel(curcirc, 5);
}

/*We are inserting a new side that does not intersect with a previous one, except possibly at the terminal vertex of the last side.*/
static void
normbound_icircs_insclean(GEN elts, GEN vcors, GEN vargs, GEN curcirc, long toins, long *found)
{
  (*found)++;
  elts[*found] = toins;/*Non-infinite side we are inserting.*/
  gel(vcors, *found) = gel(curcirc, 4);/*Terminal vertex of curcirc*/
  gel(vargs, *found) = gel(curcirc, 6);
}

/*We are performing an insertion in phase 2, i.e. we are intersecting back with the initial side.*/
static void
normbound_icircs_phase2(GEN elts, GEN vcors, GEN vargs, GEN curcirc, GEN firstab, GEN tol, long prec, long toins, long *found)
{
  (*found)++;
  elts[*found] = toins;
  gel(vcors, *found) = line_int11(gel(curcirc, 1), firstab, tol);/*Intersect with initial side*/
  gel(vargs, *found) = argmod(real_i(gel(vcors, *found)), imag_i(gel(vcors, *found)), tol, prec);/*Argument*/
}

/*Returns the hyperbolic area of the normalized boundary, which is assumed to not have any infinite sides (we keep track if they exist, and do not call this method if they do). C should be the list of [[a, b], r] in order. The area is (n-2)*Pi-sum(angles), where there are n sides.*/
static GEN
normbound_area(GEN C, long prec){
  pari_sp av = avma;
  long n = lg(C)-1, i;
  GEN area = mulsr(n-2, mppi(prec));/*(n-2)*Pi*/
  for (i=1; i<n; i++) area = subrr(area, icirc_angle(gel(C, i), gel(C, i+1), prec));
  area = subrr(area, icirc_angle(gel(C, n), gel(C, 1), prec));
  return gerepileupto(av, area);
}


/*Used to suppress warnings as build the package.*/
void warningholder()
{
  klein_to_plane(gen_0, gen_0, gen_0, 3);
  plane_to_klein(gen_0, gen_0);
}

