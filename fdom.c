/*TO DO
3. CHANGE THE LONG DECLARATIONS OUT OF FOR LOOPS
5. Check distances/area section: I don't really use this yet. It's just kind of there.
6. How to input groups between O^1 and the full positive normalizer group?
8. forqfvec with flag=2.
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

/*STATIC DECLARATIONS*/

/*SECTION 1: GEOMETRIC METHODS*/

/*1: LINES AND ARCS*/
static int angle_onarc(GEN a1, GEN a2, GEN a, GEN tol);
static GEN line_int(GEN l1, GEN l2, GEN tol);
static GEN line_int11(GEN l1, GEN l2, GEN tol);
static GEN line_line_detNULL(GEN l1, GEN l2, GEN tol);
static int onseg(GEN l, GEN p, GEN tol);
static GEN seg_int(GEN l1, GEN l2, GEN tol);

/*1: MATRIX ACTION ON GEOMETRY*/
static GEN defp(long prec);
static GEN gdat_initialize(GEN p, long prec);
static GEN klein_safe(GEN M, long prec);
static GEN uhp_safe(GEN p, long prec);

/*1: TRANSFER BETWEEN MODELS*/
static GEN disc_to_klein(GEN z);
static GEN disc_to_plane(GEN z, GEN p);
static GEN klein_to_disc(GEN z, GEN tol);
static GEN klein_to_plane(GEN z, GEN p, GEN tol);
static GEN plane_to_disc(GEN z, GEN p);
static GEN plane_to_klein(GEN z, GEN p);
static GEN m2r_to_klein(GEN M, GEN p);

/*1: DISTANCES/AREA*/
static GEN hdist_ud(GEN z1, GEN z2, long prec);

/*1: OPERATIONS ON COMPLEX REALS*/
static GEN gtocr(GEN z, long prec);
static GEN abscr(GEN z);
static GEN addcrcr(GEN z1, GEN z2);
static GEN divcrcr(GEN z1, GEN z2);
static GEN mulcrcr(GEN z1, GEN z2);
static GEN mulcrcr_conj(GEN z1, GEN z2);
static GEN mulcri(GEN z, GEN n);
static GEN mulcrr(GEN z, GEN r);
static GEN negc(GEN z);
static GEN normcr(GEN z);
static GEN subcrcr(GEN z1, GEN z2);
static GEN subrcr(GEN r, GEN z);

/*1: TOLERANCE*/
static int tolcmp(GEN x, GEN y, GEN tol);
static int toleq(GEN x, GEN y, GEN tol);
static int toleq0(GEN x, GEN tol);
static int tolsigne(GEN x, GEN tol);

/*SECTION 2: FUNDAMENTAL DOMAIN GEOMETRY*/

/*2: ISOMETRIC CIRCLES*/
static GEN icirc_angle(GEN c1, GEN c2, long prec);
static GEN icirc_klein(GEN M, GEN tol);
static GEN icirc_elt(GEN X, GEN g, GEN (*Xtoklein)(GEN, GEN), GEN gdat);
static GEN argmod(GEN x, GEN y, GEN tol);
static GEN argmod_complex(GEN c, GEN tol);

/*2: NORMALIZED BOUNDARY*/
static GEN normbound(GEN X, GEN G, GEN (*Xtoklein)(GEN, GEN), GEN gdat);
static GEN normbound_icircs(GEN C, GEN indtransfer, GEN gdat);
static long normbound_icircs_bigr(GEN C, GEN order);
static void normbound_icircs_insinfinite(GEN elts, GEN vcors, GEN vargs, GEN infinite, GEN curcirc, long *found);
static void normbound_icircs_insclean(GEN elts, GEN vcors, GEN vargs, GEN curcirc, long toins, long *found);
static void normbound_icircs_phase2(GEN elts, GEN vcors, GEN vargs, GEN curcirc, GEN firstcirc, GEN tol, long toins, long *found);
static GEN normbound_area(GEN C, long prec);
static long normbound_outside(GEN U, GEN z, GEN tol);
static int normbound_whichside(GEN side, GEN z, GEN tol);

/*2: NORMALIZED BOUNDARY APPENDING*/
static GEN normbound_append(GEN X, GEN U, GEN G, GEN (*Xtoklein)(GEN, GEN), GEN gdat);
static GEN normbound_append_icircs(GEN Uvcors, GEN Uvargs, GEN C, GEN Ctype, long rbigind, GEN gdat);

/*2: NORMALIZED BOUNDARY ANGLES*/
static int cmp_icircangle(void *nul, GEN c1, GEN c2);
static long args_find_cross(GEN args);
static long args_search(GEN args, long ind, GEN arg, GEN tol);

/*2: NORMALIZED BASIS*/
static GEN edgepairing(GEN U, GEN tol);
static GEN normbasis(GEN X, GEN U, GEN G, GEN (*Xtoklein)(GEN, GEN), GEN (*Xmul)(GEN, GEN, GEN), GEN (*Xinv)(GEN, GEN), int (*Xistriv)(GEN, GEN), GEN gdat);

/*2: NORMALIZED BOUNDARY REDUCTION*/
static GEN red_elt_decomp(GEN X, GEN U, GEN g, GEN z, GEN (*Xtoklein)(GEN, GEN), GEN (*Xmul)(GEN, GEN, GEN), GEN gdat);
static GEN red_elt(GEN X, GEN U, GEN g, GEN z, GEN (*Xtoklein)(GEN, GEN), GEN (*Xmul)(GEN, GEN, GEN), int flag, GEN gdat);

/*SECTION 3: QUATERNION ALGEBRA METHODS*/

/*3: INITIALIZE SYMMETRIC SPACE*/
static GEN afuchinit_i(GEN A, GEN O, GEN type, GEN p, long prec);
static GEN afuch_make_kleinmats(GEN A, GEN O, GEN p, long prec);
static GEN afuch_make_m2rmats(GEN A, GEN O, long prec);
static GEN Omultable(GEN A, GEN O);

/*3: ALGEBRA FUNDAMENTAL DOMAIN METHODS*/

/*3: ALGEBRA BASIC AUXILLARY METHODS*/
static GEN afuchid(GEN X);
static GEN afuchinv(GEN X, GEN g);
static int afuchistriv(GEN X, GEN g);
static GEN afuchmul(GEN X, GEN g1, GEN g2);
static GEN afuchtoklein(GEN X, GEN g);

/*3: ALGEBRA HELPER METHODS*/
static GEN algconj(GEN A, GEN x);
static GEN alggeta(GEN A);
static GEN voidalgmul(void *A, GEN x, GEN y);
static long algsplitoo(GEN A);



/*SECTION 1: GEOMETRIC METHODS*/


/*LINES, SEGMENTS, TOLERANCE, POINTS
Line	[a, b, c]
	Representing ax+by=c. We will normalize so that c=gen_1 or gen_0. It is assumed that at least one of a, b is non-zero. We also assume that a and b are of type t_REAL.
Point
	Stored as t_COMPLEX with t_REAL entries.
Segment	[a, b, c, x0, x1]
	[a, b, c] gives the line, which has start point x0 and ends at x1, which are complex. We do not allow segments going through oo. We also require x0 and x1 to have type t_COMPLEX with t_REAL components.
tol
	The tolerance, which MUST be of type t_REAL. The default choice is tol=deftol(prec). Two points are declared as equal if they are equal up to tolerance.
*/

/* GEOMETRIC DATA
We will need to deal with passing from the upper half plane model -> unit ball model -> Klein model, as well as doing computations with tolerance. For this, we will fix a "geometric data" argument:
gdat = [tol, p]
tol
	The tolerance, which should be initially set with tol=deftol(prec).
p
	The point in the upper half plane that is mapped to 0 in the unit ball model. This should be chosen to have trivial stabilizer in Gamma, otherwise issues may arise. We convert it to have components that are t_REAL.
precision
	Can be retrieved by lg(tol), so we don't store it.
*/


/*1: LINES AND ARCS*/

/*A is the arc on the unit disc from angle a1 to a2, where a1 and a2 are in [0, 2Pi). For another angle a in [0, 2Pi), this returns 1 if a==a1, 2 if a==a2, 3 if a in in the interior of the counterclockwise arc from a1 to a2, and 0 if it is outside of this arc. All angles must be t_REAL, even if they are 0.*/
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
  if (a1cmp > 0) return 0;/*Between x-axis and a1*/
  int a2cmp = tolcmp(a, a2, tol);
  if (!a2cmp) return 2;/*Equal to a2*/
  if (a2cmp > 0) return 0;/*Outside*/
  return 3;/*Inside*/
}

/*Returns the intersection of two lines. If parallel/concurrent, returns NULL. We will ensure that all output numbers are COMPLEX with REAL components. This means that 0 is stored as a real number with precision.*/
static GEN
line_int(GEN l1, GEN l2, GEN tol)
{
  pari_sp av = avma;
  GEN c1 = gel(l1, 3), c2 = gel(l2, 3);/*Get c values*/
  GEN det = line_line_detNULL(l1, l2, tol);/*ad-bc*/
  if (!det) return gc_NULL(av);/*Lines are parallel*/
  if (!signe(c1)) {/*ax+by = 0*/
    if (!signe(c2)) {/*Both must pass through 0*/
	  GEN zeror = real_0(lg(tol));/*prec=lg(tol)*/
	  return gerepilecopy(av, mkcomplex(zeror, zeror));
	}
	return gerepilecopy(av, mkcomplex(divrr(negr(gel(l1, 2)), det), divrr(gel(l1, 1), det)));/*-b/det, a/det*/
  }
  /*Next, cx+dy = 0*/
  if (!signe(c2)) return gerepilecopy(av, mkcomplex(divrr(gel(l2, 2), det), divrr(negr(gel(l2, 1)), det)));/*d/det, -c/det*/
  /*Now ax+by = cx+dy = 1*/
  GEN x = divrr(subrr(gel(l2, 2), gel(l1, 2)), det);/*(d-b)/det*/
  GEN y = divrr(subrr(gel(l1, 1), gel(l2, 1)), det);/*(a-c)/det*/
  return gerepilecopy(av, mkcomplex(x, y));
}

/*Line intersection where c1=c2=1, the typical case. l1 and l2 can have length 2 in this case.*/
static GEN
line_int11(GEN l1, GEN l2, GEN tol)
{
  pari_sp av = avma;
  GEN det = line_line_detNULL(l1, l2, tol);/*ad-bc*/
  if(!det) return gc_NULL(av);/*Lines are concurrent*/
  GEN x = divrr(subrr(gel(l2, 2), gel(l1, 2)), det);/*(d-b)/det*/
  GEN y = divrr(subrr(gel(l1, 1), gel(l2, 1)), det);/*(a-c)/det*/
  return gerepilecopy(av, mkcomplex(x, y));
}

/*Given two lines given by a1x+b1y=c1 and a2x+b2y=c2, returns a1b2-a2b1, unless this is within tolerance of 0, when we return NULL. Gerepileupto safe, leaves garbage if NULL.*/
static GEN
line_line_detNULL(GEN l1, GEN l2, GEN tol)
{
  GEN d = subrr(mulrr(gel(l1, 1), gel(l2, 2)), mulrr(gel(l1, 2), gel(l2, 1)));/*ad-bc*/
  if (toleq0(d, tol)) return NULL;/*d=0 up to tolerance*/
  return d;
}

/*p is assumed to be on the line segment defined by l; this checks if it is actually on the segment l. Returns 0 if not, 1 if p is in the interior, 2 if p is the start point, and 3 if p is the end point (all up to tolerance tol).*/
static int
onseg(GEN l, GEN p, GEN tol)
{
  pari_sp av = avma;
  GEN pstart = gel(l, 4), px = gel(p, 1);
  int xcmp1 = tolcmp(gel(pstart, 1), px, tol);/*Compare x coeff of starting point and x*/
  if (!xcmp1) {/*Same x-value*/
	GEN py = gel(p, 2);
	int ycmp1 = tolcmp(gel(pstart, 2), py, tol);
	if (!ycmp1) return gc_int(av, 2);/*We must be the starting point.*/
	GEN pend = gel(l, 5);/*At this point, the slope must be oo*/
	int ycmp2 = tolcmp(py, gel(pend, 2), tol);
	if (!ycmp2) return gc_int(av, 3);/*Since we assume that p is on the line, we must be at endpoint 2 now.*/
	if (ycmp1 == ycmp2) return gc_int(av, 1);/*Must be between*/
	return gc_int(av, 0);/*Above or below*/
  }/*Now we have distinct x-values*/
  int xcmp2 = tolcmp(px, gel(gel(l, 5), 1), tol);
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


/*1: MATRIX ACTION ON GEOMETRY*/

/*EQUATIONS FOR ACTION
UPPER HALF PLANE
	M=[a, b;c, d] in (P)GL(2, R)^+ acts on z via (az+b)/(cz+d). 
UNIT DISC
	M=[A, B;conj(B), conj(A)] in (P)SU(1, 1) acts on z in the same way as PGL. 
KLEIN
	M=[A, B] corresponding to the same (A, B) as for the unit disc action. The corresponding equation on a point in the Klein model is via 
	Mz=(A^2z+B^2conj(z)+2AB)/(|A|^2+|B|^2+A*z*conj(B)+conj(A)*conj(z)*B
	=(A(Az+B)+B(B*conj(z)+A))/(conj(B)(Az+B)+conj(A)(B*conj(z)+A)).
*/

/*Returns the default value of p, which is Pi/8+0.5*I*/
static GEN
defp(long prec)
{
  GEN p = cgetg(3, t_COMPLEX);
  gel(p, 1) = Pi2n(-3, prec);
  gel(p, 2) = real2n(-1, prec);
  return p;
}

/*Initializes gdat for a given p and precision. */
static GEN
gdat_initialize(GEN p, long prec)
{
  GEN ret = cgetg(3, t_VEC);
  gel(ret, 1) = deftol(prec);/*Default tolerance*/
  gel(ret, 2) = uhp_safe(p, prec);/*Convert p to have real components, and check that it was a valid input.*/
  return ret;
}

/*This gives the action in the Klein model, as described above.*/
GEN
klein_act_i(GEN M, GEN z)
{
  pari_sp av = avma;
  GEN A = gel(M, 1), B = gel(M, 2);
  GEN AzpB = addcrcr(mulcrcr(A, z), B);/*Az+B*/
  GEN BzcpA = addcrcr(mulcrcr_conj(B, z), A);/*B*conj(z)+A*/
  GEN num = addcrcr(mulcrcr(A, AzpB), mulcrcr(B, BzcpA));/*A(Az+B)+B(B*conj(z)+A)*/
  GEN denom = addcrcr(mulcrcr_conj(AzpB, B), mulcrcr_conj(BzcpA, A));/*(Az+B)conj(B)+(B*conj(z)+A)conj(A)*/
  return gerepilecopy(av, divcrcr(num, denom));
}

/*The safe version of klein_act_i, where we ensure z and M have the correct format.*/
GEN
klein_act(GEN M, GEN z, long prec)
{
  pari_sp av = avma;
  GEN Msafe = klein_safe(M, prec);
  GEN zsafe = gtocr(z, prec);
  return gerepileupto(av, klein_act_i(Msafe, zsafe));
}

/*Returns M=[A, B] converted to having complex entries with real components of precision prec.*/
static GEN
klein_safe(GEN M, long prec)
{
  if (typ(M) != t_VEC || lg(M) != 3) pari_err_TYPE("M must be a length 2 vector with complex entries.", M);
  GEN Msafe = cgetg(3, t_VEC);
  gel(Msafe, 1) = gtocr(gel(M, 1), prec);
  gel(Msafe, 2) = gtocr(gel(M, 2), prec);
  return Msafe;
}

/*Gives the action of a matrix in the upper half plane/unit disc model. We assume that the input/output are not infinity, which could happen with the upper half plane model.*/
GEN
pgl_act(GEN M, GEN z)
{
  pari_sp av = avma;
  GEN numer = gadd(gmul(gcoeff(M, 1, 1), z), gcoeff(M, 1, 2));
  GEN denom = gadd(gmul(gcoeff(M, 2, 1), z), gcoeff(M, 2, 2));
  return gerepileupto(av, gdiv(numer, denom));
}

/*p should be a point in the upper half plane, stored as a t_COMPLEX with t_REAL components. This converts p to this format, or raises an error if it is not convertible.*/
static GEN
uhp_safe(GEN p, long prec)
{
  if (typ(p) != t_COMPLEX) pari_err_TYPE("The point p must be a complex number with positive imaginary part.", p);
  GEN psafe = gtocr(p, prec);
  if (signe(gel(psafe, 2)) != 1) pari_err_TYPE("The point p must have positive imaginary part.", p);
  return psafe;
}


/*1: TRANSFER BETWEEN MODELS*/

/*Given a point z in the unit disc model, this transfers it to the Klein model.*/
static GEN
disc_to_klein(GEN z)
{
  pari_sp av = avma;
  GEN znorm = normcr(z);
  GEN scale = divsr(2, addsr(1, znorm));//2/(1+|z|^2)
  return gerepileupto(av, mulcrr(z, scale));/*2z/(1+|z|^2)*/
}

/*Given a point z in the unit disc model, this transfers it to the upper half plane model. The formula is conj(p)*z-p/(z-1)*/
static GEN
disc_to_plane(GEN z, GEN p)
{
  pari_sp av = avma;
  GEN num = subcrcr(mulcrcr_conj(z, p), p);/*z*conj(p)-p*/
  GEN denom = mkcomplex(subrs(gel(z, 1), 1), gel(z, 2));/*z-1*/
  return gerepileupto(av, divcrcr(num, denom));
}

/*Given a point z in the Klein model, this transfers it to the unit disc model. The formula is z/(1+sqrt(1-|z|^2)).*/
static GEN
klein_to_disc(GEN z, GEN tol)
{
  pari_sp av = avma;
  GEN znm1 = subsr(1, normcr(z));/*1-|z|^2*/
  if (toleq0(znm1, tol)) return gerepilecopy(av, z);/*z->z, so just copy and return it to avoid precision loss with sqrt(0).*/
  GEN rt = sqrtr(znm1);/*sqrt(1-|z|^2)*/
  GEN scale = invr(addsr(1, rt));//1/(1+sqrt(1-|z|^2))
  return gerepileupto(av, mulcrr(z, scale));/*z/(1+sqrt(1-|z|^2))*/
}

/*Given a point z in the Klein model, this transfers it to the upper half plane model.*/
static GEN
klein_to_plane(GEN z, GEN p, GEN tol)
{
  pari_sp av = avma;
  GEN zdisc = klein_to_disc(z, tol);
  return gerepileupto(av, disc_to_plane(zdisc, p));//Klein -> disc -> plane
}

/*Given a point z in the upper half plane model, this transfers it to the unit disc model. The formula is (z-p)/(z-conj(p))*/
static GEN
plane_to_disc(GEN z, GEN p)
{
  pari_sp av = avma;
  GEN num = subcrcr(z, p);/*z-p*/
  GEN denom = subcrcr(z, conj_i(p));/*z-conj(p)*/
  return gerepileupto(av, divcrcr(num, denom));
}

/*Given a point z in the upper half plane model, this transfers it to the Klein model.*/
static GEN
plane_to_klein(GEN z, GEN p)
{
  pari_sp av = avma;
  GEN zdisc = plane_to_disc(z, p);
  return gerepileupto(av, disc_to_klein(zdisc));/*Plane -> disc -> Klein*/
}

/*Coverts M in M(2, R) to [A, B] which corresponds to upper half plane model -action > unit disc/Klein model action when the matrices have positive determinat (det([A, B])=|A|^2-|B|^2). If M1=1/(p-conj(p))[1,-p;1,-conj(p)] and M2=[conj(p), -p;1, -1], then this is via M1*M*M2=[A, B;conj(B), conj(A)]. If M=[a, b;c, d], the explicit formula are A=(a*conj(p)-|p|^2*c+b-pd)/(p-conj(p)), and B=(-ap+p^2c-b+pd)/(p-conj(p)), with p=x+iy -> 1/(p-conj(p))=-i/2y. We ensure that A and B are of type t_COMPLEX with t_REAL coefficients.*/
static GEN
m2r_to_klein(GEN M, GEN p)
{
  pari_sp av = avma;
  GEN scale = invr(shiftr(gel(p, 2), 1));/*1/2y; -1 acts trivially*/
  GEN ampc = subrcr(gcoeff(M, 1, 1), mulcrr(p, gcoeff(M, 2, 1)));/*a-pc*/
  GEN bmpd = subrcr(gcoeff(M, 1, 2), mulcrr(p, gcoeff(M, 2, 2)));/*b-pd*/
  GEN ampc_pconj = mulcrcr_conj(ampc, p);/*(a-pc)*conj(p)*/
  GEN Apre = mulcrr(addcrcr(ampc_pconj, bmpd), scale);/*((a-pc)*conj(p)+b-pd)/2y*/
  GEN ampc_p = mulcrcr(ampc, p);/*(a-pc)*p*/
  GEN Bpre = mulcrr(negc(addcrcr(ampc_p, bmpd)), scale);/*((-a+pc)p-b+pd)/2y*/
  GEN AB = mkvec2(mkcomplex(negr(gel(Apre, 2)), gel(Apre, 1)), mkcomplex(negr(gel(Bpre, 2)), gel(Bpre, 1)));/*times i*/
  return gerepilecopy(av, AB);
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


/*1: OPERATIONS ON COMPLEX REALS*/

/*z should be a complex number with real components of type prec. This converts it to one. Clean method.*/
static GEN
gtocr(GEN z, long prec)
{
  GEN zsafe = cgetg(3, t_COMPLEX);
  if (typ(z) == t_COMPLEX) {
    gel(zsafe, 1) = gtofp(gel(z, 1), prec);
    gel(zsafe, 2) = gtofp(gel(z, 2), prec);
	return zsafe;
  }
  gel(zsafe, 1) = gtofp(z, prec);
  gel(zsafe, 2) = real_0(prec);
  return zsafe;
}

/*Returns abs(z) for a complex number z with real components. Gerepileupto safe, leaves garbage.*/
static GEN
abscr(GEN z)
{
  GEN nm=normcr(z);
  return sqrtr(nm);
}

/*Adds two complex numbers with real components, giving a complex output. Clean method.*/
static GEN
addcrcr(GEN z1, GEN z2)
{
  GEN z = cgetg(3, t_COMPLEX);
  gel(z, 1) = addrr(gel(z1, 1), gel(z2, 1));
  gel(z, 2) = addrr(gel(z1, 2), gel(z2, 2));
  return z;
}

/*Divides two complex numbers with real components, giving a complex output. Gerepileupto safe, leaves garbage.*/
static GEN
divcrcr(GEN z1, GEN z2)
{
  GEN num = mulcrcr_conj(z1, z2);/*z1/z2=(z1*conj(z2))/norm(z2)*/
  GEN den = normcr(z2);
  GEN z = cgetg(3, t_COMPLEX);
  gel(z, 1) = divrr(gel(num, 1), den);
  gel(z, 2) = divrr(gel(num, 2), den);
  return z;
}

/*Multiplies two complex numbers with real components, giving a complex output. Gerepileupto safe, leaves garbage.*/
static GEN
mulcrcr(GEN z1, GEN z2)
{
  GEN ac = mulrr(gel(z1, 1), gel(z2, 1));
  GEN ad = mulrr(gel(z1, 1), gel(z2, 2));
  GEN bc = mulrr(gel(z1, 2), gel(z2, 1));
  GEN bd = mulrr(gel(z1, 2), gel(z2, 2));
  GEN z = cgetg(3, t_COMPLEX);
  gel(z, 1) = subrr(ac, bd);
  gel(z, 2) = addrr(ad, bc);
  return z;
}

/*mulcrcr, except we do z1*conj(z2). Gerepileupto safe, leaves garbage.*/
static GEN
mulcrcr_conj(GEN z1, GEN z2)
{
  GEN ac = mulrr(gel(z1, 1), gel(z2, 1));
  GEN ad = mulrr(gel(z1, 1), gel(z2, 2));
  GEN bc = mulrr(gel(z1, 2), gel(z2, 1));
  GEN bd = mulrr(gel(z1, 2), gel(z2, 2));
  GEN z = cgetg(3, t_COMPLEX);
  gel(z, 1) = addrr(ac, bd);
  gel(z, 2) = subrr(bc, ad);
  return z;
}

/*Multiplies a complex number with real components by an integer. Clean method.*/
static GEN
mulcri(GEN z, GEN n)
{
  GEN zn = cgetg(3, t_COMPLEX);
  gel(zn, 1) = mulri(gel(z, 1), n);
  gel(zn, 2) = mulri(gel(z, 2), n);
  return zn;
}

/*Multiplies a complex number with real components by a real number. Clean method.*/
static GEN
mulcrr(GEN z, GEN r)
{
  GEN zp = cgetg(3, t_COMPLEX);
  gel(zp, 1) = mulrr(gel(z, 1), r);
  gel(zp, 2) = mulrr(gel(z, 2), r);
  return zp;
}

/*Returns -z, for a complex number with real components. Clean method.*/
static GEN
negc(GEN z)
{
  GEN zp = cgetg(3, t_COMPLEX);
  gel(zp, 1) = negr(gel(z, 1));
  gel(zp, 2) = negr(gel(z, 2));
  return zp;
}

/*Norm of a complex number with real components, giving a real output. Gerepileupto safe, leaves garbage.*/
static GEN
normcr(GEN z)
{
  GEN x = sqrr(gel(z, 1));
  GEN y = sqrr(gel(z, 2));
  return addrr(x, y);
}

/*Subtracts two complex numbers with real components, giving a complex output. Clean method.*/
static GEN
subcrcr(GEN z1, GEN z2)
{
  GEN z = cgetg(3, t_COMPLEX);
  gel(z, 1) = subrr(gel(z1, 1), gel(z2, 1));
  gel(z, 2) = subrr(gel(z1, 2), gel(z2, 2));
  return z;
}

/*Does r-z, where r is real and z is complex with real components. Clean method.*/
static GEN
subrcr(GEN r, GEN z)
{
  GEN zp = cgetg(3, t_COMPLEX);
  gel(zp, 1) = subrr(r, gel(z, 1));
  gel(zp, 2) = negr(gel(z, 2));
  return zp;
}


/*1: TOLERANCE*/

/*Returns the default tolerance given the precision, which is saying that x==y if they are equal up to half of the precision.*/
GEN
deftol(long prec)
{
  return real2n(BITS_IN_LONG/2*(2 - prec), prec);
}

/*Returns -1 if x<y, 0 if x==y, 1 if x>y (x, y are t_REAL/t_INT/t_FRAC). Accounts for the tolerance, so will deem x==y if they are equal up to tol AND at least one is inexact*/
static int
tolcmp(GEN x, GEN y, GEN tol)
{
  pari_sp av = avma;
  GEN d = gsub(x, y);
  return gc_int(av, tolsigne(d, tol));/*Return sign(x-y)*/
}

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
		  case t_REAL:
		    if (abscmprr(gel(x, i), tol) >= 0) return 0;/*Too large*/
			break;
		  case t_INT:
		    if (signe(gel(x, i))) return 0;
			break;
		  case t_FRAC:/*Fraction component, cannot be 0*/
		    return 0;
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
Let X be a structure, from which Gamma, a discrete subgroup of PSL(2, R), can be extracted. Given a vector of elements, we want to be able to compute the normalized boundary, and the normalized basis. In order to achieve this, we need to pass in various methods to deal with operations in Gamma (that should be memory clean), namely:
	i) Multiply elements of Gamma: format as GEN Xmul(GEN X, GEN g1, GEN g2).
	ii) Invert elements of Gamma: format as GEN Xinv(GEN X, GEN g). We don't actually need the inverse, as long as g*Xinv(X, g) embeds to a constant multiple of the identity we are fine. In particular, we can use conjugation in a quaternion algebra.
	iii) Embed elements of Gamma into PGU(1, 1)^+ that act on the Klein model: format as GEN Xtoklein(GEN X, GEN g). The output should be [A, B], where the corresponding element in PGU(1, 1)^+ is [A, B;conj(B), conj(A)], and |A|^2-|B|^2 is positive. If |A|^2-|B|^2!=1, this is OK, and do not rescale it, as this will be more inefficient and may cause precision loss. If you have a natural embedding into PGL(2, R)^+, then you can use m2r_to_klein, which will ensure that det(M in PGL)=|A|^2-|B|^2.
	iv) Identify if an element is trivial in Gamma: format as int Xistriv(GEN X, GEN g). Since we are working in PSL, we need to be careful that -g==g, since most representations of elements in X would be for SL. Furthermore, if we do not return the true inverse in Xinv, then we have to account for this.
	v) Pass in the identity element of Gamma and find the area of the fundamental domain. These methods are not passed in; just the values.
We do all of our computations in the Klein model.
*/


/*2: ISOMETRIC CIRCLES*/

/* ISOMETRIC CIRCLE FORMATTING
icirc
	[a, b, r, p1, p2, ang1, ang2]
a, b
	ax+by=1 is the Klein model equation of the boundary. a and b are t_REAL
	(x-a)^2+(y-b)^2=r^2 is the unit disc model equation of the boundary
r
	r^2+1=a^2+b^2 as it is orthogonal to the unit disc, and is of type t_REAL
p1/p2
	Intersection points with the unit disc; see below for the ordering. Of type t_COMPLEX with t_REAL components.
ang1
	Assume that the angle from p1 to p2 is <pi (which uniquely determines them). Then ang1 is the argument of p1, in the range [0, 2*pi). Type t_REAL
ang2
	The argument of p2, in the range [0, 2*pi). It may be less than ang1 if the arc crosses 1. Type t_REAL.
*/


/*Given two isometric circles c1 and c2, assumed to intersect, this computes the angle they form (counterclockwise from the tangent to c1 to the tangent to c2 at the intersection point).
If c1 is given by (x-a)^2+(y-b)^2=r^2 and c2 by (x-c)^2+(y-d)^2=s^2 (in the unit ball model) and theta is the angle formed, then a long computation shows that cos(theta)=(1-ac-bd)/rs.
This method is not stack clean.
*/
static GEN
icirc_angle(GEN c1, GEN c2, long prec)
{
  GEN ac = mulrr(gel(c1, 1), gel(c2, 1));
  GEN bd = mulrr(gel(c1, 2), gel(c2, 2));
  GEN omacmbd = subsr(1, addrr(ac, bd));/*1-ac-bd*/
  GEN cost = divrr(omacmbd, mulrr(gel(c1, 3), gel(c2, 3)));/*cos(theta)*/
  return gacos(cost, prec);/*acos is in the interval [0, Pi], which is what we want.*/
}

/*Given M=[A, B] with |A|^2-|B|^2=w acting on the Klein model, this returns the isometric circle associated to it. This has centre -conj(A/B), radius sqrt(w)/|B|. In the Klein model, the centre being a+bi -> xa+yb=1, and this intersects the unit disc at (a+/-br/(a^2+b^2), (b-/+ar)/(a^2+b^2)). Assumptions:
	If we pass in an element giving everything (i.e. B=0), we return 0.
	A and B are t_COMPLEX with t_REAL components
*/
static GEN
icirc_klein(GEN M, GEN tol)
{
  if (toleq0(gel(M, 2), tol)) return gen_0;/*Isometric circle is everything, don't want to call it here.*/
  pari_sp av = avma;
  long prec = lg(tol);
  GEN A = gel(M, 1), B = gel(M, 2);
  GEN centre_pre = divcrcr(A, B);/*The centre is -conj(centre_pre)*/
  GEN Bnorm = normcr(B);
  GEN w = subrr(normcr(A), Bnorm);
  GEN r = sqrtr(divrr(w, Bnorm));/*sqrt(w)/|B| is the radius*/
  GEN a = gel(centre_pre, 1), b = gel(centre_pre, 2);/*The coords of the centre.*/
  togglesign(a);/*centre_pre=x+iy, then the centre is -x+iy*/
  GEN ar = mulrr(a, r), br = mulrr(b, r);
  GEN apbr = addrr(a, br), ambr = subrr(a, br);/*a+/-br*/
  GEN bpar = addrr(b, ar), bmar = subrr(b, ar);/*b+/-ar*/
  GEN theta1 = argmod(apbr, bmar, tol);/*First point angle in [0, 2pi)*/
  GEN theta2 = argmod(ambr, bpar, tol);/*Second point angle in [0, 2pi)*/
  GEN asqbsq = addrs(sqrr(r), 1);/*a^2+b^2 = r^2+1*/
  GEN p1 = mkcomplex(divrr(apbr, asqbsq), divrr(bmar, asqbsq));/*First intersection point.*/
  GEN p2 = mkcomplex(divrr(ambr, asqbsq), divrr(bpar, asqbsq));/*Second intersection point.*/
  GEN thetadiff = subrr(theta2, theta1);
  if (signe(thetadiff) == 1) {/*If the next if block executes, theta2 is actually the "first" angle, so we swap them.*/
	if (cmprr(thetadiff, mppi(prec)) > 0) return gerepilecopy(av, mkvecn(7, a, b, r, p2, p1, theta2, theta1));
  }
  else {/*Same as above*/
	if (cmprr(thetadiff, negr(mppi(prec))) > 0) return gerepilecopy(av, mkvecn(7, a, b, r, p2, p1, theta2, theta1));
  }
  return gerepilecopy(av, mkvecn(7, a, b, r, p1, p2, theta1, theta2));/*Return the output!*/
}

/*Computes the isometric circle for g in Gamma, returning [g, M, icirc], where M gives the action of g on the Klein model.*/
static GEN
icirc_elt(GEN X, GEN g, GEN (*Xtoklein)(GEN, GEN), GEN gdat)
{
  pari_sp av = avma;
  GEN tol = gdat_get_tol(gdat);
  GEN ret = cgetg(4, t_VEC);
  gel(ret, 1) = gcopy(g);
  gel(ret, 2) = Xtoklein(X, g);
  gel(ret, 3) = icirc_klein(gel(ret, 2), tol);
  return gerepileupto(av, ret);
}

/*Returns the argument of x+iy in the range [0, 2*pi). Assumes x and y are not both 0 and are t_REAL. Gerepileupto safe, leaves garbage.*/
static GEN
argmod(GEN x, GEN y, GEN tol)
{
  long prec = lg(tol);
  int xsign = tolsigne(x, tol);/*Sign of x up to tolerance.*/
  if (xsign == 0) {/*Fixing theta when x == 0, so the line is vertical.*/
	GEN theta = Pi2n(-1, prec);/*Pi/2*/
	if (signe(y) == -1) theta = addrr(theta, mppi(prec));/*3Pi/2*/
	return theta;
  }
  GEN theta = gatan(gdiv(y, x), prec);
  if (xsign == -1) return addrr(theta, mppi(prec));/*atan lies in (-Pi/2, Pi/2), so must add by Pi to get it correct.*/
  int ysign = tolsigne(y, tol);
  if (ysign == 0) return real_0(prec);/*Convert to real number.*/
  if (ysign == -1) theta = addrr(theta, Pi2n(1, prec));/*Add 2Pi to get in the right interval*/
  return theta;
}

/*argmod, except we take c=x+iy. Gerepileupto safe, leaves garbage.*/
static GEN
argmod_complex(GEN c, GEN tol) {return argmod(gel(c, 1), gel(c, 2), tol);}


/*2: NORMALIZED BOUNDARY*/

/*NORMALIZED BOUNDARY FORMATTING
A normalized boundary is represented by U, where
U=[elts, sides, vcors, vargs, crossind, kact, area, spair, deleted]
elts
	Elements of Gamma whose isometric circles give the sides of the normalied boundary. An infinite side corresponds to the element 0. Elements cannot be stored as type t_INT.
sides
	The ith entry is the isometric circle corresponding to elts[i], stored as an icirc. Infinite side -> 0. The first side is the closest to the origin.
vcors
	Vertex coordinates, stored as a t_COMPLEX with real components. The side sides[i] has vertices vcor[i-1] and vcor[i], appearing in this order going counterclockwise about the origin.
vargs
	The argument of the corresponding vertex, normalized to lie in [0, 2*Pi), stored as t_REAL. These are stored in counterclockwise order, BUT vargs itself is not sorted: it is increasing from 1 to crossind, and from crossind+1 to the end.
crossind
	the arc from vertex crossind to crossind+1 contains the positive x-axis. If one vertex IS on this axis, then crossind+1 gives that vertex.	Alternatively, crossind is the unique vertex such that vargs[crossind]>vargs[crossind+1], taken cyclically.
kact
	The action of the corresponding element on the Klein model, for use in klein_act. Infinite side -> 0
area
	The hyperbolic area of U, which will be oo unless we are a finite index subgroup of Gamma.
spair
	Stores the side pairing of U, if it exists/has been computed. When computing the normalized boundary, this will store the indices of the sides that got deleted.
infinite
	Stores the indices of the infinite sides.
*/


/*Initializes the inputs for normbound_icircs. G is the set of elements we are forming the normalized boundary for. Returns NULL if no elements giving an isometric circle are input. Not gerepileupto safe, and leaves garbage.*/
static GEN
normbound(GEN X, GEN G, GEN (*Xtoklein)(GEN, GEN), GEN gdat)
{
  long lG = lg(G), i;
  GEN C = vectrunc_init(lG);
  GEN indtransfer = vecsmalltrunc_init(lG);/*Transfer indices of C back to G*/
  for (i = 1; i < lG; i++) {
	GEN circ = icirc_elt(X, gel(G, i), Xtoklein, gdat);
	if(gequal0(gel(circ, 3))) continue;/*No isometric circle*/
	vectrunc_append(C, circ);
	vecsmalltrunc_append(indtransfer, i);
  }
  if(lg(C) == 1) return NULL;
  return normbound_icircs(C, indtransfer, gdat);
}

/*Given C, the output of icirc_elt for each of our elements, this computes and returns the normalized boundary. Assumptions:
	C[i]=[g, M, icirc], where g=elt, M=Kleinian action, and icirc=[a, b, r, p1, p2, ang1, ang2].
	None of the entries should be infinite sides, sort this out in the method that calls this.
	C has at least one element, so the boundary is non-trivial.
	indtransfer tracks the indices of G in terms of the indices in C, for tracking the deleted elements.
	Not gerepileupto safe, and leaves garbage.
*/
static GEN
normbound_icircs(GEN C, GEN indtransfer, GEN gdat)
{
  GEN tol = gdat_get_tol(gdat);
  long prec=lg(tol);
  GEN order = gen_indexsort(C, NULL, &cmp_icircangle);/*Order the sides by initial angle.*/
  long lc = lg(C), maxsides = 2*lc, istart = normbound_icircs_bigr(C, order);/*Largest r value index (wrt order).*/
  GEN elts = cgetg(maxsides, t_VECSMALL);/*Stores indices in C of the sides. 0 represents an infinite side.*/
  GEN vcors = cgetg(maxsides, t_VEC);/*Stores coordinates of the vertices.*/
  GEN vargs = cgetg(maxsides, t_VEC);/*Stores arguments of the vertices.*/
  GEN infinite = vecsmalltrunc_init(lc);/*Stores the indices of infinite sides.*/
  GEN deleted = vecsmalltrunc_init(lc);/*Stores the indices in C of sides that don't survive.*/
  elts[1] = order[istart];
  GEN firstcirc = gmael(C, elts[1], 3);/*The equation for the fist line, used in phase 2 and for detecting phase 2 starting.*/
  gel(vcors, 1) = gel(firstcirc, 5);/*The terminal vertex of the side is the first vertex.*/
  gel(vargs, 1) = gel(firstcirc, 7);
  long found = 1, lenc = lc-1, absind;/*found=how many sides we have found up to now. This can increase and decrease.*/
  int phase2 = 0;/*Which phase we are in*/
  /*PHASE 1: inserting sides, where the end vertex does not intersect the first side. All we need to update / keep track of are:
	elts, vcors, vargs, infinite, deleted, found, absind, phase2
	PHASE 2: The same, but we have looped back around, and need to insert the last edge into the first (which will never disappear).
  */
  for (absind = 1; absind < lenc; absind++) {/*We are trying to insert the next side.*/
	long toins = order[1+(istart+absind-1)%lenc];/*The index of the side to insert.*/
	GEN curcirc = gmael(C, toins, 3);/*The isometric circle we are trying to insert*/
	GEN lastcirc = gmael(C, elts[found], 3);/*The isometric circle of the last side inserted. This is NEVER an oo side.*/
	int termloc = angle_onarc(gel(lastcirc, 6), gel(lastcirc, 7), gel(curcirc, 7), tol);/*Whether the terminal angle lies inside the last arc.*/
	switch (termloc) {
	  case 3:
	    if (angle_onarc(gel(lastcirc, 6), gel(lastcirc, 7), gel(curcirc, 6), tol)) {
		  vecsmalltrunc_append(deleted, indtransfer[toins]);
		  continue;/*The initial angle also lies here, so we are totally enveloped, and we move on.*/
		}
	  case 1:
	    phase2 = 1;/*We have looped back around and are intersecting from the right. This is also valid when case=3 and we didn't continue.*/
		/*The last side might not be the first, but since we looped all the way around there MUST be an oo side.*/
		normbound_icircs_insinfinite(elts, vcors, vargs, infinite, curcirc, &found);/*Insert oo side!*/
		normbound_icircs_phase2(elts, vcors, vargs, curcirc, firstcirc, tol, toins, &found);/*Phase 2 insertion.*/
		continue;
	  case 2:
	    vecsmalltrunc_append(deleted, indtransfer[toins]);
	    continue;/*The terminal angle is the same as the previous, so we are enveloped. It is not possible for our new side to envelop the old side.*/
	}
	/*If we make it here, then the terminal point lies outside of the last range.*/
	int initloc = angle_onarc(gel(lastcirc, 6), gel(lastcirc, 7), gel(curcirc, 6), tol);/*Whether the initial angle lies inside the last arc.*/
	switch (initloc) {
	  case 3:/*We have an intersection to consider, will do so outside the switch.*/
	    break;
	  case 1:
	    vecsmalltrunc_append(deleted, indtransfer[elts[found]]);
		absind--;/*We completely envelop the previous side. We need to delete it and redo this case.*/
		found--;
		continue;
	  case 0:
	    normbound_icircs_insinfinite(elts, vcors, vargs, infinite, curcirc, &found);/*Insert oo side!*/
	  case 2:
	    if (phase2 || angle_onarc(gel(curcirc, 6), gel(curcirc, 7), gel(firstcirc, 6), tol)) {/*Phase2 has started*/
		  phase2 = 1;
		  normbound_icircs_phase2(elts, vcors, vargs, curcirc, firstcirc, tol, toins, &found);/*Phase 2 insertion.*/
		  continue;
		}
	    normbound_icircs_insclean(elts, vcors, vargs, curcirc, toins, &found);/*New(initial)=Old(terminal), so there is no infinite side coming first (or we are coming from case=0 when we have already inserted it)*/
		continue;
	}
	/*If we make it here, our current circle intersects the last one, so we need to see if it is "better" than the previous intersection.*/
	GEN ipt = line_int11(curcirc, lastcirc, tol);/*Find the intersection point, guaranteed to be in the unit disc.*/
	GEN iptarg = argmod_complex(ipt, tol);/*Argument*/
	if (found == 1) {/*Straight up insert it; no phase 2 guaranteed.*/
	   gel(vcors, found) = ipt;/*Fix the last vertex*/
	   gel(vargs, found) = iptarg;/*Fix the last vertex argument*/
	   normbound_icircs_insclean(elts, vcors, vargs, curcirc, toins, &found);/*Clean insert*/
	   continue;
	}
	int newptdat = angle_onarc(gel(vargs, found-1), gel(vargs, found), iptarg, tol);/*The new point wrt the previous side.*/
	switch (newptdat) {
	  case 3:/*Insert it!*/
	    gel(vcors, found) = ipt;/*Fix the last vertex*/
	    gel(vargs, found) = iptarg;/*Fix the last vertex argument*/
		if (phase2 || angle_onarc(gel(curcirc, 6), gel(curcirc, 7), gel(firstcirc, 6), tol)) {/*Phase2 has started*/
		  phase2 = 1;
		  normbound_icircs_phase2(elts, vcors, vargs, curcirc, firstcirc, tol, toins, &found);/*Phase 2 insertion.*/
		  continue;
		}
	    normbound_icircs_insclean(elts, vcors, vargs, curcirc, toins, &found);/*Clean insert*/
	    continue;
	  case 2:
	    if (phase2) {
		  vecsmalltrunc_append(deleted, indtransfer[toins]);
		  continue;/*Previous vertex was intersection with firstcirc, which cannot be deleted. Thus we are enveloped.*/
		}
	    /*Now there is no need to fix previous vertex, and we don't start in phase 2*/
		if (angle_onarc(gel(curcirc, 6), gel(curcirc, 7), gel(firstcirc, 6), tol)) {/*Phase2 has started, but we can insert our side*/
		  phase2 = 1;
		  normbound_icircs_phase2(elts, vcors, vargs, curcirc, firstcirc, tol, toins, &found);/*Phase 2 insertion.*/
		  continue;
		}
		normbound_icircs_insclean(elts, vcors, vargs, curcirc, toins, &found);/*Clean insert*/
		continue;
	  case 0:
	    if (!phase2) {/*We supercede the last side, so delete and try again. We could not have deleted a side, have this trigger, and intersect the last side without superceding.*/
		  vecsmalltrunc_append(deleted, indtransfer[elts[found]]);
	      absind--;
	      found--;
	      continue;
		}
		int supercede = angle_onarc(gel(curcirc, 6), gel(curcirc, 7), gel(vargs, found - 1), tol);/*We are in phase 2, and we either miss the side to the left (so ignore and move on), or to the right (supercede previous edge).*/
		if (supercede) {
		  vecsmalltrunc_append(deleted, indtransfer[elts[found]]);
		  absind--;
		  found--;
		  continue;
		}
		vecsmalltrunc_append(deleted, indtransfer[toins]);
		continue;
	  case 1:/*We just replace the last side, we have the same vertex and vertex angle.*/
	    vecsmalltrunc_append(deleted, indtransfer[elts[found]]);
		found--;
		if(phase2 || angle_onarc(gel(curcirc, 6), gel(curcirc, 7), gel(firstcirc, 6), tol)) {/*Phase2 has started*/
		  normbound_icircs_phase2(elts, vcors, vargs, curcirc, firstcirc, tol, toins, &found);
		  continue;
		}
		normbound_icircs_insclean(elts, vcors, vargs, curcirc, toins, &found);/*Clean insert*/
	}
  }
  if (!phase2) {/*We never hit phase 2, so there is one final infinite edge to worry about.*/
	normbound_icircs_insinfinite(elts, vcors, vargs, infinite, firstcirc, &found);
  }
  /*Now we can compile everything into the return vector.*/
  long i, fp1 = found + 1;
  GEN rv = cgetg(10, t_VEC);
  GEN rv_elts = cgetg(fp1, t_VEC);
  GEN rv_sides = cgetg(fp1, t_VEC);
  GEN rv_kact = cgetg(fp1, t_VEC);
  for (i = 1; i <= found; i++) {/*Make sure we treat infinite sides correctly!*/
	gel(rv_elts, i) = elts[i] ? gmael(C, elts[i], 1) : gen_0;
	gel(rv_sides, i) = elts[i] ? gmael(C, elts[i], 3) : gen_0;
	gel(rv_kact, i) = elts[i] ? gmael(C, elts[i], 2) : gen_0;
  }
  gel(rv, 1) = rv_elts;/*The elements*/
  gel(rv, 2) = rv_sides;/*The sides*/
  gel(rv, 3) = vec_shorten(vcors, found);/*Vertex coords*/
  gel(rv, 4) = vec_shorten(vargs, found);/*Vertex arguments*/
  gel(rv, 5) = stoi(args_find_cross(gel(rv, 4)));/*Crossing point*/
  gel(rv, 6) = rv_kact;/*Kleinian action*/
  if (lg(infinite) > 1) gel(rv, 7) = mkoo();/*Infinite side means infinite area.*/
  else gel(rv, 7) = normbound_area(rv_sides, prec);
  gel(rv, 8) = deleted;/*For now, stores the deleted indices.*/
  gel(rv, 9) = infinite;/*Infinite sides*/
  return rv;
}

/*Finds the index i such that C[order[i]] has the largest r value, which is guaranteed to be on the normalized boundary.*/
static long
normbound_icircs_bigr(GEN C, GEN order)
{
  pari_sp av = avma;
  long i, best = 1;
  GEN maxr = gmael3(C, order[1], 3, 3);
  for (i = 2; i < lg(C); i++) {
	GEN newr = gmael3(C, order[i], 3, 3);
	if (cmprr(newr, maxr) <= 0) continue;/*Not as big. No need for tolerance here, it won't affect anything.*/
	best = i;
	maxr = newr;
  }
  return gc_long(av, best);
}

/*curcirc gives a new side that does not loop back around, but there is an infinite side first. This inserts the infinite side.*/
static void
normbound_icircs_insinfinite(GEN elts, GEN vcors, GEN vargs, GEN infinite, GEN curcirc, long *found)
{
  (*found)++;
  long f = *found;
  elts[f] = 0;/*Infinite side*/
  vecsmalltrunc_append(infinite, f);/*This cannot get deleted, so will be always accurate.*/
  gel(vcors, f) = gel(curcirc, 4);/*Initial vertex of curcirc.*/
  gel(vargs, f) = gel(curcirc, 6);/*Initial angle of curcirc.*/
}

/*We are inserting a new side that does not intersect back into phase 2. If it intersected with the previous side, then that vertex should be updated before calling this.*/
static void
normbound_icircs_insclean(GEN elts, GEN vcors, GEN vargs, GEN curcirc, long toins, long *found)
{
  (*found)++;
  long f = *found;
  elts[f] = toins;/*Non-infinite side we are inserting.*/
  gel(vcors, f) = gel(curcirc, 5);/*Terminal vertex of curcirc.*/
  gel(vargs, f) = gel(curcirc, 7);/*Terminal vertex of curcirc argument.*/
}

/*We are performing an insertion in phase 2, i.e. we are intersecting back with the initial side. Assume we have already dealt with updating the previous vertex, if applicable.*/
static void
normbound_icircs_phase2(GEN elts, GEN vcors, GEN vargs, GEN curcirc, GEN firstcirc, GEN tol, long toins, long *found)
{
  (*found)++;
  long f = *found;
  elts[f] = toins;
  gel(vcors, f) = line_int11(curcirc, firstcirc, tol);/*Intersect with initial side*/
  gel(vargs, f) = argmod_complex(gel(vcors, f), tol);/*Argument*/
}

/*Returns the hyperbolic area of the normalized boundary, which is assumed to not have any infinite sides (we keep track if they exist, and do not call this method if they do). C should be the list of [a, b, r] in order. The area is (n-2)*Pi-sum(angles), where there are n sides.*/
static GEN
normbound_area(GEN C, long prec)
{
  pari_sp av = avma;
  long n = lg(C) - 1, i;
  GEN area = mulsr(n - 2, mppi(prec));/*(n-2)*Pi*/
  for (i = 1; i < n; i++) area = subrr(area, icirc_angle(gel(C, i), gel(C, i + 1), prec));
  area = subrr(area, icirc_angle(gel(C, n), gel(C, 1), prec));
  return gerepileupto(av, area);
}

/*Let ind be the index of the edge that z is on when projected from the origin to the boundary (2 possibilities if it is a vertex). Returns 0 if z is in the interior of U, -ind if z is on the boundary, and ind if z is outside the boundary. Assume that the normalized boundary is non-trivial.*/
static long
normbound_outside(GEN U, GEN z, GEN tol)
{
  pari_sp av = avma;
  GEN arg = argmod_complex(z, tol);
  long sideind = args_search(normbound_get_vargs(U), normbound_get_cross(U), arg, tol);/*Find the side*/
  if (sideind < 0) sideind = -sideind;/*We line up with a vertex.*/
  GEN side = gel(normbound_get_sides(U), sideind);/*The side!*/
  int which = normbound_whichside(side, z, tol);
  if (which == -1) return gc_long(av, 0);/*Same side as 0, so inside.*/
  if (which == 0) return gc_long(av, -sideind);/*On the side*/
  return gc_long(av, sideind);/*Outside*/
}

/*Given an isometric circle side, this returns -1 if z is on the same side as the origin, 0 if it is on it, and 1 else*/
static int
normbound_whichside(GEN side, GEN z, GEN tol)
{
  pari_sp av = avma;
  if (gequal0(side)) {
	if (toleq(normcr(z), gen_1, tol)) return gc_int(av, 0);/*Unit disc point, infinite side*/
	return gc_int(av, -1);/*Interior point, infinite side*/
  }
  GEN where = subrs(addrr(mulrr(gel(z, 1), gel(side, 1)), mulrr(gel(z, 2), gel(side, 2))), 1);/*sign(where)==-1 iff where is on the same side as 0, 0 iff where is on the side.*/
  return gc_int(av, tolsigne(where, tol));
}

/*2: NORMALIZED BOUNDARY APPENDING*/

/*Initializes the inputs for normbound_append_icircs. G is the set of elements we are forming the normalized boundary for, and U is the necessarily non-trivial current boundary. Returns NULL if the normalized boundary does not change. Not gerepileupto safe, and leaves garbage.*/
static GEN
normbound_append(GEN X, GEN U, GEN G, GEN (*Xtoklein)(GEN, GEN), GEN gdat)
{
  pari_sp av = avma;
  long lG = lg(G), lU = lg(gel(U, 1)), i;
  GEN Cnew = vectrunc_init(lG);/*Start by initializing the new circles.*/
  GEN indtransfer = vecsmalltrunc_init(lG);/*Transfer indices of C back to G*/
  for (i = 1; i < lG; i++) {
	GEN circ = icirc_elt(X, gel(G, i), Xtoklein, gdat);
	if(gequal0(gel(circ, 3))) continue;/*No isometric circle*/
	vectrunc_append(Cnew, circ);
	vecsmalltrunc_append(indtransfer, i);
  }
  long lCnew = lg(Cnew);
  if (lCnew == 1) return gc_NULL(av);
  long lenU = lU-1;/*Now we figure out the order of U, before merging them.*/
  long Ustart = normbound_get_cross(U)%lenU+1;/*The side which crosses the origin. We do a linear probe to find the smallest initial angle.*/
  long Unext = Ustart%lenU+1;
  GEN Usides = normbound_get_sides(U), tol = gdat_get_tol(gdat);
  for (;;) {
	if (gequal0(gel(Usides, Ustart))) {Ustart = Unext; break;}/*The next side is not infinite, and must be the first one.*/
	if (gequal0(gel(Usides, Unext))) {Ustart = Unext%lenU + 1; break;}/*Must move past the infinite side.*/
	if (tolcmp(gmael(Usides, Ustart, 6), gmael(Usides, Unext, 6), tol) >= 0) {Ustart = Unext; break;}/*We have crossed over.*/
	Ustart = Unext;
	Unext = Unext%lenU + 1;
  }
  GEN Uelts = normbound_get_elts(U), Ukact = normbound_get_kact(U);
  GEN order = gen_indexsort(Cnew, NULL, &cmp_icircangle);/*Order the new sides by initial angle.*/
  long Cbestr = normbound_icircs_bigr(Cnew, order);/*Biggest r value of the new sides*/
  long rbigind;/*Stores the largest r value index, which is either U[1], or Cnew[Cbestr].*/
  if (tolcmp(gmael3(Cnew, order[Cbestr], 3, 3), gmael(Usides, 1, 3), tol) > 0) rbigind = -Cbestr;/*Negative to indicate C*/
  else rbigind = 0;/*it's 1, but we store it as 0, and swap it to the index that U[1] appears in when found.*/
  GEN C = vectrunc_init(lCnew + lU);/*Tracks the sorted isometric circle triples.*/
  GEN Ctype = vecsmalltrunc_init(lCnew + lU);/*Tracks how the C-index corresponds to U or G. negative=new and positive=old.*/
  long Cnewind, Uind;
  if (tolcmp(gmael3(Cnew, order[1], 3, 6), gmael(Usides, Ustart, 6), tol) >= 0) {/*Start with U.*/
    vectrunc_append(C, mkvec3(gel(Uelts, Ustart), gel(Ukact, Ustart), gel(Usides, Ustart)));
	vecsmalltrunc_append(Ctype, Ustart);/*Old side*/
	if (!rbigind && Ustart == 1) rbigind = 1;
	Cnewind = 1;
	Uind = Ustart%lenU + 1;
	if (gequal0(gel(Usides, Uind))) Uind = Uind%lenU + 1;/*Move past infinite side.*/
  }
  else {/*Start with C*/
	vectrunc_append(C, gel(Cnew, order[1]));
	vecsmalltrunc_append(Ctype, -indtransfer[order[1]]);/*New side*/
	if (rbigind == -1) rbigind = 1;
	Cnewind = 2;
	Uind = Ustart;
  }
  for (;;) {/*Now we merge the two into one list, keeping track of old vs new sides.*/
	if (lg(C)>Cnewind && Uind == Ustart) {/*We have finished with the U's (first inequality tracks if we have added any U's).*/
	  while (Cnewind < lCnew) {/*Just C's left.*/
	    vectrunc_append(C, gel(Cnew, order[Cnewind]));
	    vecsmalltrunc_append(Ctype, -indtransfer[order[Cnewind]]);
		if (rbigind == -Cnewind) rbigind = lg(C) - 1;
	    Cnewind++;
	  }
	  break;
	}
	if (Cnewind == lCnew) {/*Done with C's but not the U's.*/
	  do {
		vectrunc_append(C, mkvec3(gel(Uelts, Uind), gel(Ukact, Uind), gel(Usides, Uind)));
		vecsmalltrunc_append(Ctype, Uind);
		if (!rbigind && Uind == 1) rbigind = lg(C) - 1;
		Uind = Uind%lenU + 1;
	    if (gequal0(gel(Usides, Uind))) Uind = Uind%lenU + 1;/*Move past infinite side.*/
	  } while (Uind != Ustart);
	  break;
	}
	/*We still have U's and C's left.*/
	if (cmprr(gmael3(Cnew, order[Cnewind], 3, 6), gmael(Usides, Uind, 6)) >= 0) {/*No need for tolcmp, we don't care if =.*/
	  vectrunc_append(C, mkvec3(gel(Uelts, Uind), gel(Ukact, Uind), gel(Usides, Uind)));
	  vecsmalltrunc_append(Ctype, Uind);
	  if (!rbigind && Uind == 1) rbigind = lg(C) - 1;
	  Uind = Uind%lenU + 1;
	  if (gequal0(gel(Usides, Uind))) Uind = Uind%lenU + 1;/*Move past infinite side.*/
	  continue;
	}
	vectrunc_append(C, gel(Cnew, order[Cnewind]));
	vecsmalltrunc_append(Ctype, -indtransfer[order[Cnewind]]);
	if (rbigind == -Cnewind) rbigind = lg(C) - 1;
	Cnewind++;
  }
  return normbound_append_icircs(normbound_get_vcors(U), normbound_get_vargs(U), C, Ctype, rbigind, gdat);
}

/*Does normbound_icircs, except we already have a normalized boundary U that we add to. We prep this method with normbound_append. This method is very similar to normbound_icircs.
Ucors: the previously computed coordinates of vertices
Uargs: the previuosly computed arguments of vertices
C: Same as for normbound
Ctype[i] = ind>0 if C[i] = U[1][ind], and -ind if C[i] = G[-ind] where G was input in normbound_append.
If we do NOT change the boundary, we just return NULL. This is for convenience with normbasis.
Not gerepileupto safe, and leaves garbage.
*/
static GEN
normbound_append_icircs(GEN Uvcors, GEN Uvargs, GEN C, GEN Ctype, long rbigind, GEN gdat)
{
  pari_sp av = avma;
  GEN tol = gdat_get_tol(gdat);
  long prec=lg(tol);
  long lc = lg(C), maxsides = 2*lc, lenU = lg(Uvcors) - 1;
  GEN elts = cgetg(maxsides, t_VECSMALL);/*Stores indices in C of the sides. 0 represents an infinite side.*/
  GEN vcors = cgetg(maxsides, t_VEC);/*Stores coordinates of the vertices.*/
  GEN vargs = cgetg(maxsides, t_VEC);/*Stores arguments of the vertices.*/
  GEN infinite = vecsmalltrunc_init(lc);/*Stores the indices of infinite sides.*/
  GEN deleted = vecsmalltrunc_init(lc);/*Stores the indices in C of sides that don't survive.*/
  elts[1] = rbigind;
  GEN firstcirc = gmael(C, elts[1], 3);/*The equation for the fist line, used in phase 2 and for detecting phase 2 starting.*/
  gel(vcors, 1) = gel(firstcirc, 5);/*The terminal vertex of the side is the first vertex.*/
  gel(vargs, 1) = gel(firstcirc, 7);
  long found = 1, lenc = lc - 1, absind;/*found=how many sides we have found up to now. This can increase and decrease.*/
  int phase2 = 0;/*Which phase we are in, and if there are infinite sides or not*/
  long newU = 0;/*The first index of a new side, used at the end to see if we are the same as before or not.*/
  if (Ctype[rbigind] < 0) newU = 1;/*We use this to detect if the normalized boundary changed or not.*/
  /*PHASE 1: inserting sides, where the end vertex does not intersect the first side. All we need to update / keep track of are:
	elts, vcors, vargs, infinite, deleted, found, absind, phase2, infinitesides, newU
	PHASE 2: The same, but we have looped back around, and need to insert the last edge into the first (which will never disappear).
  */
  for (absind = 1; absind < lenc; absind++) {/*We are trying to insert the next side.*/
	long toins = 1 + (rbigind + absind - 1)%lenc;/*The index of the side to insert.*/
	GEN curcirc = gmael(C, toins, 3);/*The isometric circle we are trying to insert*/
	GEN lastcirc = gmael(C, elts[found], 3);/*The isometric circle of the last side inserted. This is NEVER an oo side.*/
	long t1 = Ctype[toins], t2 = Ctype[elts[found]];/*Whether the sides are new/old, and the corresponding indices if old.*/
	GEN ipt, iptarg;
	if (!phase2 && t1 > 0 && t2 > 0) {/*Two old sides, and not in phase 2!!*/
	  long diff = smodss(t1 - t2, lenU);
	  if (diff > 1) {/*There was an infinite side here, and it remains!*/
		normbound_icircs_insinfinite(elts, vcors, vargs, infinite, curcirc, &found);/*Insert oo side!*/
		if (phase2 || angle_onarc(gel(curcirc, 6), gel(curcirc, 7), gel(firstcirc, 6), tol)) {/*Phase 2 has started*/
		  phase2 = 1;
		  normbound_icircs_phase2(elts, vcors, vargs, curcirc, firstcirc, tol, toins, &found);/*Phase 2 insertion.*/
		  continue;
		}
		normbound_icircs_insclean(elts, vcors, vargs, curcirc, toins, &found);/*Clean insert, no phase 2.*/
		continue;
	  }
	  /*Now our old side intersected the last old side. We can recover the intersection point and angle. This part is very similar to when we have an intersection of the last and current sides as seen below.*/
	  ipt = gel(Uvcors, t2);
	  iptarg = gel(Uvargs, t2);
	  if (found == 1) {/*Only 2 sides: this and the old one, they intersect, and we can't be in phase 2.*/
	    gel(vcors, found) = ipt;/*Fix the last vertex*/
	    gel(vargs, found) = iptarg;/*Fix the last vertex argument*/
		normbound_icircs_insclean(elts, vcors, vargs, curcirc, toins, &found);/*Clean insert*/
		continue;
	  }
	  long efm1 = elts[found - 1];
	  if (!efm1 || Ctype[efm1] > 0) {/*Previous vertex was infinite/old, so we have the exact same structure and can just insert it.*/
	    gel(vcors, found) = ipt;/*Fix the last vertex*/
	    gel(vargs, found) = iptarg;/*Fix the last vertex argument*/
		if (phase2 || angle_onarc(gel(curcirc, 6), gel(curcirc, 7), gel(firstcirc, 6), tol)) {/*Phase2 has started*/
		  phase2 = 1;
		  normbound_icircs_phase2(elts, vcors, vargs, curcirc, firstcirc, tol, toins, &found);/*Phase 2 insertion.*/
		  continue;
		}
	    normbound_icircs_insclean(elts, vcors, vargs, curcirc, toins, &found);/*Clean insert*/
	    continue;
	  }
	  /*Previous vertex was new, so we might envelop the last side. The rest is the same as normal, so we use a GOTO statement.*/
	  goto GENERICINT;
	}
	/*One of the two is new, so we do the normal thing.*/
	int termloc = angle_onarc(gel(lastcirc, 6), gel(lastcirc, 7), gel(curcirc, 7), tol);/*Whether the terminal angle lies inside the last arc.*/
	switch (termloc) {
	  case 3:
	    if (angle_onarc(gel(lastcirc, 6), gel(lastcirc, 7), gel(curcirc, 6), tol)) {
		  vecsmalltrunc_append(deleted, t1);
		  continue;/*The initial angle also lies here, so we are totally enveloped, and we move on.*/
		}
	  case 1:
	    phase2 = 1;/*We have looped back around and are intersecting from the right. This is also valid when case=3 and we didn't continue.*/
		normbound_icircs_insinfinite(elts, vcors, vargs, infinite, curcirc, &found);/*Insert oo side!*/
		normbound_icircs_phase2(elts, vcors, vargs, curcirc, firstcirc, tol, toins, &found);/*Phase 2 insertion.*/
		if (!newU && t1 < 0) newU = found;/*This is the first new circle.*/
		continue;
	  case 2:
	    vecsmalltrunc_append(deleted, t1);
	    continue;/*The terminal angle is the same as the previous, so we are enveloped. It is not possible for our new side to envelop the old side.*/
	}
	/*If we make it here, then the terminal point lies outside of the last range.*/
	int initloc = angle_onarc(gel(lastcirc, 6), gel(lastcirc, 7), gel(curcirc, 6), tol);/*Whether the initial angle lies inside the last arc.*/
	switch (initloc) {
	  case 3:/*We have an intersection to consider, will do so outside the switch.*/
	    break;
	  case 1:
	    if (newU == found) newU = 0;/*The last side was the only new one, and we are deleting it.*/
		vecsmalltrunc_append(deleted, Ctype[elts[found]]);
		absind--;/*We completely envelop the previous side. We need to delete it and redo this case.*/
		found--;
		continue;
	  case 0:
	    normbound_icircs_insinfinite(elts, vcors, vargs, infinite, curcirc, &found);/*Insert oo side!*/
	  case 2:
	    if (phase2 || angle_onarc(gel(curcirc, 6), gel(curcirc, 7), gel(firstcirc, 6), tol)) {/*Phase2 has started*/
		  phase2 = 1;
		  normbound_icircs_phase2(elts, vcors, vargs, curcirc, firstcirc, tol, toins, &found);/*Phase 2 insertion.*/
		  if (!newU && t1 < 0) newU = found;/*This is the first new circle.*/
		  continue;
		}
	    normbound_icircs_insclean(elts, vcors, vargs, curcirc, toins, &found);/*New(initial)=Old(terminal), so there is no infinite side coming first (or we are coming from case=0 when we have already inserted it)*/
		if (!newU && t1 < 0) newU = found;/*This is the first new circle.*/
		continue;
	}
	/*If we make it here, our current circle intersects the last one, so we need to see if it is "better" than the previous intersection.*/
	ipt = line_int11(curcirc, lastcirc, tol);/*Find the intersection point, guaranteed to be in the unit disc.*/
	iptarg = argmod_complex(ipt, tol);/*Argument*/
	if (found == 1) {/*Straight up insert it; no phase 2 guaranteed.*/
	   gel(vcors, found) = ipt;/*Fix the last vertex*/
	   gel(vargs, found) = iptarg;/*Fix the last vertex argument*/
	   normbound_icircs_insclean(elts, vcors, vargs, curcirc, toins, &found);/*Clean insert*/
	   if (!newU && t1 < 0) newU = found;/*This is the first new circle.*/
	   continue;
	}
	GENERICINT:;
	int newptdat = angle_onarc(gel(vargs, found-1), gel(vargs, found), iptarg, tol);/*The new point wrt the previous side.*/
	switch (newptdat) {
	  case 3:/*Insert it!*/
	    gel(vcors, found) = ipt;/*Fix the last vertex*/
	    gel(vargs, found) = iptarg;/*Fix the last vertex argument*/
		if (phase2 || angle_onarc(gel(curcirc, 6), gel(curcirc, 7), gel(firstcirc, 6), tol)) {/*Phase2 has started*/
		  phase2 = 1;
		  normbound_icircs_phase2(elts, vcors, vargs, curcirc, firstcirc, tol, toins, &found);/*Phase 2 insertion.*/
		  if (!newU && t1 < 0) newU = found;/*This is the first new circle.*/
		  continue;
		}
	    normbound_icircs_insclean(elts, vcors, vargs, curcirc, toins, &found);/*Clean insert*/
		if (!newU && t1 < 0) newU = found;/*This is the first new circle.*/
	    continue;
	  case 2:
	    if (phase2) {
		  vecsmalltrunc_append(deleted, t1);
		  continue;/*Previous vertex was intersection with firstcirc, which cannot be deleted. Thus we are enveloped.*/
		}
	    /*Now there is no need to fix previous vertex, and we don't start in phase 2*/
		if (angle_onarc(gel(curcirc, 6), gel(curcirc, 7), gel(firstcirc, 6), tol)) {/*Phase2 has started, but we can insert our side*/
		  phase2 = 1;
		  normbound_icircs_phase2(elts, vcors, vargs, curcirc, firstcirc, tol, toins, &found);/*Phase 2 insertion.*/
		  if (!newU && t1 < 0) newU = found;/*This is the first new circle.*/
		  continue;
		}
		normbound_icircs_insclean(elts, vcors, vargs, curcirc, toins, &found);/*Clean insert*/
		if (!newU && t1 < 0) newU = found;/*This is the first new circle.*/
		continue;
	  case 0:
	    if (!phase2) {/*We supercede the last side, so delete and try again. We could not have deleted a side, have this trigger, and intersect the last side without superceding.*/
		  if (newU == found) newU = 0;/*The last side was the only new one, and we are deleting it.*/
		  vecsmalltrunc_append(deleted, Ctype[elts[found]]);
	      absind--;
	      found--;
	      continue;
		}
		int supercede = angle_onarc(gel(curcirc, 6), gel(curcirc, 7), gel(vargs, found-1), tol);/*We are in phase 2, and we either miss the side to the left (so ignore and move on), or to the right (supercede previous edge).*/
		if (supercede) {
		  if (newU == found) newU = 0;/*The last side was the only new one, and we are deleting it.*/
		  vecsmalltrunc_append(deleted, Ctype[elts[found]]);
		  absind--;
		  found--;
		  continue;
		}
		vecsmalltrunc_append(deleted, t1);
		continue;
	  case 1:/*We just replace the last side, we have the same vertex and vertex angle.*/
	    if (newU == found) newU = 0;/*The last side was the only new one, and we are deleting it.*/
		vecsmalltrunc_append(deleted, Ctype[elts[found]]);
		if (!newU && t1 < 0) newU = found;/*This is the first new circle.*/
		found--;
		if(phase2 || angle_onarc(gel(curcirc, 6), gel(curcirc, 7), gel(firstcirc, 6), tol)) {/*Phase2 has started*/
		  normbound_icircs_phase2(elts, vcors, vargs, curcirc, firstcirc, tol, toins, &found);
		  continue;
		}
		normbound_icircs_insclean(elts, vcors, vargs, curcirc, toins, &found);/*Clean insert*/
	}
  }
  if (!newU) return gc_NULL(av);/*We added no new sides, so return NULL.*/
  if (!phase2) {/*We never hit phase 2, so there is one final infinite edge to worry about.*/
	normbound_icircs_insinfinite(elts, vcors, vargs, infinite, firstcirc, &found);
  }
  /*Now we can compile everything into the return vector.*/
  long i, fp1 = found + 1;
  GEN rv = cgetg(10, t_VEC);
  GEN rv_elts = cgetg(fp1, t_VEC);
  GEN rv_sides = cgetg(fp1, t_VEC);
  GEN rv_kact = cgetg(fp1, t_VEC);
  for (i = 1; i <= found; i++) {/*Make sure we treat infinite sides correctly!*/
	gel(rv_elts, i) = elts[i] ? gmael(C, elts[i], 1) : gen_0;
	gel(rv_sides, i) = elts[i] ? gmael(C, elts[i], 3) : gen_0;
	gel(rv_kact, i) = elts[i] ? gmael(C, elts[i], 2) : gen_0;
  }
  gel(rv, 1) = rv_elts;/*The elements*/
  gel(rv, 2) = rv_sides;/*The sides*/
  gel(rv, 3) = vec_shorten(vcors, found);/*Vertex coords*/
  gel(rv, 4) = vec_shorten(vargs, found);/*Vertex arguments*/
  gel(rv, 5) = stoi(args_find_cross(gel(rv, 4)));/*Crossing point*/
  gel(rv, 6) = rv_kact;/*Kleinian action*/
  if (lg(infinite) > 1) gel(rv, 7) = mkoo();/*Infinite side means infinite area.*/
  else gel(rv, 7) = normbound_area(rv_sides, prec);
  gel(rv, 8) = deleted;/*For now, stores the deleted indices.*/
  gel(rv, 9) = infinite;/*Infinite sides*/
  return rv;
}


/*2: NORMALIZED BOUNDARY ANGLES*/

/*Used for sorting C by initial angles.*/
static int
cmp_icircangle(void *nul, GEN c1, GEN c2){return cmprr(gmael(c1, 3, 6), gmael(c2, 3, 6));}

/*Assuming there is a unique index i such that args[i]>args[i+1], we return i. No need for tolerance. We assume that #args>=2.*/
static long
args_find_cross(GEN args)
{
  long l1 = 1, l2 = lg(args) - 1;
  int c = cmprr(gel(args, l1), gel(args, l2));
  if (c < 0) return l2;/*Sorted already, so it only loops back at the end.*/
  while (l2-l1 > 1) {
	long l = (l1 + l2)>>1;/*floor((l1+l2)/2)*/
	int c = cmprr(gel(args, l1), gel(args, l));
	if (c < 0) l1 = l;
	else l2 = l;
  }
  return l1;
}

/*args is a partially sorted list of angles in [0, 2*Pi): ind is the unique index with args[ind]>args[ind+1]. If arg==args[i] up to tolerance, we return -i. Otherwise, we return the index where it should be inserted. If this is at the end, we return 1, since it could be inserted there equally well. We assume that #args>=2.*/
static long
args_search(GEN args, long ind, GEN arg, GEN tol)
{
  long na = lg(args)-1, l1, l2;/*l1 and l2 track the min and max indices possible.*/
  int c1 = tolcmp(gel(args, 1), arg, tol);/*Compare to entry 1.*/
  if (c1 == 0) return -1;/*Equal*/
  if (c1 < 0) {/*In block from 1 to ind*/
    l1 = 1;
	int cind = tolcmp(gel(args, ind), arg, tol);
    if (cind < 0) return (ind%na) + 1;/*Make sure we overflow back to 1 if ind=na. This also covers if ind=1.*/
    if (cind == 0) return -ind;
	l2 = ind;/*We are strictly between args[l1] and args[l2], and l1<l2.*/
  }
  else {/*In block from ind+1 to na.*/
	if (ind == na) return 1;/*Strictly increasing, so we can insert at the start.*/
	l1 = ind + 1;
	int cind = tolcmp(gel(args, l1), arg, tol);
	if (cind > 0) return l1;
	if (cind == 0) return -l1;
	int cend = tolcmp(gel(args, na), arg, tol);
	if (cend < 0) return 1;/*We can insert it at the start. This also covers if l1=na.*/
	if (cend == 0) return -na;
	l2 = na;/*We are strictly between args[l1] and args[l2], and l1<l2.*/
  }
  while (l2 - l1 > 1) {
	long l = (l1 + l2)>>1;/*floor((l1+l2)/2)*/
	int c = tolcmp(gel(args, l), arg, tol);
	if (c > 0) {l2 = l;continue;}
	if (c == 0) return -l;
	l1 = l;
  }
  return l2;
}


/*2: NORMALIZED BASIS*/

/*Returns the edge pairing, as VECSMALL v where v[i]=j means i is paired with j (infinite sides are by convention paired with themselves). If not all sides can be paired, instead returns [v1, v2, ...] where vi=[gind, vind, side, location] is a VECSMALL. gind is the index of the unpaired side, vind is the corresponding index of an unpaired vertex (so gind==vind or gind==vind+1 mod # of sides). The point gv is on side side, and location=-1 means it is inside U, 0 on the boundary, and 1 is outside U. If v is infinite, then location<=0 necessarily.*/
static GEN
edgepairing(GEN U, GEN tol)
{
  pari_sp av = avma;
  GEN vcors = normbound_get_vcors(U);/*Vertex coordinates*/
  GEN vargs = normbound_get_vargs(U);/*Vertex arguments*/
  long cross = normbound_get_cross(U);/*crossing point, used in args_search*/
  GEN kact = normbound_get_kact(U);/*Action of the sides*/
  GEN sides = normbound_get_sides(U);/*Sides*/
  long lU = lg(kact), lenU = lU - 1, i;/*Number of sides*/
  GEN unpair = vectrunc_init(2*lU - 1), pair = zero_zv(lenU);/*Unpair stores the unpaired edges, pair stores the paired edges*/
  for (i = 1; i < lU; i++) {/*Try to pair the ith side.*/
    GEN act = gel(kact, i);
	if (gequal0(act)) {pair[i] = i;continue;}/*oo side, go next (we say it is paired with itself)*/
	if (pair[i]) continue;/*We have already paired this side, so move on.*/
	GEN v1 = klein_act_i(act, gel(vcors, i));/*Image of the ith vertex.*/
	GEN v1arg = argmod_complex(v1, tol);/*Argument*/
	long v1pair = args_search(vargs, cross, v1arg, tol), v1loc, foundv1;
	if (v1pair < 0) {/*We found the angle. In theory the points might not be equal, so we check this now. It suffices to check if v1 is on the side.*/
	  v1pair = -v1pair;
	  v1loc = normbound_whichside(gel(sides, v1pair), v1, tol);
	  if (v1loc) foundv1 = 0;/*Sadly, they are not equal.*/
	  else foundv1 = 1;
	}
	else {/*We didn't find it.*/
	  foundv1 = 0;
	  v1loc = normbound_whichside(gel(sides, v1pair), v1, tol);
	}
	long i2;
	if (i == 1) i2 = lenU;
	else i2 = i-1;/*i2 is the previous vertex, also on the ith side.*/
	GEN v2 = klein_act_i(act, gel(vcors, i2));/*Image of the second (i-1 st) vertex.*/
	GEN v2arg = argmod_complex(v2, tol);/*Argument*/
	long v2pair = args_search(vargs, cross, v2arg, tol), v2loc, foundv2;
	if (v2pair < 0) {/*We found the angle. In theory the points might not be equal, so we check this now. It suffices to check if v2 is on the side.*/
	  v2pair = -v2pair;
	  v2loc = normbound_whichside(gel(sides, v2pair), v2, tol);
	  if (v2loc) foundv2 = 0;/*Sadly, they are not equal.*/
	  else foundv2 = 1;
	}
	else {/*We didn't find it.*/
	  foundv2 = 0;
	  v2loc = normbound_whichside(gel(sides, v2pair), v2, tol);
	}
	if (foundv1) {/*v1 paired*/
	  if (foundv2) {/*Both vertices were paired!*/
		pair[i] = v2pair;
		pair[v2pair] = i;
		continue;
	  }
	  vectrunc_append(unpair, mkvecsmall4(i, i2, v2pair, v2loc));/*Add the data for v2.*/
	  continue;
	}
	/*v1 is unpaired*/
	vectrunc_append(unpair, mkvecsmall4(i, i, v1pair, v1loc));/*Add the data for v1.*/
	if (!foundv2) vectrunc_append(unpair, mkvecsmall4(i, i2, v2pair, v2loc));/*Add the data for v2.*/
  }
  if (lg(unpair) == 1) return gerepilecopy(av, pair);/*All sides paired*/
  return gerepilecopy(av, unpair);
}

/*Returns the normalized basis of G. Can pass in U as a normalized basis to append to, or NULL to just start with G. If no normalized basis, returns NULL.*/
static GEN
normbasis(GEN X, GEN U, GEN G, GEN (*Xtoklein)(GEN, GEN), GEN (*Xmul)(GEN, GEN, GEN), GEN (*Xinv)(GEN, GEN), int (*Xistriv)(GEN, GEN), GEN gdat)
{
  pari_sp av = avma;
  GEN tol = gdat_get_tol(gdat);
  long prec = lg(tol);
  GEN origin = gtocr(gen_0, prec);
  long lG = lg(G), i, j = lG;
  GEN Gstart = cgetg(2*lG - 1, t_VEC), Gadd;/*Make the vector with G and its inverses.*/
  for (i = 1; i < lG; i++) {
	gel(Gstart, i) = gel(G, i);
	gel(Gstart, j) = Xinv(X, gel(G, i));
	j++;
  }
  /*Start by initializing U and Gadd by reducing the elements of Gstart*/
  if (!U) {/*No U is passed in.*/
	U = normbound(X, Gstart, Xtoklein, gdat);
	if (!U) return gc_NULL(av);/*No valid isometric circles.*/
	GEN del = normbound_get_spair(U), g;/*The deleted sides*/
	long ldel = lg(del);
	Gadd = vectrunc_init(ldel);
	for (i = 1; i<ldel; i++) {
	  g = red_elt(X, U, Xinv(X, gel(Gstart, del[i])), origin, Xtoklein, Xmul, 0, gdat);/*Reduce x^-1 for each deleted side x.*/
	  if (Xistriv(X, g)) continue;/*Reduces to 1, move on.*/
	  vectrunc_append(Gadd, Xinv(X, g));/*Add g^-1 to our list*/
	}
  }
  else {/*U is passed in, so just initialize Gadd by reducing each element.*/
    lG = lg(Gstart);
	Gadd = vectrunc_init(lG);
	GEN g;
	for (i = 1; i < lG; i++) {
	  g = red_elt(X, U, gel(Gstart, i), origin, Xtoklein, Xmul, 0, gdat);/*Reduce x. No need to do the inverse, as our list includes all inverses of elements so we just hit it later.*/
	  if (Xistriv(X, g)) continue;/*Reduces to 1, move on.*/
	  vectrunc_append(Gadd, Xinv(X, g));/*Add g^-1 to our list*/
	}
  }
  for (;;) {/*We have passed to U, and must see if it has changed or not. Gadd is the set of elements we last tried to add.*/
    while (lg(Gadd) > 1) {/*Continue reducing elements and recomputing the boundary until stable.*/
	  GEN Unew = normbound_append(X, U, Gadd, Xtoklein, gdat);
	  if (!Unew) {/*The boundary did not change, so we must check every element of Gadd against U.*/
		lG = lg(Gadd);
		GEN Gnew = vectrunc_init(lG), g;
		for (i = 1; i < lG; i++) {
		  g = red_elt(X, U, Xinv(X, gel(Gadd, i)), origin, Xtoklein, Xmul, 0, gdat);/*Reduce g.*/
		  if (Xistriv(X, g)) continue;/*Reduces to 1, move on.*/
		  vectrunc_append(Gnew, Xinv(X, g));/*Add g^-1 to our list*/
		}
		Gadd = Gnew;
		continue;/*Try again!*/
	  }
	  GEN del = normbound_get_spair(Unew);/*The deleted sides*/
	  long ldel = lg(del), ind;
	  GEN Gnew = vectrunc_init(ldel), Uelts = normbound_get_elts(U), g;
	  for (i = 1; i<ldel; i++) {
		ind = del[i];
		if (ind > 0) g = red_elt(X, Unew, Xinv(X, gel(Uelts, ind)), origin, Xtoklein, Xmul, 0, gdat);/*Reduce element from U*/
		else g = red_elt(X, Unew, Xinv(X, gel(Gadd, -ind)), origin, Xtoklein, Xmul, 0, gdat);/*Reduce element from Gadd*/
	    if (Xistriv(X, g)) continue;/*Reduces to 1, move on.*/
	    vectrunc_append(Gnew, Xinv(X, g));/*Add g^-1 to our list*/
	  }
	  Gadd = Gnew;
	  U = Unew;/*Update U*/
	}
	/*At this point, <G> = <U[1]>_(no-inverse). Now we have to check if there is a side pairing.*/
	GEN unpair = edgepairing(U, tol);
	if (typ(unpair) == t_VECSMALL) {/*All paired up!*/
	  gel(U, 8) = unpair;
	  return gerepilecopy(av, U);
	}
	GEN vcors = normbound_get_vcors(U);/*Vertex coordinates.*/
	GEN kact = normbound_get_kact(U);
	GEN elts = normbound_get_elts(U);
	long lunp = lg(unpair), lenU = lg(elts)-1;
	Gadd = vectrunc_init(lunp);/*For each entry [gind, vind, gv side, location] in unpair, we will find a new element intersecting inside U.*/
	for (i = 1; i < lunp; i++) {
	  GEN dat = gel(unpair, i);
	  GEN v = gel(vcors, dat[2]);
	  long gind = dat[1];/*The index of g*/
	  if (!toleq(normcr(v), gen_1, tol)) {/*Interior vertex.*/
	    if (dat[4] == -1) {/*gv is INSIDE U, so since gv is on I(g^-1), we add in g^-1.*/
		  vectrunc_append(Gadd, Xinv(X, gel(elts, gind)));
		  continue;
		}
	    if (dat[4] == 0) {/*gv is on another side. We need to see if it is I(g^-1) or not.*/
		  long gvind = dat[3];/*The side that gv is on.*/
		  GEN endptimg = klein_act_i(gel(kact, gind), gel(vcors, smodss(gind - 2, lenU) + 1));/*Find the image of the first vertex of the side gind.*/
		  if (toleq(endptimg, gel(vcors, gvind), tol)) {/*gv is on I(g^-1) which is part of U.*/
			/*Let the element of the side after gdind be w (so v is the intersection). Then I(wg^-1) contains gv, and will give gv as an intersection vertex, so we add in wg^-1.*/
			vectrunc_append(Gadd, Xmul(X, gel(elts, gind%lenU + 1), Xinv(X, gel(elts, gind))));
		  }
		  else {/*gv is on a side that is NOT I(g^-1), so adding in g^-1 gives us a new side of U.*/
			vectrunc_append(Gadd, Xinv(X, gel(elts, gind)));
		  }
		  continue;
		}
	    /*Now dat[4]=1, i.e. gv is outside of U. We reduce g with respect to v to get Wg, where Wgv is strictly closer to the origin than v is. Therefore I(Wg) encloses v, so will delete the vertex v and change U. This is the case described by John Voight.*/
		vectrunc_append(Gadd, red_elt(X, U, stoi(gind), v, Xtoklein, Xmul, 0, gdat));/*Add the reduction in.*/
		continue;
	  }
	  /*We have an infinite vertex.*/
	  if (gequal0(gel(elts, dat[3]))) {/*gv is on an infinite side, so g^-1 gives a new side of U.*/
		vectrunc_append(Gadd, Xinv(X, gel(elts, gind)));
		continue;
	  }
	  /*gv is strictly inside I(w). Any interior point v' on I(g) close enough to v has gv' inside I(w) as well, so wgv' is closer to the origin, hence v' is inside I(wg). Thus we can add in wg to Gadd and envelop the vertex v.*/
	  vectrunc_append(Gadd, Xmul(X, gel(elts, dat[3]), gel(elts, gind)));
	}
  }
}



/*2: NORMALIZED BOUNDARY REDUCTION*/

/*Reduces g with respect to z, i.e. finds g' such that g'(gz) is inside U (or on the boundary), and returns [g'g, decomp], where decomp is the Vecsmall of indices so that g'=algmulvec(A, U[1], decomp).*/
static GEN
red_elt_decomp(GEN X, GEN U, GEN g, GEN z, GEN (*Xtoklein)(GEN, GEN), GEN (*Xmul)(GEN, GEN, GEN), GEN gdat)
{
  pari_sp av = avma;
  GEN tol = gdat_get_tol(gdat);
  GEN zorig = z;
  z = klein_act_i(Xtoklein(X, g), z);/*Starting point*/
  GEN elts = normbound_get_elts(U);
  GEN kact = normbound_get_kact(U);
  long ind = 1, maxind = 32;
  GEN decomp = cgetg(maxind + 1, t_VECSMALL);
  for (;;) {
	if (ind%20 == 0) {/*Recompute z every 20 moves to account for precision issues.*/
	  GEN gact = Xtoklein(X, g);
	  z = klein_act_i(gact, zorig);
	}
	long outside = normbound_outside(U, z, tol);
	if (outside <= 0) {/*We are inside or on the boundary.*/
	  GEN gact = Xtoklein(X, g);
	  z = klein_act_i(gact, zorig);
	  outside = normbound_outside(U, z, tol);
	  if (outside <= 0) break;/*We recompile to combat loss of precision, which CAN and WILL happen.*/
	}
    z = klein_act_i(gel(kact, outside), z);/*Act on z.*/
    g = Xmul(X, gel(elts, outside), g);/*Multiply on the left of g.*/
	if (ind > maxind) {/*Make decomp longer*/
	  maxind = maxind<<1;/*Double it.*/
	  decomp = vecsmall_lengthen(decomp, maxind);
	}
	decomp[ind] = outside;
	ind++;
  }
  decomp = vecsmall_shorten(decomp, ind - 1);/*Shorten it. We also need to reverse it, since decomp goes in the reverse order to what it needs to be.*/
  return gerepilecopy(av, mkvec2(g, vecsmall_reverse(decomp)));
}

/*red_elt_decomp, except we return g'g if flag=0, g'gz if flag=1, and [g'g,g'gz] if flag=2. Can pass g as t_INT so that g=U[1][g].*/
static GEN
red_elt(GEN X, GEN U, GEN g, GEN z, GEN (*Xtoklein)(GEN, GEN), GEN (*Xmul)(GEN, GEN, GEN), int flag, GEN gdat)
{
  pari_sp av = avma;
  GEN tol = gdat_get_tol(gdat);
  GEN elts = normbound_get_elts(U);
  GEN kact = normbound_get_kact(U);
  GEN zorig = z;
  GEN gklein;
  if (typ(g) == t_INT) {/*Retrieve it.*/
	long gind = itos(g);
	gklein = gel(kact, gind);
	g = gel(elts, gind);
  }
  else gklein = Xtoklein(X, g);/*Compute it.*/
  z = klein_act_i(gklein, z);/*Starting point.*/
  long ind = 1;
  for (;;) {
	if (ind%20 == 0) {/*Recompute z every 20 moves to account for precision issues.*/
	  GEN gact = Xtoklein(X, g);
	  z = klein_act_i(gact, zorig);
	}
	long outside = normbound_outside(U, z, tol);
	if (outside <= 0) {/*We are inside or on the boundary.*/
	  GEN gact = Xtoklein(X, g);
	  z = klein_act_i(gact, zorig);
	  outside = normbound_outside(U, z, tol);
	  if (outside <= 0) break;/*We recompile to combat loss of precision, which CAN and WILL happen.*/
	}
    z = klein_act_i(gel(kact, outside), z);/*Act on z.*/
    g = Xmul(X, gel(elts, outside), g);/*Multiply on the left of g.*/
	ind++;
  }
  if (flag == 0) return gerepilecopy(av, g);
  if (flag == 1) return gerepilecopy(av, z);
  return gerepilecopy(av, mkvec2(g, z));
}


/*SECTION 3: QUATERNION ALGEBRA METHODS*/


/*ALGEBRA REQUIREMENTS: UPDATE THIS!
Let
	F be a totally real number field
	A a quaternion algebra over F split at a unique real place
	O be an order in A
We can compute the fundamental domain of groups that live between O^1 and N{B^{\times}}(O)^+, i.e. berween the units of norm 1, and the elements of the normalizer with positive norm at the unique split real place.
Inputs to most methods are named "X", which represents the algebra A, the order O, data to describe the exact group we are computing, and various other pieces of data that will be useful. You should first initialize this with afuchinit.
*/


/*3: INITIALIZE ARITHMETIC FUCHSIAN GROUPS*/


/*ARITHMETIC FUCHSIAN GROUPS FORMATTING
An arithmetic Fuchsian group initialization is represented by X. In general, elements of the algebra A are represented in terms of the basis of O, NOT in terms of the basis of A. If O is the identity, then these notions are identical. We have:
X
	[A, [O, Oinv], Oconj, Omultable, chol, embmats, type, gdat, fdom, pres, sig]
A
	The algebra
O, Oinv
	The order and its inverse, given as a matrix whose columns generate the order (with respect to the stored order in A).
Oconj
	Conjugates of the elements of O. To conjugate an element x of A, it suffices to compute Oconj*x.
Omultable
	To multiply elements of O more quickly.
chol
	Cholesky decomposition of the norm form on O, used in algnorm_chol to compute norms quickly.
embmats
	O[,i] is sent to embmats[i] which acts on the unit disc/Klein model.
type
	Which symmetric space we want to compute.
gdat
	Geometric data.
fdom
	Fundamental domain, if computed.
pres
	Presentation, if computed.
sig
	Signature, if computed.
*/

/*
TO DO: CHOL, TYPE
type=0 means O^1, the default.*/

/*Initializes the arithmetic Fuchsian group of the given inputs, ready to compute a fundamental domain. Not gerepileupto suitable, and leaves garbage.*/
static GEN
afuchinit_i(GEN A, GEN O, GEN type, GEN p, long prec)
{
  if(!O) O = matid(lg(alg_get_basis(A)) - 1);
  GEN gdat = gdat_initialize(p, prec);
  GEN AX = cgetg(12, t_VEC);
  gel(AX, 1) = A;
  gel(AX, 2) = mkvec2(O, QM_inv(O));
  long lgO = lg(O), i;
  gel(AX, 3) = cgetg(lgO, t_MAT);
  for (i = 1; i < lgO; i++) gmael(AX, 3, i) = algconj(A, gel(O, i));
  gel(AX, 4) = Omultable(A, O);
  gel(AX, 5) = gen_0;/*TO DO*/
  gel(AX, 6) = afuch_make_kleinmats(A, O, gdat_get_p(gdat), prec);/*Make sure p is safe.*/
  gel(AX, 7) = type;/*TO DO*/
  gel(AX, 8) = gdat;
  gel(AX, 9) = gen_0;
  gel(AX, 10) = gen_0;
  gel(AX, 11) = gen_0;
  return AX;
}

/*Clean initialization of the symmetric space. Can pass p as NULL and will set it to the default.*/
GEN
afuchinit(GEN A, GEN O, GEN type, GEN p, long prec)
{
  pari_sp av = avma;
  if(!p) p = defp(prec);
  return gerepilecopy(av, afuchinit_i(A, O, type, p, prec));
}

/*Returns a vector v of pairs [A, B] such that O[,i] is sent to v[i] acting on the Klein model..*/
static GEN
afuch_make_kleinmats(GEN A, GEN O, GEN p, long prec)
{
  pari_sp av = avma;
  GEN m2r = afuch_make_m2rmats(A, O, prec);
  long l, i;
  GEN klein = cgetg_copy(m2r, &l);
  for (i = 1; i < l; i++) {
	gel(klein, i) = m2r_to_klein(gel(m2r, i), p);
  }
  return gerepileupto(av, klein);
}

/*Returns a vector v of matrices in M(2, R) such that O[,i] is sent to v[i].*/
static GEN
afuch_make_m2rmats(GEN A, GEN O, long prec)
{
  pari_sp av = avma;
  long split = algsplitoo(A);/*The split real place*/
  if (split == 0) pari_err_TYPE("Quaternion algebra has 0 or >=2 split real infinite places, not valid for fundamental domains.", A);
  GEN K = alg_get_center(A);/*The centre, i.e K where A=(a,b/K)*/
  GEN Kroot = gel(nf_get_roots(K), split);/*The split root*/
  long Kvar = nf_get_varn(K);
  GEN a = alggeta(A), b = lift(alg_get_b(A));/*A=(a, B/K).*/
  GEN aval = poleval(a, Kroot);
  GEN bval = poleval(b, Kroot), rt;
  int apos;/*Tracks if a>0 at the split place or not, as this determines which embedding we will take.*/
  if (signe(aval) == 1) {apos = 1; rt = gsqrt(aval, prec);}
  else {apos = 0; rt = gsqrt(bval, prec);}
  long lO = lg(O), i, j;
  GEN mats = cgetg(lO, t_VEC);/*To store the 2x2 matrices.*/
  for (i = 1; i < lO; i++) {
	GEN x = algbasisto1ijk(A, gel(O, i));/*ith basis element in the form e+fi+gj+hk*/
	for (j=1; j<=4; j++) gel(x, j) = gsubst(gel(x, j), Kvar, Kroot);/*Evaluate it at Kroot. Will be real if K!=Q, else will be rational.*/
	if (apos) {/*e+fi+gj+hk -> [e+fsqrt(a), b(g+hsqrt(a));g-hsqrt(a), e-fsqrt(a)*/
	  GEN frta = gmul(gel(x, 2), rt);/*f*sqrt(a)*/
	  GEN hrta = gmul(gel(x, 4), rt);/*h*sqrt(a)*/
	  GEN topl = gtofp(gadd(gel(x, 1), frta), prec);/*e+f*sqrt(a)*/
	  GEN topr = gtofp(gmul(gadd(gel(x, 3), hrta), bval), prec);/*b(g+h*sqrt(a))*/
	  GEN botl = gtofp(gsub(gel(x, 3), hrta), prec);/*g-h*sqrt(a)*/
	  GEN botr = gtofp(gsub(gel(x, 1), frta), prec);/*e-f*sqrt(a)*/
	  gel(mats, i) = mkmat22(topl, topr, botl, botr);
	  continue;
	}
	/*e+fi+gj+hk -> [e+g*sqrt(b), a(f-h*sqrt(b));f+h*sqrt(b), e-g*sqrt(b)]*/
	GEN grtb = gmul(gel(x, 3), rt);/*g*sqrt(b)*/
	GEN hrtb = gmul(gel(x, 4), rt);/*h*sqrt(b)*/
	GEN topl = gtofp(gadd(gel(x, 1), grtb), prec);/*e+g*sqrt(b)*/
	GEN topr = gtofp(gmul(gsub(gel(x, 2), hrtb), aval), prec);/*a(f-h*sqrt(b))*/
	GEN botl = gtofp(gadd(gel(x, 2), hrtb), prec);/*f+h*sqrt(b)*/
	GEN botr = gtofp(gsub(gel(x, 1), grtb), prec);/*e-g*sqrt(b)*/
	gel(mats, i) = mkmat22(topl, topr, botl, botr);
  }
  return gerepilecopy(av, mats);
}

/*Initializes a table to compute multiplication of elements written in terms of O's basis more efficiently.*/
static GEN
Omultable(GEN A, GEN O)
{
  long n = lg(O), i, j;
  GEN M = cgetg(n, t_VEC);
  for (i = 1; i < n; i++) {
	gel(M, i) = cgetg(n, t_MAT);
	for (j = 1; j < n; j++) {
	  gmael(M, i, j) = algmul(A, gel(O, i), gel(O, j));
	}
  }
  return M;
}


/*3: ALGEBRA FUNDAMENTAL DOMAIN METHODS*/

/*Returns the isometric circle of an element of A.*/
GEN
afuchicirc(GEN X, GEN g)
{
  pari_sp av = avma;
  GEN gdat = afuch_get_gdat(X);
  GEN ginO = QM_QC_mul(afuch_get_orderinv(X), g);/*Convert into O's basis.*/
  GEN icirc_all = icirc_elt(X, ginO, &afuchtoklein, gdat);
  return gerepilecopy(av, gel(icirc_all, 3));
}

/*Returns the element representing the action of g on the Klein model.*/
GEN
afuchklein(GEN X, GEN g)
{
  pari_sp av = avma;
  GEN ginO = QM_QC_mul(afuch_get_orderinv(X), g);/*Convert into O's basis.*/
  return gerepileupto(av, afuchtoklein(X, ginO));
}

/*Returns the normalized boundary of the set of elements G in A. ASSUMES O=IDENTITY*/
GEN
afuchnormbasis(GEN X, GEN G)
{
  pari_sp av = avma;
  GEN gdat = afuch_get_gdat(X);
  GEN U = normbasis(X, NULL, G, &afuchtoklein, &afuchmul, &afuchinv, &afuchistriv, gdat);
  if (U) return gerepilecopy(av, U);
  return gc_const(av, gen_0);/*No normalized basis, return 0.*/
}

/*Returns the normalized boundary of the set of elements G in A. ASSUMES O=IDENTITY*/
GEN
afuchnormbound(GEN X, GEN G)
{
  pari_sp av = avma;
  GEN gdat = afuch_get_gdat(X);
  GEN U = normbound(X, G, &afuchtoklein, gdat);
  if (U) return gerepilecopy(av, U);
  return gc_const(av, gen_0);/*No normalized boundary, return 0.*/
}

/*Computes the normalized boundary of U union G, where U is an already computed normalized boundary. ASSUMES O=IDENTITY*/
GEN
afuchnormbound_append(GEN X, GEN U, GEN G)
{
  pari_sp av = avma;
  GEN gdat = afuch_get_gdat(X);
  GEN Unew = normbound_append(X, U, G, &afuchtoklein, gdat);
  if (Unew) return gerepilecopy(av, Unew);
  return gc_const(av, gen_0);
}

/*Reduces gz to the normalized boundary, returning [g'g, g'gz, decomp] where g'gz is inside the normalized boundary. Can supply g=NULL, which defaults it to the identity element. ASSUMES O=IDENTITY*/
GEN
afuchredelt(GEN X, GEN U, GEN g, GEN z)
{
  pari_sp av = avma;
  GEN gdat = afuch_get_gdat(X);
  GEN tol = gdat_get_tol(gdat);
  GEN zsafe = gtocr(z, lg(tol));/*Make sure z has the appropriate formatting.*/
  if(!g) g = afuchid(X);
  GEN r = red_elt_decomp(X, U, g, zsafe, &afuchtoklein, &afuchmul, gdat);/*[g'g, decomp].*/
  GEN gact = afuchtoklein(X, gel(r, 1));
  GEN zimg = klein_act_i(gact, zsafe);
  return gerepilecopy(av, mkvec3(gel(r, 1), zimg, gel(r, 2)));
}


/*3: ALGEBRA BASIC AUXILLARY METHODS*/

/*Returns the identity element*/
static GEN
afuchid(GEN X){return col_ei(lg(alg_get_tracebasis(afuch_get_alg(X)))-1, 1);}

/*alginv formatted for the input of an afuch, for use in the geometry section. Since we work in O, entries are all in Z.*/
static GEN
afuchinv(GEN X, GEN g){return ZM_ZC_mul(afuch_get_orderconj(X), g);}

/*Returns 1 if g is a scalar, i.e. g==conj(g).*/
static int
afuchistriv(GEN X, GEN g)
{
  pari_sp av = avma;
  GEN gconj = afuchinv(X, g);/*It's actually the conjugate, not the inverse.*/
  return gc_int(av, ZV_equal(g, gconj));
}

/*algmul formatted for the input of an afuch, for use in the geometry section.*/
static GEN
afuchmul(GEN X, GEN g1, GEN g2)
{
  pari_sp av = avma;
  GEN T = afuch_get_ordermultable(X);
  GEN g = ZC_Z_mul(ZM_ZC_mul(gel(T, 1), g2), gel(g1, 1));
  long lg1 = lg(g1), i;
  for (i = 2; i < lg1; i++) {
	g = ZC_add(g, ZC_Z_mul(ZM_ZC_mul(gel(T, i), g2), gel(g1, i)));
  }
  return gerepileupto(av, g);
}

/*Given an element g of O (of non-zero norm) written in basis form, this returns the image of g in acting on the Klein model.*/
static GEN
afuchtoklein(GEN X, GEN g)
{
  pari_sp av = avma;
  GEN embs = afuch_get_embmats(X);
  GEN A = mulcri(gmael(embs, 1, 1), gel(g, 1));
  GEN B = mulcri(gmael(embs, 1, 2), gel(g, 2));
  long lg = lg(g), i;
  for (i = 2; i < lg; i++) {
	A = addcrcr(A, mulcri(gmael(embs, i, 1), gel(g, i)));
	B = addcrcr(B, mulcri(gmael(embs, i, 2), gel(g, i)));
  }
  return gerepilecopy(av, mkvec2(A, B));
}


/*3: ALGEBRA HELPER METHODS*/

/*Given an element in the algebra representation of A, returns [e, f, g, h], where x=e+fi+gj+hk. e, f, g, h live in the centre of A.*/
GEN
algalgto1ijk(GEN A, GEN x)
{
  pari_sp av = avma;
  x = liftall(x);/*Make sure there are no mods.*/
  GEN L = alg_get_splittingfield(A);/*L=F(i)*/
  long Lvar = rnf_get_varn(L);
  GEN e = polcoef_i(gel(x, 1), 0, Lvar);
  GEN f = polcoef_i(gel(x, 1), 1, Lvar);
  GEN g = polcoef_i(gel(x, 2), 0, Lvar);
  GEN mh = polcoef_i(gel(x, 2), 1, Lvar);
  return gerepilecopy(av, mkvec4(e, f, g, gneg(mh)));
}

/*Given an element in the basis representation of A, returns [e, f, g, h], where x=e+fi+gj+hk. e, f, g, h live in the centre of A.*/
GEN
algbasisto1ijk(GEN A, GEN x)
{
  pari_sp av = avma;
  GEN xalg = algbasistoalg(A, x);
  return gerepileupto(av, algalgto1ijk(A, xalg));
}

/*Returns the conjugate of the element x in basis form. Not particularly fast.*/
static GEN
algconj(GEN A, GEN x)
{
  pari_sp av = avma;
  GEN tr = algtrace(A, x, 0);
  GEN trinA = algalgtobasis(A, mkcol2(tr, gen_0));/*Move it to A*/
  return gerepileupto(av, RgC_sub(trinA, x));
}

/*Given a quaternion algebra, return a,*/
static GEN
alggeta(GEN A)
{
  pari_sp av=avma;
  GEN L = alg_get_splittingfield(A);/*L=K(sqrt(a)).*/
  long Lvar = rnf_get_varn(L);/*Variable number for L*/
  return gerepileupto(av, gneg(polcoef_i(rnf_get_pol(L), Lvar, 0)));/*Defining polynomial is x^2-a, so retrieve a.*/
}

/*Returns G[L[1]]*G[L[2]]*...*G[L[n]], where L is a vecsmall or vec*/
GEN
algmulvec(GEN A, GEN G, GEN L)
{
  pari_sp av = avma;
  long n = lg(L), i;
  if (n == 1) return gerepilecopy(av, gel(alg_get_basis(A), 1));/*The identity*/
  GEN elts = cgetg(n, t_VEC);
  if (typ(L) == t_VECSMALL) {
    for (i = 1; i < n; i++) {
      long ind = L[i];
      if (ind > 0) gel(elts, i) = gel(G, ind);
      else gel(elts, i) = alginv(A, gel(G, -ind));
    }
  }
  else if (typ(L) == t_VEC) {
    for (i = 1; i < n; i++) {
      long ind = itos(gel(L, i));
      if (ind > 0) gel(elts, i) = gel(G, ind);
      else gel(elts, i) = alginv(A, gel(G, -ind));
    }
  }
  else pari_err_TYPE("L needs to be a vector or vecsmall of indices", L);
  return gerepileupto(av, gen_product(elts, &A, &voidalgmul));
}

/*Formats algmul for use in gen_product, which is used in algmulvec.*/
static GEN
voidalgmul(void *A, GEN x, GEN y){return algmul(*((GEN*)A), x, y);}

/*Returns the vector of finite ramified places of the algebra A.*/
GEN
algramifiedplacesf(GEN A)
{
  pari_sp av = avma;
  GEN hass = alg_get_hasse_f(A);/*Shallow*/
  long nhass = lg(gel(hass, 2)), i;
  GEN rp = vectrunc_init(nhass);
  for (i = 1; i < nhass; i++) {
    if (gel(hass, 2)[i] == 0) continue;/*Unramified*/
    vectrunc_append(rp, gmael(hass, 1, i));/*Ramified*/
  }
  return gerepilecopy(av, rp);
}

/*If the algebra A has a unique split infinite place, this returns the index of that place. Otherwise, returns 0.*/
static long
algsplitoo(GEN A)
{
  GEN infram = alg_get_hasse_i(A);/*shallow*/
  long split=0, linf = lg(infram), i;
  for (i = 1; i < linf; i++){/*Finding the split place*/
    if (infram[i] == 0){/*Split place*/
      if(split > 0) return 0;
      split = i;
    }
  }
  return split;/*No garbage!!*/
}


