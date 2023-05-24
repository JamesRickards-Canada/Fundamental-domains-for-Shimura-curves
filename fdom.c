/*TO DO
1. algorderdisc can be very slow in some cases. Maybe randomize the choice of i1 -> i4?
2. Do we want the debug level here to be the same as for algebras? Currently it is.
3. In afuch_moreprec, when alg_hilbert is updated to allow for denominators, this method can be simplified.
4. my_alg_changeorder may have a better successor with the updated quaternion algebra methods.
5. Don't check if in normalizer for primes dividing the discriminant or for unit norms.
6. Add testing for changing p.
7. When initialize by a, b allows for denominators, fix afuchlist to not use this.
8. Geodesics that intersect a vertex.

POSSIBLE FUTURE ADDITIONS:
1. Parallelization of element enumeration along with partial domain computations.
2. Methods of Imbert (see Voight's original paper) to find a minimal presentation in canonical form.
3. Computation of cohomology groups.
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

/*Possibly make it's own debugging level?*/
#define DEBUGLEVEL DEBUGLEVEL_alg

/*STATIC DECLARATIONS*/

/*SECTION 1: GEOMETRIC METHODS*/

/*1: LINES AND ARCS*/
static int angle_onarc(GEN a1, GEN a2, GEN a, GEN tol);
static GEN line_from_crrcrr(GEN z1, GEN z2, GEN tol);
static GEN line_int(GEN l1, GEN l2, GEN tol);
static GEN line_int11(GEN l1, GEN l2, GEN tol);
static GEN line_line_detNULL(GEN l1, GEN l2, GEN tol);

/*1: MATRIX ACTION ON GEOMETRY*/
static GEN defp(long prec);
static GEN gdat_initialize(GEN p, long prec);
static GEN klein_safe(GEN M, long prec);
static GEN uhp_safe(GEN p, long prec);

/*1: TRANSFER BETWEEN MODELS*/
static GEN m2r_to_klein(GEN M, GEN p);

/*1: DISTANCES/AREA*/
static GEN hdiscradius(GEN area, long prec);
static GEN hdiscrandom_arc(GEN R, GEN ang1, GEN ang2, long prec);

/*1: OPERATIONS ON COMPLEX REALS*/
static GEN gtocr(GEN z, long prec);
static GEN addcrcr(GEN z1, GEN z2);
static GEN divcrcr(GEN z1, GEN z2);
static GEN mulcrcr(GEN z1, GEN z2);
static GEN mulcrcr_conj(GEN z1, GEN z2);
static GEN mulcri(GEN z, GEN n);
static GEN mulcrr(GEN z, GEN r);
static GEN negc(GEN z);
static GEN normcr(GEN z);
static GEN rM_upper_r_mul(GEN M, GEN r);
static GEN rM_upper_add(GEN M1, GEN M2);
static GEN RgM_upper_add(GEN M1, GEN M2);
static GEN subcrcr(GEN z1, GEN z2);
static GEN subrcr(GEN r, GEN z);

/*1: TOLERANCE*/
static GEN deftol(long prec);
static GEN deflowtol(long prec);
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
static int cmp_icircangle(void *tol, GEN c1, GEN c2);
static long args_find_cross(GEN args);
static long args_search(GEN args, long ind, GEN arg, GEN tol);

/*2: NORMALIZED BASIS*/
static GEN edgepairing(GEN U, GEN tol);
static GEN normbasis(GEN X, GEN U, GEN G, GEN (*Xtoklein)(GEN, GEN), GEN (*Xmul)(GEN, GEN, GEN), GEN (*Xinv)(GEN, GEN), int (*Xistriv)(GEN, GEN), GEN gdat);

/*2: NORMALIZED BOUNDARY REDUCTION*/
static GEN red_elt_decomp(GEN X, GEN U, GEN g, GEN z, GEN (*Xtoklein)(GEN, GEN), GEN (*Xmul)(GEN, GEN, GEN), GEN gdat);
static GEN red_elt(GEN X, GEN U, GEN g, GEN z, GEN (*Xtoklein)(GEN, GEN), GEN (*Xmul)(GEN, GEN, GEN), int flag, GEN gdat);

/*2: CYCLES AND SIGNATURE*/
static GEN minimalcycles(GEN pair);
static GEN minimalcycles_bytype(GEN X, GEN U, GEN Xid, GEN (*Xmul)(GEN, GEN, GEN), int (*Xisparabolic)(GEN, GEN), int (*Xistriv)(GEN, GEN));
static GEN signature(GEN X, GEN U, GEN Xid, GEN (*Xmul)(GEN, GEN, GEN), int (*Xisparabolic)(GEN, GEN), int (*Xistriv)(GEN, GEN));

/*2: GEODESICS*/
static long fdom_intersect_sidesmidpt(long i1, long i2, long n);
static GEN fdom_intersect(GEN U, GEN geod, GEN tol, long s1);
static GEN geodesic_klein(GEN X, GEN g, GEN (*Xtoklein)(GEN, GEN), GEN tol);
static GEN geodesic_fdom(GEN X, GEN U, GEN g, GEN Xid, GEN (*Xtoklein)(GEN, GEN), GEN (*Xmul)(GEN, GEN, GEN), GEN (*Xinv)(GEN, GEN), GEN gdat);

/*2: PRESENTATION*/
static GEN presentation(GEN X, GEN U, GEN Xid, GEN (*Xmul)(GEN, GEN, GEN), int (*Xisparabolic)(GEN, GEN), int (*Xistriv)(GEN, GEN));
static void presentation_update(GEN words, long ind, GEN repl);
static GEN word_collapse(GEN word);
static GEN word_inv(GEN word);
static GEN word_substitute(GEN word, long ind, GEN repl, GEN invrepl);
static GEN word(GEN X, GEN U, GEN P, GEN g, GEN (*Xtoklein)(GEN, GEN), GEN (*Xmul)(GEN, GEN, GEN), GEN (*Xinv)(GEN, GEN), int (*Xistriv)(GEN, GEN), GEN gdat);

/*SECTION 3: QUATERNION ALGEBRA METHODS*/

/*3: INITIALIZE SYMMETRIC SPACE*/
static void afuch_moreprec(GEN X, long inc);
static GEN afuch_make_kleinmats(GEN A, GEN O, GEN p, long prec);
static GEN afuch_make_m2rmats(GEN A, GEN O, long prec);
static GEN afuch_make_qfmats(GEN kleinmats);
static GEN afuch_make_traceqf(GEN X, GEN nm, GEN Onorm);
static GEN Omultable(GEN A, GEN O, GEN Oinv);
static GEN Onorm_makechol(GEN F, GEN Onorm);
static GEN Onorm_makemat(GEN A, GEN O, GEN AOconj);
static GEN Onorm_toreal(GEN A, GEN Onorm);

/*3: ALGEBRA FUNDAMENTAL DOMAIN CONSTANTS*/
static GEN afuchO1area(GEN A, GEN O, GEN Olevel_fact, long computeprec, long prec);
static GEN afuchbestC(GEN A, GEN O, GEN Olevel_nofact, long prec);
static GEN afuchfdomdat_init(GEN A, GEN O, long prec);

/*3: ALGEBRA FUNDAMENTAL DOMAIN METHODS*/
static GEN afuchfdom_i(GEN X, GEN *startingset);
static int nextsub(GEN S, long n);

/*3: NON NORM 1 METHODS*/
static GEN afuch_makeunitelts(GEN X);
static GEN afuch_makeALelts(GEN X);
static GEN afuch_makenormelts(GEN X);
static GEN AL_make_norms(GEN B, long split, GEN ideals, long prec);
static GEN normalizer_make_norms(GEN B, long split, GEN ideals, long prec);

/*3: ALGEBRA BASIC AUXILLARY METHODS*/
static GEN afuchconj(GEN X, GEN g);
static GEN afuchid(GEN X);
static int afuchinnormalizer(GEN X, GEN g);
static int afuchisparabolic(GEN X, GEN g);
static int afuchistriv(GEN X, GEN g);
static GEN afuchmul(GEN X, GEN g1, GEN g2);
static GEN afuchnorm_fast(GEN X, GEN g);
static GEN afuchnorm_chol(GEN F, GEN chol, GEN g);
static GEN afuchnorm_mat(GEN F, GEN Onorm, GEN g);
static GEN afuchnorm_real(GEN X, GEN g);
static GEN afuchtoklein(GEN X, GEN g);
static GEN afuchtrace(GEN X, GEN g);

/*3: FINDING ELEMENTS*/
static GEN afuch_make_qf(GEN X, GEN nm, GEN z, GEN tracepart, GEN realnm);
static GEN afuchfindelts(GEN X, GEN nm, GEN z, GEN C, long maxelts, GEN tracepart, GEN realnm);
static GEN afuchfindoneelt_i(GEN X, GEN nm, GEN C);

/*3: ALGEBRA HELPER METHODS*/
static GEN algconj(GEN A, GEN x);
static GEN alggeta(GEN A);
static GEN voidalgmul(void *A, GEN x, GEN y);
static GEN algramifiedplacesf(GEN A);
static long algsplitoo(GEN A);

/*3: ALGEBRA ORDER METHODS*/
static GEN algd(GEN A, GEN a);
static GEN algorderdisc(GEN A, GEN O, int reduced, int factored);

/*SECTION 4: FINCKE POHST FOR FLAG=2 WITH PRUNING*/

/*4: SUPPORTING METHODS TO FINCKE POHST*/
static GEN clonefill(GEN S, long s, long t);
static int mpgreaterthan(GEN x, GEN y);
static GEN norm_aux(GEN xk, GEN yk, GEN zk, GEN vk);

/*4: MAIN FINCKE POHST METHODS*/
static GEN smallvectors_prune(GEN q, GEN C, GEN prune);

/*SECTION 5: CHANGING ORDER METHODS, WHICH WERE DELETED FROM LIBPARI*/
static GEN my_alg_changeorder(GEN al, GEN ord);
static GEN elementabsmultable(GEN mt, GEN x);
static GEN elementabsmultable_Fp(GEN mt, GEN x, GEN p);
static GEN algbasismultable(GEN al, GEN x);
static GEN algtracebasis(GEN al);
static GEN elementabsmultable_Z(GEN mt, GEN x);
static GEN FpM_trace(GEN x, GEN p);
static GEN ZM_trace(GEN x);

/*MAIN BODY*/


/*SECTION 1: GEOMETRIC METHODS*/

/*LINES, SEGMENTS, TOLERANCE, POINTS
Line    [a, b, c]
    Representing ax+by=c. We will normalize so that c=gen_1 or gen_0. It is assumed that at least one of a, b is non-zero. We also assume that a and b are of type t_REAL.
Point
    Stored as t_COMPLEX with t_REAL entries.
Segment [a, b, c, x0, x1]
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

/*Given two complex numbers with real components, this returns the line [a, b, c] through them (signifying ax+by=c), where c=0 or 1.*/
static GEN
line_from_crrcrr(GEN z1, GEN z2, GEN tol)
{
  pari_sp av = avma;
  GEN det = line_line_detNULL(z1, z2, tol);/*Still works for our inputs!*/
  if (!det) {/*det 0, hence 0*/
    return gerepilecopy(av, mkvec3(gel(z1, 2), negr(gel(z1, 1)), gen_0));
  }
  GEN scale = invr(det);
  GEN a = mulrr(subrr(gel(z2, 2), gel(z1, 2)), scale);
  GEN b = mulrr(subrr(gel(z1, 1), gel(z2, 1)), scale);
  return gerepilecopy(av, mkvec3(a, b, gen_1));
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

/*Given two lines given by a1x+b1y=c1 and a2x+b2y=c2, returns a1b2-a2b1, unless this is within tolerance of 0, when we return NULL. Gerepileupto safe, leaves garbage.*/
static GEN
line_line_detNULL(GEN l1, GEN l2, GEN tol)
{
  GEN d = subrr(mulrr(gel(l1, 1), gel(l2, 2)), mulrr(gel(l1, 2), gel(l2, 1)));/*ad-bc*/
  if (toleq0(d, tol)) return NULL;/*d=0 up to tolerance*/
  return d;
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
GEN
disc_to_klein(GEN z)
{
  pari_sp av = avma;
  GEN znorm = normcr(z);
  GEN scale = divsr(2, addsr(1, znorm));//2/(1+|z|^2)
  return gerepileupto(av, mulcrr(z, scale));/*2z/(1+|z|^2)*/
}

/*Given a point z in the unit disc model, this transfers it to the upper half plane model. The formula is conj(p)*z-p/(z-1)*/
GEN
disc_to_plane(GEN z, GEN p)
{
  pari_sp av = avma;
  GEN num = subcrcr(mulcrcr_conj(z, p), p);/*z*conj(p)-p*/
  GEN denom = mkcomplex(subrs(gel(z, 1), 1), gel(z, 2));/*z-1*/
  return gerepileupto(av, divcrcr(num, denom));
}

/*Given a point z in the Klein model, this transfers it to the unit disc model. The formula is z/(1+sqrt(1-|z|^2)).*/
GEN
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
GEN
klein_to_plane(GEN z, GEN p, GEN tol)
{
  pari_sp av = avma;
  GEN zdisc = klein_to_disc(z, tol);
  return gerepileupto(av, disc_to_plane(zdisc, p));//Klein -> disc -> plane
}

/*Given a point z in the upper half plane model, this transfers it to the unit disc model. The formula is (z-p)/(z-conj(p))*/
GEN
plane_to_disc(GEN z, GEN p)
{
  pari_sp av = avma;
  GEN num = subcrcr(z, p);/*z-p*/
  GEN denom = subcrcr(z, conj_i(p));/*z-conj(p)*/
  return gerepileupto(av, divcrcr(num, denom));
}

/*Given a point z in the upper half plane model, this transfers it to the Klein model.*/
GEN
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

/*Given the area (t_REAL) of a hyperbolic disc, this returns the radius. The formula is area=4*Pi*sinh(R/2)^2, or R=2arcsinh(sqrt(area/4Pi))*/
static GEN
hdiscradius(GEN area, long prec)
{
  pari_sp av = avma;
  GEN aover4pi = divrr(area, Pi2n(2, prec));/*area/4Pi*/
  GEN rt = sqrtr(aover4pi);
  GEN half = gasinh(rt, prec);
  return gerepileupto(av, shiftr(half, 1));
}

/*Returns a random point z in the Klein model, uniform inside the ball of radius R. See page 19 of Page (before section 2.5).*/
GEN
hdiscrandom(GEN R, long prec)
{
  pari_sp av = avma;
  GEN zbound = expIPiR(shiftr(randomr(prec), 1), prec);/*Random boundary point. Now we need to scale by a random hyperbolic distance in [0, R]*/
  GEN x = mulrr(gsinh(gshift(R, -1), prec), sqrtr(randomr(prec)));
  /*We need to return zbound*tanh(2asinh(x))=zbound*2x*sqrt(x^2+1)/(2*x^2+1).*/
  GEN xsqr = sqrr(x);
  GEN num = shiftr(mulrr(x, sqrtr(addrs(xsqr, 1))), 1);/*2*x*sqrt(x^2+1)*/
  GEN denom = addrs(shiftr(xsqr, 1), 1);/*2x^2+1*/
  return gerepileupto(av, gmul(zbound, divrr(num, denom)));
}

/*Returns a random point z in the unit disc, uniform inside the ball of radius R, with argument uniform in [ang1, ang2]. See page 19 of Page (before section 2.5).*/
static GEN
hdiscrandom_arc(GEN R, GEN ang1, GEN ang2, long prec)
{
  pari_sp av = avma;
  GEN arg = gadd(ang1, gmul(randomr(prec), gsub(ang2, ang1)));/*Random angle in [ang1, ang2]*/
  GEN zbound = expIr(arg);/*The boundary point. Now we need to scale by a random hyperbolic distance in [0, R]*/
  GEN x = mulrr(gsinh(gshift(R, -1), prec), sqrtr(randomr(prec)));
  /*We need to return tanh(2asinh(x))=2x*sqrt(x^2+1)/(2*x^2+1).*/
  GEN xsqr = sqrr(x);
  GEN num = shiftr(mulrr(x, sqrtr(addrs(xsqr, 1))), 1);/*2*x*sqrt(x^2+1)*/
  GEN denom = addrs(shiftr(xsqr, 1), 1);/*2x^2+1*/
  return gerepileupto(av, gmul(zbound, divrr(num, denom)));
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

/*Multiplies a square upper triangular matrix M with upper triangle having real entries by a real number r. Clean method.*/
static GEN
rM_upper_r_mul(GEN M, GEN r)
{
  long n = lg(M), i, j;
  GEN P = cgetg(n, t_MAT);
  for (i = 1; i < n; i++) {
    gel(P, i) = cgetg(n, t_COL);
    for (j = 1; j <= i; j++) gmael(P, i, j) = mulrr(gcoeff(M, j, i), r);
    for (j = i + 1; j < n; j++) gmael(P, i, j) = gen_0;
  }
  return P;
}

/*Adds two square upper triangular matrices M1, M2 with upper triangle having real entries. Clean method.*/
static GEN
rM_upper_add(GEN M1, GEN M2)
{
  long n = lg(M1), i, j;
  GEN P = cgetg(n, t_MAT);
  for (i = 1; i < n; i++) {
    gel(P, i) = cgetg(n, t_COL);
    for (j = 1; j <= i; j++) gmael(P, i, j) = addrr(gcoeff(M1, j, i), gcoeff(M2, j, i));
    for (j = i + 1; j < n; j++) gmael(P, i, j) = gen_0;
  }
  return P;
}

/*Adds two square upper triangular matrices. Clean method.*/
static GEN
RgM_upper_add(GEN M1, GEN M2)
{
  long n = lg(M1), i, j;
  GEN P = cgetg(n, t_MAT);
  for (i = 1; i < n; i++) {
    gel(P, i) = cgetg(n, t_COL);
    for (j = 1; j <= i; j++) gmael(P, i, j) = gadd(gcoeff(M1, j, i), gcoeff(M2, j, i));
    for (j = i + 1; j < n; j++) gmael(P, i, j) = gen_0;
  }
  return P;
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
static GEN
deftol(long prec)
{
  return real2n((BITS_IN_LONG >> 1)*(2 - prec), prec);
}

/*Lower tolerance. Useful if we may have more precision loss and are OK with the larger range (e.g. we have a slow routine to check exact equality, and use this to sift down to a smaller set).*/
static GEN
deflowtol(long prec)
{
  return real2n(-16, prec);
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
    iii) Embed elements of Gamma into PGU(1, 1)^+ that act on the Klein model: format as GEN Xtoklein(GEN X, GEN g). The output should be [A, B], where the corresponding element in PGU(1, 1)^+ is [A, B;conj(B), conj(A)], and |A|^2-|B|^2 is positive. If |A|^2-|B|^2!=1, this is OK, and do not rescale it, as this will be more inefficient and may cause precision loss. If you have a natural embedding into PGL(2, R)^+, then you can use m2r_to_klein, which will ensure that det(M in PGL)=|A|^2-|B|^2. Also ensure that A and B are of type t_COMPLEX with t_REAL components, even if one of them is real or even integral.
    iv) Identify if an element is trivial in Gamma: format as int Xistriv(GEN X, GEN g). Since we are working in PSL, we need to be careful that -g==g, since most representations of elements in X would be for SL. Furthermore, if we do not return the true inverse in Xinv, then we have to account for this.
    v) Determine if an element of X is parabolic or not: format at int Xisparabolic(GEN X, GEN g) (typically you can compare trd(g)^2 to +/-4*nrd(g)).
    vi) Input an element representing the trivial element of X.
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
U=[elts, sides, vcors, vargs, crossind, kact, area, spair, infinite]
elts
    Elements of Gamma whose isometric circles give the sides of the normalied boundary. An infinite side corresponds to the element 0. Elements cannot be stored as type t_INT.
sides
    The ith entry is the isometric circle corresponding to elts[i], stored as an icirc. Infinite side -> 0. The first side is the closest to the origin.
vcors
    Vertex coordinates, stored as a t_COMPLEX with real components. The side sides[i] has vertices vcor[i-1] and vcor[i], appearing in this order going counterclockwise about the origin.
vargs
    The argument of the corresponding vertex, normalized to lie in [0, 2*Pi), stored as t_REAL. These are stored in counterclockwise order, BUT vargs itself is not sorted: it is increasing from 1 to crossind, and from crossind+1 to the end.
crossind
    the arc from vertex crossind to crossind+1 contains the positive x-axis. If one vertex IS on this axis, then crossind+1 gives that vertex.  Alternatively, crossind is the unique vertex such that vargs[crossind]>vargs[crossind+1], taken cyclically.
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
  GEN order = gen_indexsort(C, &tol, &cmp_icircangle);/*Order the sides by initial angle.*/
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
    long toins = order[1 + (istart + absind - 1)%lenc];/*The index of the side to insert.*/
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
    int newptdat = angle_onarc(gel(vargs, found - 1), gel(vargs, found), iptarg, tol);/*The new point wrt the previous side.*/
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
        int ignore = angle_onarc(gel(vargs, found), gel(lastcirc, 7), iptarg, tol);/*We are in phase 2, and we either miss the side to the left (so ignore and move on), or to the right (supercede previous edge).*/
        if (!ignore) {
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
          phase2 = 1;
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

/*Given an isometric circle side, this returns -1 if z is on the same side as the origin, 0 if it is on it, and 1 else. Also used in fdom_intersect, since this works for the inputs there too.*/
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
  long lG = lg(G), lU = lg(normbound_get_elts(U)), i;
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
  long lenU = lU - 1;/*Now we figure out the order of U, before merging them.*/
  long Ustart = normbound_get_cross(U)%lenU + 1;/*The side which crosses the origin. We do a linear probe to find the smallest initial angle.*/
  long Unext = Ustart%lenU + 1;
  GEN Usides = normbound_get_sides(U), tol = gdat_get_tol(gdat);
  for (;;) {
    if (gequal0(gel(Usides, Ustart))) { Ustart = Unext; break; }/*The next side is not infinite, and must be the first one.*/
    if (gequal0(gel(Usides, Unext))) { Ustart = Unext%lenU + 1; break; }/*Must move past the infinite side.*/
    if (tolcmp(gmael(Usides, Ustart, 6), gmael(Usides, Unext, 6), tol) >= 0) { Ustart = Unext; break; }/*We have crossed over.*/
    Ustart = Unext;
    Unext = Unext%lenU + 1;
  }
  GEN Uelts = normbound_get_elts(U), Ukact = normbound_get_kact(U);
  GEN order = gen_indexsort(Cnew, &tol, &cmp_icircangle);/*Order the new sides by initial angle.*/
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
    int newptdat = angle_onarc(gel(vargs, found - 1), gel(vargs, found), iptarg, tol);/*The new point wrt the previous side.*/
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
        int ignore = angle_onarc(gel(vargs, found), gel(lastcirc, 7), iptarg, tol);/*We are in phase 2, and we either miss the side to the left (so ignore and move on), or to the right (supercede previous edge).*/
        if (!ignore) {
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
          phase2 = 1;
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

/*Used for sorting C by initial angles, then terminal angles, by radius (larger r comes first).*/
static int
cmp_icircangle(void *tol, GEN c1, GEN c2)
{
  int comp1 = tolcmp(gmael(c1, 3, 6), gmael(c2, 3, 6), *((GEN*)tol));
  if (comp1) return comp1;
  return tolcmp(gmael(c2, 3, 3), gmael(c1, 3, 3), *((GEN*)tol));/*If initial points are =, we want the larger r to come first, as this will envelop the other side.*/
}

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
  long na = lg(args) - 1, l1, l2;/*l1 and l2 track the min and max indices possible.*/
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

/*Returns the edge pairing, as VECSMALL v where v[i]=j means i is paired with j (infinite sides are by convention paired with themselves). If not all sides can be paired, instead returns [v1, v2, ...] where vi=[gind, vind, side, location] is a VECSMALL. gind is the index of the unpaired side, vind is the corresponding index of an unpaired vertex (so gind==vind or gind==vind+1 mod # of sides). The point gv projects to side side, and location=-1 means it is inside U, 0 on the boundary, and 1 is outside U. If v is infinite, then location<=0 necessarily.*/
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
    if (gequal0(act)) { pair[i] = i; continue; }/*oo side, go next (we say it is paired with itself)*/
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
  GEN unpair = gen_0, unpairtol = tol;/*Adjusting tolerance for edge pairing, in case we lose too much tolerance.*/
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
    GEN newunpair = edgepairing(U, unpairtol);
    if (typ(newunpair) == t_VECSMALL) {/*All paired up!*/
      gel(U, 8) = newunpair;
      return gerepilecopy(av, U);
    }
    /*Typically our tolerance is good enough, but in rare cases (e.g.: F=nfinit(y^3 - y^2 - 4*y + 1); A=alginit(F, [16*y^2 - 48*y - 5, 264*y^2 - 880*y - 821]); elt1=[292266712979,11224175772635,19761120255325,2858856075355,9063381984802,3524132002710,20505994302672,18181504506668,13053164032091,16489355539407,5374462579812,-22448351545240]~; elt2=[-477262530497,-18328732946074,-32269300054888,-4668423819984,-14800223317337,-5754798895964,-33485656406450,-29689836243799,-21315414378062,-26926609156642,-8776331674561,36657465892164]~) the isometic circles are so close to being equal that we lose just enough tolerance.*/
    if (gequal(unpair, newunpair)) {/*In this case, we likely have hit an infinite loop with a true side pairing that was missed due to tolerance. We decrease the tolerance.*/
      unpairtol = shiftr(unpairtol, 1);
      newunpair = edgepairing(U, unpairtol);
      if (typ(newunpair) == t_VECSMALL) {/*All paired up!*/
        gel(U, 8) = newunpair;
        return gerepilecopy(av, U);
      }
    }
    unpair = newunpair;
    GEN vcors = normbound_get_vcors(U);/*Vertex coordinates.*/
    GEN elts = normbound_get_elts(U);
    long lunp = lg(unpair), lenU = lg(elts) - 1;
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
          if (Xistriv(X, Xmul(X, gel(elts, gind), gel(elts, gvind)))) {/*gv is on I(g^-1) which is part of U.*/
            /*Let the element of the side intersecting gind at v gdind be w. Then I(wg^-1) contains gv, and will give gv as an intersection vertex, so we add in wg^-1.*/
            long wind;
            if (gind == dat[2]) wind = gind%lenU + 1;
            else wind = dat[2];
            vectrunc_append(Gadd, Xmul(X, gel(elts, wind), gel(elts, gvind)));/*gel(elts, gvind) = g^-1, or at least they have the same isometric circle*/
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
      maxind <<= 1;/*Double it.*/
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


/*2: CYCLES AND SIGNATURE*/

/*Returns the set of minimal cycles of the side pairing pair. A cycle is a vecsmall [i1, i2,..., in] so that the cycle is v_i1, v_i2, ..., v_in. A cycle [-i] means that the "vertex" on side i is is a one element cycle (happens when a side is fixed).*/
static GEN
minimalcycles(GEN pair)
{
  pari_sp av = avma;
  long np1 = lg(pair), n = np1 - 1, vleft = np1, i;/*Number of sides/vertices (not counting vertices that occur on the middle of a side).*/
  GEN vind = const_vecsmall(n, 1);/*Tracking if the vertices have run out or not*/
  GEN cycles = vectrunc_init(2*np1), cyc;/*Max number of cycles, since each side could have a middle vertex. In reality the number is probably much smaller, but this is safe.*/
  long startind = 1, ind;
  for (i = 1; i < np1; i++){/*We sort the fixed sides first, as later on we would miss the ones that get removed before checking.*/
    if (pair[i] == i) vectrunc_append(cycles, mkvecsmall(-i));/*Middle of the side is fixed.*/
  }
  do {
    cyc = vecsmalltrunc_init(vleft);
    vecsmalltrunc_append(cyc, startind);/*Starting the cycle.*/
    vind[startind] = 0;
    vleft--;
    ind = smodss(pair[startind] - 2, n) + 1;/*Hit it with the side pairing and subtract 1 to reach the paired vertex.*/
    while (ind != startind) {/*Move along the cycle.*/
      vind[ind] = 0;
      vleft--;/*One less vertex*/
      vecsmalltrunc_append(cyc, ind);/*Append it*/
      ind = smodss(pair[ind] - 2, n) + 1;/*Update*/
    }
    vectrunc_append(cycles, cyc);/*New cycle.*/
    while (startind < np1) {/*Finding the next vertex we haven't eliminated.*/
      startind++;
      if (vind[startind] == 1) break;
    }
  }
  while(startind < np1);
  return gerepilecopy(av, cycles);
}

/*Returns [cycles, types], where cycles[i] has type types[i]. Type 0=parabolic, 1=accidental, m>=2=elliptic of order m. It is returned with the types sorted, i.e. parabolic cycles first, then accidental, then elliptic.*/
static GEN
minimalcycles_bytype(GEN X, GEN U, GEN Xid, GEN (*Xmul)(GEN, GEN, GEN), int (*Xisparabolic)(GEN, GEN), int (*Xistriv)(GEN, GEN))
{
  pari_sp av = avma;
  GEN G = normbound_get_elts(U);
  GEN cycles = minimalcycles(normbound_get_spair(U));/*Min cycles*/
  long lgcyc = lg(cycles), i, j;
  GEN types = cgetg(lgcyc, t_VECSMALL);/*The types.*/
  for (i = 1; i < lgcyc; i++) {
    GEN cyc = gel(cycles, i), g;
    if (cyc[1] < 0) {/*Minimal cycle that was the middle of a side. We update it to be positive now.*/
      cyc[1] = -cyc[1];
      g = gel(G, cyc[1]);
    }
    else {
      g = Xid;
      for (j = 1; j < lg(cyc); j++) g = Xmul(X, gel(G, cyc[j]), g);/*Multiply on the left by G[cyc[j]]*/
    }
    if (Xistriv(X, g)) { types[i] = 1; continue; }/*Accidental cycle, continue on.*/
    if (Xisparabolic(X, g)) { types[i] = 0; continue; }/*Parabolic cycle.*/
    long ord=1;
    GEN gpower=g;
    do {/*Finding the order of g*/
      ord++;
      gpower = Xmul(X, g, gpower);
    }
    while (!Xistriv(X, gpower));
    types[i] = ord;
  }
  GEN ordering = vecsmall_indexsort(types);
  return gerepilecopy(av, mkvec2(vecpermute(cycles, ordering), vecsmallpermute(types, ordering)));/*The return, [cycles, types]*/
}

/*Computes the signature of the fundamental domain U. The return in [g, V, s], where g is the genus, V=[m1,m2,...,mt] (vecsmall) are the orders of the elliptic cycles (all >=2), and s is the number of parabolic cycles. The signature is normally written as (g;m1,m2,...,mt;s).*/
static GEN
signature(GEN X, GEN U, GEN Xid, GEN (*Xmul)(GEN, GEN, GEN), int (*Xisparabolic)(GEN, GEN), int (*Xistriv)(GEN, GEN))
{
  pari_sp av = avma;
  GEN mcyc = minimalcycles_bytype(X, U, Xid, Xmul, Xisparabolic, Xistriv);/*The minimal cycles and their types.*/
  long nfixed = 0, lgcyc = lg(gel(mcyc, 1)), i;/*The number of fixed sides, and number of cycles+1*/
  for (i = 1; i < lgcyc; i++) {
    if (lg(gmael(mcyc, 1, i)) == 2 && gel(mcyc, 2)[i] == 2) nfixed++;
  }/*Fixed sides MUST correspond to length 1 cycles with order 2.*/
  GEN elts = normbound_get_elts(U);
  long genus = (((lg(elts) - 1 + nfixed) >> 1) - lgcyc + 2) >> 1;
  /*2*g+t+s+1_{t+s>0}=min number of generators. Initially, we have (n+k)/2, where k is the number of sides fixed by the element (nfixed) and n is the number of sides of the fdom (b/c we take one of g and g^(-1) every time). Then each accidental cycle AFTER THE FIRST removes exactly one generator. This gives the formula (if there are no accidental cycles then we are off by 1/2, but okay as >> 1 rounds down in C. We are actually overcounting by 1 if t+s=0 and there are no accidental cycles, but this is impossible as there is >=1 cycle.*/
  long s = 0, firstell = lgcyc;/*s counts the number of parabolic cycles, and firstell is the first index of an elliptic cycle.*/
  int foundlastpar = 0;
  for (i = 1; i < lgcyc; i++) {
    if (!foundlastpar) {
      if (gel(mcyc, 2)[i] == 0) continue;
      s = i - 1;/*Found the last parabolic.*/
      foundlastpar = 1;
    }
    if (gel(mcyc, 2)[i] == 1) continue;
    firstell = i;/*Found the last accidental.*/
    break;
  }
  long lgell = lgcyc - firstell + 1;
  GEN rvec=cgetg(4, t_VEC);/*[g, V, s]*/
  gel(rvec, 1) = stoi(genus);
  gel(rvec, 2) = cgetg(lgell, t_VECSMALL);
  for (i = 1; i < lgell; i++) {
    gel(rvec, 2)[i] = gel(mcyc, 2)[firstell];
    firstell++;
  }
  gel(rvec, 3) = stoi(s);
  return gerepileupto(av, rvec);
}


/*2: GEODESICS*/

/*Assumes 1<=i1,i2<=n. Finds the index halfway between i1 and i2, where we loop around at n.*/
static long
fdom_intersect_sidesmidpt(long i1, long i2, long n)
{
  if (i2 >= i1) return (i1 + i2) >> 1;
  long temp = (i1 + i2 + n) >> 1;
  if (temp > n) return temp - n;
  return temp;
}

/*Assumes l1=[a, b, 1], returns a*z[1]+b*z[2]-1, which gives the relative distance to the line l1 from the point z*/
INLINE GEN
fdom_intersect_dist(GEN l1, GEN z)
{
  return subrs(addrr(mulrr(gel(z, 1), gel(l1, 1)), mulrr(gel(z, 2), gel(l1, 2))), 1);
}

/*Given a fundamental domain U and a geodesic geod that properly intersects it, we find the two intersections of geod with U, returning [s1, s2, v1, v2, [a, b, c]], where geod[1] to geod[2] hits v1 on side s1, then v2 on side s2, and the equation of the geodesic is ax+by=c=0 or 1. If s1 is passed as non-zero, we assume that this is the first side hit. Note that in the return, s1 and s2 are GENS of type t_INT, not longs.*/
static GEN
fdom_intersect(GEN U, GEN geod, GEN tol, long s1)
{
  pari_sp av = avma;
  GEN sides = normbound_get_sides(U);
  GEN geodeq = gel(geod, 3);
  if (!signe(gel(geodeq, 3))) {/*Line through origin.*/
    if (!s1) s1 = normbound_outside(U, gel(geod, 1), tol);
    long s2 = normbound_outside(U, gel(geod, 2), tol);/*These are automatically correct.*/
    GEN v1 = line_int(geodeq, gel(sides, s1), tol);
    GEN v2 = line_int(geodeq, gel(sides, s2), tol);
    return gerepilecopy(av, mkvec5(stoi(s1), stoi(s2), v1, v2, geodeq));
  }
  GEN verts = normbound_get_vcors(U);
  int forwards, prec = lg(tol);/*Tracks whether the arc goes from v1 to v2 or v2 to v1 counterclockwise around the circle.*/
  GEN arg1 = argmod_complex(gel(geod, 1), tol), arg2 = argmod_complex(gel(geod, 2), tol);/*Arguments of the vertices of geod*/
  if (cmprr(arg1, arg2) < 0) {/*arg2 > arg1*/
    GEN diff = subrr(arg2, arg1);
    if (cmprr(diff, mppi(prec)) < 0) forwards = 1;/*arg1 -> arg2 length <pi, so correct*/
    else forwards = 0;
  }
  else {/*arg1 > arg2*/
    GEN diff = subrr(arg1, arg2);
    if (cmprr(diff, mppi(prec)) < 0) forwards = 0;/*arg2 -> arg1 length <pi, so backwards*/
    else forwards = 1;
  }
  long n = lg(sides) - 1;
  if (!s1) s1 = normbound_outside(U, gel(geod, 1), tol);
  long s2 = normbound_outside(U, gel(geod, 2), tol);
  long i1, i2, i3 = 0, i4;
  if (forwards) { i1 = s1; i2 = s2; }
  else {i1 = s2; i2 = s1; }/*Going from side i1 to i2 counterclockwise now.*/
  i1 = smodss(i1 - 2, n) + 1;/*Now, the line from vertex i1 to i2 "cuts off" our geodesic.*/
  for (;;) {/*We loop with i1 -> i2 cutting off our geodesic. We eventually find a vertex on the other side, when we stop.*/
    long i = fdom_intersect_sidesmidpt(i1, i2, n);
    int whichside = normbound_whichside(geodeq, gel(verts, i), tol);
    if (!whichside) pari_err(e_MISC,"TO DO: GEODESIC INTERSECTING VERTEX.");
    if (whichside == 1) { i3 = i; i4 = i; break; }/*Crossed over*/
    GEN dist1 = fdom_intersect_dist(geodeq, gel(verts, i));
    long ip = (i%n) + 1;
    GEN dist2 = fdom_intersect_dist(geodeq, gel(verts, ip));
    int ips = tolsigne(dist2, tol);
    if (ips == 1) { i3 = ip; i4 = ip; break; }/*i+1 crosses over.*/
    if (!ips) pari_err(e_MISC, "TO DO: GEODESIC INTERSECTING VERTEX.");
    if (cmprr(dist1, dist2) > 0) {/*vertex i is closer than i+1, so we replace i2 with i*/
      i2 = i;
    }
    else i1 = ip;
  }
  while (smodss(i3 - i1, n) > 1) {/*Sort out the first side*/
    long i = fdom_intersect_sidesmidpt(i1, i3, n);
    int whichside = normbound_whichside(geodeq, gel(verts, i), tol);
    if (!whichside) pari_err(e_MISC,"TO DO: GEODESIC INTERSECTING VERTEX.");
    if (whichside == 1) i3 = i;
    else i1 = i;
  }
  while (smodss(i2 - i4, n) > 1) {/*Sort out the second side*/
    long i = fdom_intersect_sidesmidpt(i4, i2, n);
    int whichside = normbound_whichside(geodeq, gel(verts, i), tol);
    if (!whichside) pari_err(e_MISC,"TO DO: GEODESIC INTERSECTING VERTEX.");
    if (whichside == 1) i4 = i;
    else i2 = i;
  }
  if (forwards) { s1 = i3; s2 = i2; }
  else { s1 = i2; s2 = i3; }/*Now we have our sides!*/
  GEN v1 = line_int11(geodeq, gel(sides, s1), tol);
  GEN v2 = line_int11(geodeq, gel(sides, s2), tol);
  return gerepilecopy(av, mkvec5(stoi(s1), stoi(s2), v1, v2, geodeq));
}

/*For a hyperbolic element g, this returns [p1, p2, [a, b, c]], where v1 and v2 are the unit circle such that the root geodesic of g goes from p1 to p2, and the equation of this is ax+by=c=0 or 1. If the action in the Klein model is [A, B], then the formula is the first (attracting, i.e. p2) root is (Imag(A)*I+sqrt(Real(A)^2-1))/conj(B) IF trd(g)>0. If trd(g)<0, then change the +sqrt(...) to -sqrt(...). The other root is the same with the plus/minus swapped.*/
static GEN
geodesic_klein(GEN X, GEN g, GEN (*Xtoklein)(GEN, GEN), GEN tol)
{
  pari_sp av = avma;
  GEN M = Xtoklein(X, g);/*Move to the Klein model.*/
  GEN rA = gmael(M, 1, 1), iA = gmael(M, 1, 2);/*real/imaginary parts of A*/
  GEN toroot = subrs(sqrr(rA), 1);
  if (tolsigne(toroot, tol) <= 0 ) pari_err_TYPE("The element g should be hyperbolic.", g);
  GEN rt = sqrtr(toroot);/*sqrt(rA^2-1)*/
  GEN root1 = mkcomplex(rt, iA);/*Choosing + the root*/
  GEN root2 = mkcomplex(negr(rt), iA);/*Choosing - the root*/
  GEN Bconj = mkcomplex(gmael(M, 2, 1), negr(gmael(M, 2, 2)));/*Conjugate of B*/
  root1 = divcrcr(root1, Bconj);
  root2 = divcrcr(root2, Bconj);
  GEN eq = line_from_crrcrr(root1, root2, tol);/*Line equation*/
  if (signe(rA) > 0) return gerepilecopy(av, mkvec3(root2, root1, eq));/*Positive trace*/
  return gerepilecopy(av, mkvec3(root1, root2, eq));/*Negative trace*/
  
  
}

/*Finds the image of the root geodesic of g in the fundamental domain. The return is a vector of [g, s1, s2, v1, v2, [a, b, c]], where each component runs from vertex v1 on side s1 to vertex v2 on side s2, which has equation ax+by=c=0 or 1. The components are listed in order.*/
static GEN
geodesic_fdom(GEN X, GEN U, GEN g, GEN Xid, GEN (*Xtoklein)(GEN, GEN), GEN (*Xmul)(GEN, GEN, GEN), GEN (*Xinv)(GEN, GEN), GEN gdat)
{
  pari_sp av = avma;
  GEN tol = gdat_get_tol(gdat);
  GEN spair = normbound_get_spair(U), elts = normbound_get_elts(U);
  GEN gstart = geodesic_klein(X, g, Xtoklein, tol);/*Find the initial geodesic.*/
  GEN midpt = addcrcr(gel(gstart, 1), gel(gstart, 2));
  gel(midpt, 1) = shiftr(gel(midpt, 1), -1);
  gel(midpt, 2) = shiftr(gel(midpt, 2), -1);
  GEN red = red_elt(X, U, Xid, midpt, Xtoklein, Xmul, 0, gdat);/*Reduce this to the interior.*/
  g = Xmul(X, Xmul(X, red, g), Xinv(X, red));/*Conjugate g by red, so now the root geodesic of g goes through the interior.*/
  GEN geod = geodesic_klein(X, g, Xtoklein, tol);/*Update the initial geodesic.*/
  GEN sidestart = fdom_intersect(U, geod, tol, 0);/*First intersection*/
  long maxsides = 20, nsides = 1;
  GEN sides = cgetg(maxsides + 1, t_VEC);
  gel(sides, 1) = vec_prepend(sidestart, g);/*Initialize the container for the sides.*/
  GEN side = sidestart;
  for (;;) {
    long end = itos(gel(side, 2));/*The side we ended on, */
    long s1 = spair[end];/*This must be the next side.*/
    GEN toconj = gel(elts, end);
    g = Xmul(X, Xmul(X, toconj, g), Xinv(X, toconj));/*Update g*/
    geod = geodesic_klein(X, g, Xtoklein, tol);/*New geodesic*/
    side = fdom_intersect(U, geod, tol, s1);/*Intersections*/
    if (toleq(gel(sidestart, 3), gel(side, 3), tol) && toleq(gel(sidestart, 4), gel(side, 4), tol)) break;/*Back to the start, so done!*/
    nsides++;/*Now we are an actually new side.*/
    if (nsides > maxsides) {
      maxsides <<= 1;/*Double the sides.*/
      sides = vec_lengthen(sides, maxsides);
    }
    gel(sides, nsides) = vec_prepend(side, g);
  }
  setlg(sides, nsides + 1);
  return gerepilecopy(av, sides);
}


/*2: PRESENTATION*/

/*
A word is a vecsmall, which is taken in reference to a list of elements G. A negative entry denotes taking the inverse of the corresponding index, so the word [1, -2, 5] is G[1]*G[2]^-1*G[5]. The fundamental domain gives us a set of generators, and we will manipulate this into a minimal set of generators with relations (stored as words).
*/

/*Returns the group presentation of the fundamental domain U. The return is a vector V, where the V[1] is the list of generators (a subset of U[1]). V[2] is the vector of relations. Finally the V[3] is a vector whose ith entry is the representation of the word U[1][i] in terms of V[1]. Each term in V[2] and V[3] is formatted as [i1, i2, ..., ir], which corresponds to g_{|i1|}^(sign(i1))*...*g_{|ir|}^sign(ir). Thus, one of the generators is given the entry [ind], and its inverse is given [-ind].*/
static GEN
presentation(GEN X, GEN U, GEN Xid, GEN (*Xmul)(GEN, GEN, GEN), int (*Xisparabolic)(GEN, GEN), int (*Xistriv)(GEN, GEN))
{
  pari_sp av = avma;
  /*Initially, we start with using the original indices (of U[1]), and transfer over at the very end.*/
  long lgelts, i, j;
  GEN elts = normbound_get_elts(U);
  GEN words = cgetg_copy(elts, &lgelts);/*Stores elts[i] as a word in terms of the minimal generators.*/
  for (i = 1; i < lgelts; i++) gel(words, i) = mkvecsmall(i);/*Initially, each entry is just itself.*/
  GEN mcyc = minimalcycles_bytype(X, U, Xid, Xmul, Xisparabolic, Xistriv);/*Minimal cycles by type.*/
  GEN cyc = gel(mcyc, 1), cyctype = gel(mcyc, 2), pairing = normbound_get_spair(U);/*Cycles, types, pairing*/
  long lgcyc = lg(cyc), ngens = 0;
  for (i = 1; i < lgcyc; i++) {/*The relations are in reverse.*/
    if (cyctype[i] == 0) gel(cyc, i) = mkvecsmall(0);/*Parabolic elements give no relations (for Shimura curves, they should only exist with M_2(Q))*/
    GEN cy1 = vecsmall_reverse(gel(cyc, i)), cy = cy1;/*The element, whose cyctype[i]th power is 1.*/
    for (j = 2; j <= cyctype[i]; j++) cy = vecsmall_concat(cy, cy1);/*Concat to make the power right.*/
    gel(cyc, i) = cy;/*The actual relations are now stored in cyc.*/
  }
  GEN H=cgetg(lgelts, t_VECSMALL);/*We start by taking only one of g or g^(-1) for all g. In general, H tracks which elements are left in the generating set. Of the pair that appears, we just keep the first one we find.*/
  for (i = 1; i < lgelts; i++){/*Initially, H[i]=1 if g=g^(-1), and for exactly one of [g, g^(-1)] for all other g.*/
    long pairedind = pairing[i];/*What side i is paired to.*/
    if (pairedind >= i) { H[i]=1; ngens++; continue; }/*No need to update words, we are keeping this generator for now.*/
    H[i] = 0;/*H[i]=0 and update the word to g^-1, g is the paired index.*/
    gel(words, i)[1] = -pairedind;
    presentation_update(cyc, i, gel(words, i));/*Update the cycles.*/
  }
  /*At this point, we have replaced all inverses.*/
  long accind = 1;/*accind tracks the first accidental cycle.*/
  while (accind < lgcyc && cyctype[accind] == 0) accind++;/*Find the first accidental cycle.*/
  long ellind = accind;
  while (ellind < lgcyc && cyctype[ellind] == 1) ellind++;/*Find the first elliptic cycle.*/
  long naccident = ellind-accind;/*How many accidental cycles*/
  GEN r = gen_0;/*The accidental relation, if it exists.*/
  if (naccident) r = gel(cyc, accind);/*The first accidental relation.*/
  if (naccident > 1) {/*More than one relation. We go through the relations, updating r by solving for elements it has in common in another relation.*/
    ngens = ngens - naccident + 1;/*Updating the number of generators.*/
    long lastrel = ellind - 1;/*The last relation we consider. We swap this with the relation we find every step and decrease it, so we only need to consider relations from accind+1 to lastrel at each step.*/
    long indrep=1, torep;/*indrep stores the index replaced in r. On the next pass, we may as well start from there, as we have already checked the previous indices for replacement. No collapsing occurs, since each index and its inverse appear at most once in the relations (and cannot disappear, unless we cancel them).*/
    GEN repind = gen_0, cycle;/*repind stores [j, l, m], where term j in relation r is replaced by using term m of relation l.*/
    /*Now we look for a common term.*/
    for (i = 1; i < naccident; i++) {/*Each step we solve the relation.*/
      if (gc_needed(av, 2)) {
        gerepileall(av, 5, &cyc, &cyctype, &H, &words, &r);
      }
      long lr = lg(r);
      for (j = indrep; j < lr; j++) {/*Trying to replace index j. We stop once we find one replacable.*/
        torep = r[j];/*What we try to replace.*/
        long mtorep = -torep, l, m;/*The inverse.*/
        for (l = accind + 1; l <= lastrel; l++) {/*Looking at cycle l*/
          for(m = 1; m < lg(gel(cyc, l)); m++) {/*element m of cycle l*/
            long thisind = gel(cyc, l)[m];
            if(thisind != torep && thisind != mtorep) continue;
            if (torep > 0) H[torep] = 0;/*Make sure it has been deleted from H.*/
            else H[mtorep] = 0;
            repind = mkvecsmall3(j, l, m);
            l = lastrel+1;/*Break loop*/
            j = lr;/*Break loop*/
            break;/*Break loop*/
          }
        }
      }
      /*Now, repind gives us all the information we need to do the replacement.*/
      cycle = gel(cyc, repind[2]);/*The cycle being deleted.*/
      long lcyc = lg(cycle);
      GEN repl = cgetg(lcyc-1, t_VECSMALL);/*cycle[repind[3]]. we have a1*...*an=1 -> a_m=a_{m-1}^(-1)*...*a_1^(-1)*a_n^(-1)*a_{n-1}^(-1)*...*a_{m+1}^(-1).*/
      long irepl = 1;
      for (j = repind[3] - 1; j > 0; j--) { repl[irepl] = -cycle[j]; irepl++; }
      for (j = lcyc - 1; j > repind[3]; j--) { repl[irepl] = -cycle[j]; irepl++; }/*Replace index cycle[repind[3]] with the word repl.*/
      presentation_update(words, cycle[repind[3]], repl);/*Replace the words.*/
      r = word_substitute(r, cycle[repind[3]], repl, NULL);/*Replace it in r. No need to replace it in the other relations, since it cannot appear.*/
      indrep = repind[1];/*The index we replaced.*/
      gel(cyc, repind[2]) = gel(cyc, lastrel);/*Replace the relation we replaced with the last one.*/
      lastrel--;/*One less relation.*/
    }
  }
  /*Now we are almost done. If there is an elliptic cycle, we can do one final replacement.*/
  if (ellind < lgcyc) {
    ngens--;/*One less generator.*/
    long lr = lg(r), irel = 0, rind = 0;/*r[rind] shares an index with cyc[irel].*/
    for (i = ellind; i < lgcyc; i++) {/*Testing elliptic relation i. In the generic case, all elliptic relations have length 1 (and all accidental have length 3). However, this does NOT need to be the case if our point p happens to be from a certain set (c.f. Voight Prop 5.4). For an explicit example, F=nfinit(y);A=alginit(F, [33, -1]);X=afuchinit(A);p=I/2*/
      long maxj = (lg(gel(cyc, i)) - 1)/cyctype[i];/*Since the elliptic relation is repeated cyctype[i] times, we only need to check the first grouping of indices.*/
      for (j = 1; j <= maxj; j++) {
        long ind = gel(cyc, i)[j], mind = -ind, rinc;
        for (rinc = 1; rinc < lr; rinc++) {/*There MUST be an accidental cycle, so r is non-zero.*/
          if (r[rinc] != ind && r[rinc] != mind) continue;/*Not equal*/
          irel = i;
          rind = rinc;
          j = maxj + 1;/*break*/
          i = lgcyc;/*break*/
          break;/*break*/
        }
      }
    }
    /*Now we solve for the index with r.*/
    if (r[rind] > 0) H[r[rind]] = 0;/*Deleting from H.*/
    else H[-r[rind]] = 0;
    GEN repl = cgetg(lr - 1, t_VECSMALL);
    long irepl = 1;
    for (j = rind - 1; j > 0; j--) { repl[irepl] = -r[j]; irepl++; }
    for (j = lr - 1; j > rind; j--) { repl[irepl] = -r[j]; irepl++; }/*The replacment for r.*/
    presentation_update(words, r[rind], repl);/*Replace the words.*/
    gel(cyc, irel) = word_substitute(gel(cyc, irel), r[rind], repl, NULL);/*Update the elliptic relation.*/
  }
  /*Ok, time to finish. H tracks which indices remain, words track the words, and we have some relations. We have to change the indices to in terms of the remaining generators.*/
  GEN gens = cgetg(ngens + 1, t_VEC);/*The generators.*/
  long igens = 1, Hind = 1;
  for (i = 1; i < lgelts; i++) {
    if (H[i] == 1) {
      H[i] = Hind;
      Hind++;
      gel(gens, igens) = gel(elts, i);
      igens++;
    }
  }/*Index i is replaced with index H[i].*/
  /*Let's do words first.*/
  for (i = 1; i < lgelts; i++) {
    for (j = 1; j < lg(gel(words, i)); j++) {
      long k = gel(words, i)[j];
      if (k > 0) gel(words, i)[j] = H[k];
      else gel(words, i)[j] = -H[-k];
    }
  }
  /*Now, relations*/
  GEN relations;
  if (ellind < lgcyc) {/*There were elliptic relations, and this is all that remains (since the final accidental relation was combined into an elliptic one).*/
    long nell = lgcyc - ellind + 1;/*Number of relations;*/
    relations = cgetg(nell, t_VEC);
    long relind = 1;
    for (i = ellind; i < lgcyc; i++) {
      long lr;
      gel(relations, relind) = cgetg_copy(gel(cyc, i), &lr);
      for (j = 1; j < lr; j++) {
        long k = gel(cyc, i)[j];
        if (k > 0) gel(relations, relind)[j] = H[k];
        else gel(relations, relind)[j] = -H[-k];
      }
      relind++;
    }
  }
  else {/*No elliptic relations, so just r leftover.*/
    long lr = lg(r);
    for (i = 1; i < lr; i++) {
      long k = r[i];
      if(k > 0) r[i] = H[k];
      else r[i] = -H[-k];
    }
    relations = mkvec(r);
  }
  return gerepilecopy(av, mkvec3(gens, relations, words));
}

/*Perform word_substitute(word, ind, repl) on each word in words, and updates it. This is OK since word is a GEN vector, so we are updating the memory addresses.*/
static void
presentation_update(GEN words, long ind, GEN repl)
{
  GEN invrepl = word_inv(repl);
  long lwords = lg(words), i;
  for (i = 1; i < lwords; i++) gel(words, i)=word_substitute(gel(words, i), ind, repl, invrepl);/*Substitute!*/
}

/*Given a word, "collapses" it down, i.e. deletes consecutive indices of the form i, -i, and repeats until the word is reduced.*/
static GEN
word_collapse(GEN word)
{
  pari_sp av = avma;
  long lgword = lg(word), n = lgword - 1, last1 = 0, newlg = lgword, i;/*last1 tracks the last value of 1 in H that we have tracked. newlg tracks the new length.*/
  GEN H = const_vecsmall(n, 1);/*H tracks if we keep each index in the collapsed word.*/
  for (i = 1; i < n; i++) {/*We go through the word.*/
    if (word[i] != -word[i+1]) {last1 = i;continue;}/*Nothing to do here, move on. H[i]=1 is necessarily true.*/
    H[i] = 0;
    H[i+1] = 0;/*Deleting i and i+1*/
    newlg -= 2;/*2 less entries*/
    if (last1 == 0) {i++;continue;}/*All previous entries are already deleted, along with i+1, so just move on to i+2.*/
    long next1 = i + 2;/*Next entry with 1*/
    while (last1 > 0 && next1 < lgword) {
      if (word[last1] != -word[next1]) break;/*Maximum collapsing achieved.*/
      H[last1] = 0;
      H[next1] = 0;
      newlg -= 2;/*2 less entries*/
      while (last1 > 0 && H[last1] == 0) last1--;/*Going backwards*/
      next1++;/*Going forwards. We are guaranteed that next1=lgword or H[next1]=1.*/
    }
    i = next1 - 1;/*We do i++ at the end of the loop, so we want to go on to i=next1 next.*/
  }
  GEN newword = cgetg(newlg, t_VECSMALL);
  long nind = 1;
  for (i = 1; i < lgword; i++) {
    if (H[i] == 0) continue;/*Go on*/
    newword[nind] = word[i];
    nind++;
  }
  return gerepileupto(av, newword);
}

/*Inverse of the word.*/
static GEN
word_inv(GEN word)
{
  long lword, i;
  GEN invword = cgetg_copy(word, &lword);/*The inverse of word*/
  long invind = lword;
  for (i = 1; i < lword; i++) {/*Negate and reverse word*/
    invind--;
    invword[invind] = -word[i];
  }
  return invword;
}

/*We replace all instances of the index ind in word with repl, and then collapse down to a reduced word (with regards to the free group on the indices). We return the new word. Thus, substitute_word([1,2,3,-4,5], 3, [-2, -2, 5, 4])=[1, -2, 5, 5]. invrepl can be passed as NULL.*/
static GEN
word_substitute(GEN word, long ind, GEN repl, GEN invrepl)
{
  pari_sp av = avma;
  long lrepl = lg(repl), i;/*length of the replacement word*/
  if (!invrepl) invrepl = word_inv(repl);
  long nreplace = 0, mind = -ind, lgword;/*nreplace counts how many replacement we need to perform, and mind is -ind.*/
  GEN replaceplaces = cgetg_copy(word, &lgword);/*index 1 means replace this with repl,-1 with invrepl, and 0 means keep.*/
  for (i = 1; i < lgword; i++) {
    if (word[i] == ind) { replaceplaces[i] = 1; nreplace++; }
    else if (word[i] == mind) { replaceplaces[i] = -1; nreplace++; }
    else replaceplaces[i] = 0;
  }
  if (nreplace==0) return gerepileupto(av, word_collapse(word));/*Nothing to replace! Just collapse the word down.*/
  long newwordlg = lgword + nreplace*(lrepl - 2);/*The new lg*/
  GEN newword = cgetg(newwordlg, t_VECSMALL);/*The new word*/
  long newind = 1, j;
  for (i = 1; i < lgword; i++){/*Make the new word*/
    switch (replaceplaces[i]) {
      case 0:/*Move this index to the new word.*/
        newword[newind] = word[i];
        newind++;
        break;
      case 1:
        for (j = 1; j < lrepl; j++) {
          newword[newind] = repl[j];
          newind++;
        }
        break;
      default:/*Case -1, add in invrepl.*/
        for (j = 1; j < lrepl; j++) {
          newword[newind] = invrepl[j];
          newind++;
        }
    }
  }
  return gerepileupto(av, word_collapse(newword));/*Collapse the finalword.*/
}

/*Writes g as a word in terms of the presentation.*/
static GEN
word(GEN X, GEN U, GEN P, GEN g, GEN (*Xtoklein)(GEN, GEN), GEN (*Xmul)(GEN, GEN, GEN), GEN (*Xinv)(GEN, GEN), int (*Xistriv)(GEN, GEN), GEN gdat)
{
  pari_sp av = avma;
  GEN tol = gdat_get_tol(gdat);
  long prec = lg(tol);
  GEN origin = gtocr(gen_0, prec);
  GEN ginv = Xinv(X, g);/*g^-1*/
  GEN gred = red_elt_decomp(X, U, ginv, origin, Xtoklein, Xmul, gdat);/*Reduction of g^-1, which gives a word for g.*/
  if (!Xistriv(X, gel(gred, 1))) pari_err(e_MISC, "We could not reduce the element to the identity. Increase the precision perhaps?");
  GEN oldword = gel(gred, 2);/*g as a word in U[1]. We must move it to a word in P[1].*/
  GEN Ptrans = gel(P, 3);/*U[1] in terms of P.*/
  long newlg = 1, lold = lg(oldword), i;
  for (i = 1; i < lold; i++) newlg = newlg + lg(gel(Ptrans, oldword[i])) - 1;/*Getting the new length.*/
  GEN newword = cgetg(newlg, t_VECSMALL);
  long inew = 1;
  for (i = 1; i < lold; i++) {
    GEN repl = gel(Ptrans, oldword[i]);/*What to sub in.*/
    long lrepl = lg(repl), j;
    for (j = 1; j < lrepl; j++) { newword[inew] = repl[j]; inew++; }
  }
  return gerepileupto(av, word_collapse(newword));
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
An arithmetic Fuchsian group initialization is represented by X, and is stored as a lazy vector. In general, elements of the algebra A are represented in terms of the basis of O, NOT in terms of the basis of A. If O is the identity, then these notions are identical. We have:
X
    [O, Oinv, Oconj, Omultable, Onormdat, type [A, Onormreal, kleinmats, qfmats, gdat, fdomdat, fdom, sig, pres, savedelts, normalizernorms]]

O, Oinv
    The order and its inverse, given as a matrix whose columns generate the order (with respect to the stored order in A).
Oconj
    Conjugates of the elements of O. To conjugate an element x of A, it suffices to compute Oconj*x.
Omultable
    To multiply elements of O more quickly.
Onormdat
    Decomposition used to compute norms of elements more quickly.
type
    Which symmetric space we want to compute. 0 is O^1, 1 is totally positive reduced norm, 2 is AL(O), 3 is N_{B^{\times}}^+(O).
A
    The algebra 
Onormreal
    Used to compute real approximations to norms of elements, which is much more efficient when deg(F)>1.
kleinmats
    O[,i] is sent to embmats[i] which acts on the unit disc/Klein model.
qfmats
    For supply into afuch_make_qf: they help make the quadratic form Q_{z, 0}(g), where if g has norm 1, then Q_{z, 0}(g)=cosh(d(gz, 0))+n-1 (n=deg(F), F is the centre of A).
gdat
    Geometric data.
fdomdat
    [area, C, R, epsilon, passes]: Constants used specifically for computing the fundamental domain. Area is the area of the domain, C is the optimal constant for finding elements of norm 1, R is the radius used to pick random points, epsilon is the growth rate of R (in case we stagnate), and passes is the "expected" number of times we generate new points.
fdom
    Fundamental domain, if computed.
sig
    Signature, if computed.
presentation
    Presentation, if computed.
savedelts
    [O1elts, totposelts, ALelts, normelts]. If type=3, we initialize this, which saves the data required to compute any fundamental domain between O^1 and N_{B^{\times}}^+(O). If type <3, we don't bother initializing this. If you want to compute a variety of types for the same algebra, you should initialize type=3 first.
normalizernorms
    Only initialized if you call afuchnormalizernorms or initialize with type=3. This is the collection [n1, n2, n3], where ni is a collection of norms possible in the positive normalizer of the Eichler order O. n1 is a set of generators for the unit norms, n2 is for the Atkin-Lehner norms, and n3 finishes off the entire normalizer group.
*/

/*Clean initialization of data for the fundamental domain. Can pass p=NULL and will set it to the default, O=NULL gives the stored maximal order, and type=NULL gives norm 1 group. flag>0 means we also initialize the fundamental domain, and flag=2 means we do the signature and presentation as well.*/
GEN
afuchinit(GEN A, GEN O, GEN type, int flag, long prec)
{
  pari_sp av = avma;
  GEN F = alg_get_center(A);
  long nF = nf_get_degree(F);
  if (!O) O = matid(nF << 2);
  if (!type) type = gen_0;
  GEN p = defp(prec);
  GEN gdat = gdat_initialize(p, prec);
  GEN Oinv, AX = obj_init(6, 11);
  gel(AX, 1) = O;
  gel(AX, 2) = Oinv = QM_inv(O);/*Inverse*/
  long lgO = lg(O), i, j;
  GEN AOconj = cgetg(lgO, t_MAT);/*Conjugates of basis of O, written in A.*/
  for (i = 1; i < lgO; i++) gel(AOconj, i) = algconj(A, gel(O, i));
  gel(AX, 3) = QM_mul(Oinv, AOconj);/*Matrix of conjugates; write it in terms of O.*/
  gel(AX, 4) = Omultable(A, O, Oinv);
  GEN Onorm = Onorm_makemat(A, O, AOconj);
  if (nF >= 5) {/*More efficient to use the Cholesky norm.*/
    gel(AX, 5) = Onorm_makechol(F, Onorm);
    if (gequal0(gel(AX, 5))) gel(AX, 5) = gcopy(Onorm);/*I don't think this can ever trigger, but just in case we can't get the Cholesky version, we use Onorm as backup.*/
  }
  else gel(AX, 5) = gcopy(Onorm);/*Testing showed that Cholesky was faster if deg(F)>=5, and Onorm if deg(F)<=4. Both were much faster than algnorm.*/
  gel(AX, 6) = gcopy(type);
  obj_insert(AX, afuch_A, A);
  if (nF > 1) obj_insert(AX, afuch_ONORMREAL, Onorm_toreal(A, Onorm));/*No need to store approximate if n=1.*/
  GEN kleinmats = afuch_make_kleinmats(A, O, gdat_get_p(gdat), prec);
  obj_insert(AX, afuch_KLEINMATS, kleinmats);/*Make sure p is safe.*/
  GEN qfm = afuch_make_qfmats(kleinmats);
  for (i = 1; i < lgO; i++) {/*Making the traces in Onorm.*/
    for (j = i; j < lgO; j++) gcoeff(Onorm, i, j) = nftrace(F, gcoeff(Onorm, i, j));
  }
  gel(qfm, 5) = Onorm;
  obj_insert(AX, afuch_QFMATS, qfm);
  obj_insert(AX, afuch_GDAT, gdat);
  obj_insert(AX, afuch_FDOMDAT, afuchfdomdat_init(A, O, prec));
  if (flag) {
    afuchfdom(AX);
    if (flag == 2) { afuchsignature(AX); afuchpresentation(AX); }
  }
  return gerepilecopy(av, AX);
}

/*Updates X to the new value of p.*/
void
afuch_changep(GEN X, GEN p)
{
  pari_sp av = avma;
  GEN A = afuch_get_alg(X);
  GEN O = afuch_get_O(X);
  long prec = afuch_get_prec(X);
  p = uhp_safe(p, prec);
  GEN kleinmats = afuch_make_kleinmats(A, O, p, prec);
  obj_insert(X, afuch_KLEINMATS, kleinmats);
  GEN qfm = afuch_make_qfmats(kleinmats);
  gel(qfm, 5) = gcopy(gel(afuch_get_qfmats(X), 5));/*This part was exact, so we can copy it over.*/
  obj_insert(X, afuch_QFMATS, qfm);
  GEN gdat = gcopy(afuch_get_gdat(X));
  gel(gdat, 2) = p;
  obj_insert(X, afuch_GDAT, gdat);
  GEN U = afuch_get_fdom(X);
  if (U) {/*Recompute the fundamental domain.*/
    U = normbasis(X, NULL, normbound_get_elts(U), &afuchtoklein, &afuchmul, &afuchconj, &afuchistriv, gdat);
    obj_insert(X, afuch_FDOM, U);
    GEN pres = afuch_get_pres(X);
    if (pres) {/*Recompute the presentation.*/
      pres = presentation(X, U, afuchid(X), &afuchmul, &afuchisparabolic, &afuchistriv);
      obj_insert(X, afuch_PRES, pres);
    }
  }
  set_avma(av);
}

/*Updates X to have precision prec+inc. Does not initialize the fundamental domain, signature, or precision: this will typically be triggered when we find that we don't have enough precision for qfminim when computing the fundamental domain.*/
static void
afuch_moreprec(GEN X, long inc)
{
  pari_sp av = avma, av2;
  if (inc < 0) inc = 1;
  GEN old_gdat = afuch_get_gdat(X);
  long old_prec = lg(gdat_get_tol(old_gdat));
  av2 = avma;
  GEN old_A = afuch_get_alg(X);
  GEN old_nf = alg_get_center(old_A);
  long nF = nf_get_degree(old_nf);
  GEN old_rnf = alg_get_splittingfield(old_A);
  GEN aut = alg_get_aut(old_A);
  GEN hi = alg_get_hasse_i(old_A);
  GEN hf = alg_get_hasse_f(old_A);
  GEN b = alg_get_b(old_A);/*The basic data*/
  long new_prec = old_prec + inc;/*The new precision*/
  GEN new_nf = nfinit(nf_get_pol(old_nf), new_prec);/*Recompile the number field*/
  GEN new_rnf = rnfinit(new_nf, rnf_get_pol(old_rnf));/*Recompile the rnf*/
  GEN bden = Q_denom(lift(b));
  GEN new_A;/*The new algebra.*/
  if (!isint1(bden)) {/*b has a denominator. We use alg_complete, but this MAY change the value of b, so we run it until it does not.*/
    do {
      new_A = alg_complete(new_rnf, aut, hf, hi, 0);/*Remake the algebra.*/
    } while (!gequal(alg_get_b(new_A), b));/*Sometimes this changes b, so we loop until it doesn't.*/
  }
  else new_A = alg_cyclic(new_rnf, aut, b, 0);/*Initialize with alg_cyclic*/
  GEN basis_want = alg_get_basis(old_A);
  GEN basis_dontwant = alg_get_invbasis(new_A);
  GEN basis_change = QM_mul(basis_dontwant, basis_want);
  new_A = my_alg_changeorder(new_A, basis_change);/*Correct the basis to the old one*/
  new_A = gerepileupto(av2, new_A);
  /*We have the new algebra*/
  obj_insert(X, afuch_A, new_A);
  GEN O = afuch_get_O(X);
  GEN AOconj = QM_mul(O, afuch_get_Oconj(X));
  GEN Onorm = Onorm_makemat(new_A, O, AOconj);
  if (nF > 1) obj_insert(X, afuch_ONORMREAL, Onorm_toreal(new_A, Onorm));/*No need to store approximate if n=1.*/
  GEN p = gtocr(gdat_get_p(old_gdat), new_prec);
  GEN gdat = gdat_initialize(p, new_prec);
  GEN kleinmats = afuch_make_kleinmats(new_A, O, p, new_prec);
  obj_insert(X, afuch_KLEINMATS, kleinmats);
  GEN qfm = afuch_make_qfmats(kleinmats);
  gel(qfm, 5) = gcopy(gel(afuch_get_qfmats(X), 5));/*This part was exact, so we can copy it over.*/
  obj_insert(X, afuch_QFMATS, qfm);
  obj_insert(X, afuch_GDAT, gdat);
  obj_insert(X, afuch_FDOMDAT, afuchfdomdat_init(new_A, O, new_prec));
  set_avma(av);
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

/*Q_{z, 0}(g)=Q_1(g)+Tr_{F/Q}(nrd(g)) is used to find elements such that gz is near 0 (see Definition 4.2 of my paper). Let c=(1/sqrt(1-|z|^2)+1), let w=z/(1+sqrt(1-|z|^2)), and let w=w_1+w_2i. This method returns a vector of symmetric matrices [Msq, M1, M2, Mconst, 0] such that Q_1(g)=c(Msq*(w1^2+w2^2)+M1*w1+M2*w2+Mconst), which enables for efficient computation of Q_1 given z. We only store the upper half, since they are symmetric. We leave the entry 0 since we will insert the other half of the quadratic form, Tr_{F/Q}(nrd(g)), seperately.*/
static GEN
afuch_make_qfmats(GEN kleinmats)
{
  pari_sp av = avma;
  long lk = lg(kleinmats), n = lk - 1, i, j;/*Let kleinmats[i]=[Ei, Fi].*/
  GEN Msq = zeromatcopy(n, n);
  GEN M1 = zeromatcopy(n, n);
  GEN M2 = zeromatcopy(n, n);
  GEN Mconst = zeromatcopy(n, n);
  for (i = 1; i < lk; i++) {
    GEN rEi = gmael3(kleinmats, i, 1, 1), imEi = gmael3(kleinmats, i, 1, 2);/*real and imaginary parts of Ei*/
    GEN rFi = gmael3(kleinmats, i, 2, 1), imFi = gmael3(kleinmats, i, 2, 2);/*real and imaginary parts of Fi*/
    for (j = i; j < lk; j++) {
      GEN rEj = gmael3(kleinmats, j, 1, 1), imEj = gmael3(kleinmats, j, 1, 2);/*real and imaginary parts of Ej*/
      GEN rFj = gmael3(kleinmats, j, 2, 1), imFj = gmael3(kleinmats, j, 2, 2);/*real and imaginary parts of Fj*/
      GEN t = addrr(mulrr(rEi, rEj), mulrr(imEi, imEj));
      gcoeff(Msq, i, j) = t;/*rEi*rEj+imEi+imEj*/
      t = addrr(addrr(mulrr(rEi, rFj), mulrr(rEj, rFi)), addrr(mulrr(imEi, imFj), mulrr(imEj, imFi)));
      gcoeff(M1, i, j) = t;/*rEi*rFj+rEj*rFi+imEi*imFj+imEj*imFi*/
      t = subrr(addrr(mulrr(rEi, imFj), mulrr(rEj, imFi)), addrr(mulrr(imEi, rFj), mulrr(imEj, rFi)));
      gcoeff(M2, i, j) = t;/*rEi*imFj+rEj*imFi-imEi*rFj-imEj*rFi*/
      t = addrr(mulrr(rFi, rFj), mulrr(imFi, imFj));
      gcoeff(Mconst, i, j) = t;/*rFi*rFj+imFi*imFj*/
    }
  }
  return gerepilecopy(av, mkvec5(Msq, M1, M2, Mconst, gen_0));
}

/*Makes the upper triangular matrix, whose symmetrization corresponds to Tr_{F/Q)(nrd(g)/nm). Can pass in Onorm, which is the matrix for nrd(g) on its own.*/
static GEN
afuch_make_traceqf(GEN X, GEN nm, GEN Onorm)
{
  pari_sp av = avma;
  long lgO, i, j;
  GEN A = afuch_get_alg(X), F = alg_get_center(A);
  if (!Onorm) {/*Finding Onorm if not passed in.*/
    GEN O = afuch_get_O(X);
    GEN AOconj = cgetg_copy(O, &lgO);/*Conjugates of basis of O, written in A.*/
    for (i = 1; i < lgO; i++) gel(AOconj, i) = algconj(A, gel(O, i));
    Onorm = Onorm_makemat(A, O, AOconj);
  }
  GEN nminv = nfinv(F, nm);
  GEN M = cgetg_copy(Onorm, &lgO);
  for (i = 1; i < lgO; i++) {/*Making the traces.*/
    gel(M, i) = cgetg(lgO, t_COL);
    for (j = 1; j <= i; j++) gcoeff(M, j, i) = nftrace(F, nfmul(F, gcoeff(Onorm, j, i), nminv));
    for (j = i + 1; j < lgO; j++) gcoeff(M, j, i) = gen_0;
  }
  return gerepileupto(av, M);
}

/*Returns a new Fuchsian group with the same algebra and order, but a new type. This entry will point to some of the data on the heap for the original one, further enforcing the requirement that you do not access and change the components directly.*/
GEN
afuch_newtype(GEN X, GEN type)
{
  pari_sp av = avma;
  GEN newX = obj_init(6, 11);
  gel(newX, 1) = afuch_get_O(X);
  gel(newX, 2) = afuch_get_Oinv(X);
  gel(newX, 3) = afuch_get_Oconj(X);
  gel(newX, 4) = afuch_get_Omultable(X);
  gel(newX, 5) = afuch_get_Onormdat(X);
  gel(newX, 6) = type;
  obj_insert(newX, afuch_A, afuch_get_alg(X));
  GEN Onormreal = afuch_get_Onormreal(X);
  if (Onormreal) obj_insert(newX, afuch_ONORMREAL, Onormreal);
  obj_insert(newX, afuch_KLEINMATS, afuch_get_kleinmats(X));
  obj_insert(newX, afuch_QFMATS, afuch_get_qfmats(X));
  obj_insert(newX, afuch_GDAT, afuch_get_gdat(X));
  obj_insert(newX, afuch_FDOMDAT, obj_check(X, afuch_FDOMDAT));
  GEN saved = afuch_get_savedelts(X);
  if (saved) obj_insert(newX, afuch_SAVEDELTS, saved);
  GEN NN = afuch_get_normalizernorms(X);
  if (NN) obj_insert(newX, afuch_NORMALIZERNORMS, NN);
  GEN U = afuch_get_fdom(X);
  if (U) {
    afuchfdom(newX);
    GEN S = afuch_get_sig(X);
    if (S) afuchsignature(newX);
    GEN P = afuch_get_pres(X);
    if (P) afuchpresentation(newX);
  }
  return gerepilecopy(av, newX);
}

/*Initializes a table to compute multiplication of elements written in terms of O's basis more efficiently.*/
static GEN
Omultable(GEN A, GEN O, GEN Oinv)
{
  long n = lg(O), i, j;
  GEN M = cgetg(n, t_VEC);
  for (i = 1; i < n; i++) {
    gel(M, i) = cgetg(n, t_MAT);
    for (j = 1; j < n; j++) {
      gmael(M, i, j) = QM_QC_mul(Oinv, algmul(A, gel(O, i), gel(O, j)));/*Move it back to O*/
    }
  }
  return M;
}

/*If possible, this initializes the Cholesky decomposition of the norm form on O, which can be evaluated with afuchnorm_chol. This is typically faster than afuchnorm_mat if deg(F)>=5, and slower if deg(F)<=4. Returns 0 if not possible to make.*/
static GEN
Onorm_makechol(GEN F, GEN Onorm)
{
  pari_sp av = avma;
  long n = lg(Onorm) - 1, i, j, k;/*Onorm is nxn*/
  GEN M = RgM_shallowcopy(Onorm);
  GEN inds = vecsmalltrunc_init(5);/*At most 4 places where we get something non-zero.*/
  GEN vecs = vectrunc_init(5);
  for (i = 1; i <= n; i++) {
    if (gequal0(gcoeff(M, i, i))) {
      for (j = i + 1; j <= n; j++) if (!gequal0(gcoeff(M, i, j))) return gc_const(av, gen_0);/*Bad! Not sure this will ever trigger though.*/
      continue;
    }
    GEN v = cgetg(n + 1, t_VEC);
    for (j = 1; j < i; j++) gel(v, j) = gen_0;
    gel(v, i) = gcoeff(M, i, i);
    for (j = i + 1; j <= n; j++) gel(v, j) = nfdiv(F, gcoeff(M, i, j), gel(v, i));/*M[i,j]=M[i,j]/M[i,i]*/
    for (j = i + 1; j <= n; j++) {
      for (k = j; k <= n; k++) gcoeff(M, j, k)=nfsub(F, gcoeff(M, j, k), nfmul(F, gcoeff(M, i, j), gel(v, k)));/*M[j,k]=M[j,k]-M[i,i]*M[i,j]*M[i,k];*/
    }
    vecsmalltrunc_append(inds, i);
    vectrunc_append(vecs, v);
  }
  return gerepilecopy(av, mkvec2(inds, vecs));
}

/*Returns the upper half part of the symmetric matrix M (with coefficients in centre(A)) such that nrd(g in basis O)=g~*M*g. Ensures that this upper half part are written in terms of the basis representation for F.*/
static GEN
Onorm_makemat(GEN A, GEN O, GEN AOconj)
{
  GEN F = alg_get_center(A);
  long lO = lg(O), i, j;
  GEN M = cgetg(lO, t_MAT);
  for (i = 1; i < lO; i++) {
    gel(M, i) = cgetg(lO, t_COL);
    for (j = 1; j < i; j++) gmael(M, i, j) = algtobasis(F, gdivgs(lift(algtrace(A, algmul(A, gel(O, i), gel(AOconj, j)), 0)), 2));
    gmael(M, i, i) = algtobasis(F, algnorm(A, gel(O, i), 0));
    for (j = i + 1; j < lO; j++) gmael(M, i, j) = gen_0;
  }
  return M;
}

/*Given Onorm, which gives the norm, evaluates the entries at the unique split place, returning the matrix with real entries. This will compute the approximate norm much more quickly when F!=Q.*/
static GEN
Onorm_toreal(GEN A, GEN Onorm)
{
  pari_sp av = avma;
  long split = algsplitoo(A), n = lg(Onorm) - 1, i, j;
  GEN F = alg_get_center(A);
  long Fvar = nf_get_varn(F);
  GEN rt = gel(nf_get_roots(F), split);
  GEN Onormreal = zeromatcopy(n, n);
  for (i = 1; i <= n; i++) {
    for (j = i; j <=n; j++) gcoeff(Onormreal, i, j) = gsubst(lift(basistoalg(F, gcoeff(Onorm, i, j))), Fvar, rt);
  }
  return gerepilecopy(av, Onormreal);
}


/*3: ALGEBRA FUNDAMENTAL DOMAIN CONSTANTS*/

/*Returns the area of the fundamental domain of Q, computed to computeprec (Olevel is the output of algorderlevel(A, O, 1)). Assumes the order is Eichler; see Voight Theorem 39.1.8 to update for the theorem. We also submit the old precision since the loss of precision causes errors later.*/
static GEN
afuchO1area(GEN A, GEN O, GEN Olevel_fact, long computeprec, long prec)
{
  pari_sp av = avma;
  long bits = bit_accuracy(computeprec);
  GEN F = alg_get_center(A), pol = nf_get_pol(F);
  GEN zetaval = lfun(pol, gen_2, bits);/*zeta_F(2)*/
  GEN rams = algramifiedplacesf(A);
  long np = lg(rams), i;
  GEN norm=gen_1;
  for (i = 1; i < np; i++) norm = mulii(norm, subis(idealnorm(F, gel(rams, i)), 1));/*Product of N(p)-1 over finite p ramifying in A*/
  GEN elevpart = gen_1;
  np = lg(gel(Olevel_fact, 1));
  for (i = 1; i < np; i++) {/*We have an Eichler part for each i that triggers.*/
    GEN Np = idealnorm(F, gcoeff(Olevel_fact, i, 1));/*Norm of the prime*/
    GEN Npexp = gcoeff(Olevel_fact, i, 2);/*Exponent*/
    if (equali1(Npexp)) {
      elevpart = mulii(elevpart, addis(Np, 1));/*Times N(p)+1*/
      continue;
    }
    GEN Npem1 = powii(Np, gsubgs(Npexp, 1));/*Np^{e-1}*/
    GEN Npe = mulii(Npem1, Np);/*Np^e*/
    elevpart = mulii(elevpart, addii(Npe, Npem1));
  }/*Product of N(p)^e*(1+1/N(p)) over p^e||level.*/
  GEN ar = gmul(gpow(nfdisc(pol), mkfracss(3, 2), computeprec), norm);/*d_F^(3/2)*phi(D)*/
  ar = gmul(ar, elevpart);/*d_F^(3/2)*phi(D)*psi(M)*/
  long n = nf_get_degree(F), twon = n << 1;
  ar = gmul(ar, zetaval);/*zeta_F(2)*d_F^(3/2)*phi(D)*/
  ar = mpshift(ar, 3 - twon);
  ar = gmul(ar, gpowgs(mppi(computeprec), 1 - twon));/*2^(3-2n)*Pi^(1-2n)*zeta_F(2)*d_F^(3/2)*phi(D)*/
  return gerepileupto(av, gtofp(ar, prec));
}

/*Generate the optimal C_n value for efficiently finding elements OF NORM 1. See the comment inside afuchfindonelt_i for how to change C when looking for non-norm 1 elements.*/
static GEN
afuchbestC(GEN A, GEN O, GEN Olevel_nofact, long prec)
{
  pari_sp av = avma;
  GEN F = alg_get_center(A);
  long n = nf_get_degree(F);
  GEN Adisc = algdiscnorm(A);/*Norm to Q of disc(A)*/
  if (!gequal1(Olevel_nofact)) Adisc = mulii(Adisc, idealnorm(F, Olevel_nofact));/*Incorporating the norm to Q of the level.*/
  GEN discpart = gmul(nf_get_disc(F), gsqrt(Adisc, prec));/*disc(F)*sqrt(Adisc)*/
  GEN discpartroot = gpow(discpart, gdivgs(gen_1, n), prec);/*discpart^(1/n)=disc(F)^(1/n)*algdisc^(1/2n)*/
  GEN npart;
  double npart_d[9] = {0, 2.5, 1.325, 1.21, 1.21, 1.35, 1.4, 1.44, 1.5};
  if (n <= 8) npart = gtofp(dbltor(npart_d[n]), prec);
  else npart = gtofp(gadd(dbltor(1.5), gmulsg(n - 8, dbltor(0.05))), prec);/*Seems to be a reasonably good choice, though hard to say for sure.*/
  GEN best = gerepileupto(av, gmul(npart, discpartroot));/*npart*disc(F)^(1/n)*N_F/Q(algebra disc)^(1/2n)*/
  if (gcmpgs(best, n) <= 0) best = gerepileupto(av, gaddsg(n, gen_2));/*Make sure best>n. If it is not, then we just add 2 (I doubt this will ever occur, but maybe in a super edge case).*/
  return best;
}

/*Generates the chosen constants used in computing the fundamental domain.*/
static GEN
afuchfdomdat_init(GEN A, GEN O, long prec)
{
  pari_sp av = avma;
  GEN F = alg_get_center(A);
  GEN Olevel_fact = algorderlevel(A, O, 1);
  GEN area = afuchO1area(A, O, Olevel_fact, 3, prec);
  GEN Olevel_nofact = idealfactorback(F, Olevel_fact, NULL, 0);
  GEN C = afuchbestC(A, O, Olevel_nofact, prec);
  GEN gamma = gtofp(dbltor(2.1), prec);
  GEN R = hdiscradius(gpow(area, gamma, prec), prec);/*Setting R*/
  GEN epsilon = mkfracss(1, 6), passes;
  if (nf_get_degree(alg_get_center(A)) == 1) passes = gen_2;
  else passes = stoi(8);
  return gerepilecopy(av, mkvec5(area, C, R, epsilon, passes));
}


/*3: ALGEBRA FUNDAMENTAL DOMAIN METHODS*/

/*Returns the algebra in X, which may have been computed to higher accuracy.*/
GEN
afuchalg(GEN X)
{
  return gcopy(afuch_get_alg(X));
}

/*Returns the area of X.*/
GEN
afucharea(GEN X)
{
  pari_sp av = avma;
  GEN U = afuch_get_fdom(X);
  if (!U) { afuchfdom(X); U = afuch_get_fdom(X); }
  return gerepilecopy(av, normbound_get_area(U));
}

/*Returns the elements giving the side for the fundamental domain of X.*/
GEN
afuchelts(GEN X)
{
  pari_sp av = avma;
  GEN U = afuch_get_fdom(X);
  if (!U) { afuchfdom(X); U = afuch_get_fdom(X); }
  GEN elts = normbound_get_elts(U);
  GEN O = afuch_get_O(X);
  if (gequal1(O)) return gerepilecopy(av, elts);
  long i, le;
  GEN newelts = cgetg_copy(elts, &le);
  for (i = 1; i < le; i++) gel(newelts, i) = QM_QC_mul(O, gel(elts, i));
  return gerepileupto(av, newelts);
}

/*Computes the fundamental domain for O^1, DEBUGLEVEL allows extra input to be displayed. Returns NULL if precision too low. Can pass in a starting set. This is useful in case we do some computations then have too low precision, we don't lose the computations.*/
static GEN
afuchfdom_i(GEN X, GEN *startingset)
{
  pari_sp av = avma;
  GEN gdat = afuch_get_gdat(X);
  long prec = lg(gdat_get_tol(gdat));
  GEN twopi = Pi2n(1, prec);
  GEN area = afuch_get_O1area(X);
  GEN areabound = addrr(area, shiftr(area, -1));/*area*1.5*/
  GEN C = afuch_get_bestC(X);/*Optimally chosen value of C*/
  GEN R = afuch_get_R(X);/*Starting radius for finding points. Will increase if not sufficient.*/
  GEN epsilon = afuch_get_epsilon(X);
  long n = nf_get_degree(alg_get_center(afuch_get_alg(X)));
  GEN passes = afuch_get_passes(X);/*Based on heuristics, we choose the number of points at each stage to expect to finish after this number of passes. Well, that would be true if the probability that we find a point that is there is 1, which it isn't since we are pruning.*/
  long N = 1 + itos(gceil(gdiv(gsqr(area), gmul(gmul(shiftr(twopi, 2), gsubgs(C, n)), passes))));/*Area^2/(8*Pi*(C-n)*#Passes)*/
  if (N < 3) N = 3;/*Make sure N>=2; I think it basically always should be, but this guarantees it.*/
  if (DEBUGLEVEL > 0) pari_printf("Initial constants:\n   C=%P.8f\n   N=%d\n   R=%P.8f\nTarget Area: %P.8f\n\n", C, N - 1, R, area);
  pari_sp av_mid = avma;
  GEN firstelts = *startingset;
  if (!firstelts) {/*No starting set passed.*/
    firstelts = afuchfindelts(X, gen_1, gtocr(gen_0, prec), C, 1, NULL, NULL);/*This may find an element with large radius, a good start.*/
    if (!firstelts) return gc_NULL(av);/*Precision too low.*/
  }
  GEN U = NULL;
  if (lg(firstelts) > 1) U = normbound(X, firstelts, &afuchtoklein, gdat);/*Start off the boundary.*/
  long nU = U ? lg(normbound_get_elts(U)) - 1 : 0;/*Number of sides.*/
  long pass = 0;
  for (;;) {
    pass++;
    if (DEBUGLEVEL > 0) pari_printf("Pass %d with %d random points in the ball of radius %P.8f\n", pass, N - 1, R);
    GEN elts = vectrunc_init(N);
    GEN oosides = nU ? normbound_get_infinite(U) : NULL;
    long noo = nU ? lg(oosides) - 1 : 0, i;
    if (noo) {/*Looking near infinite sides*/
      if (DEBUGLEVEL > 0) pari_printf("%d infinite sides\n", noo);
      long iside = 0;
      GEN Uvargs = normbound_get_vargs(U);
      for (i = 1; i < N; i++){
        iside++;
        if (iside > noo) iside = 1;
        long sind = oosides[iside];
        GEN ang2 = gel(Uvargs, sind), ang1;
        if (sind == 1) ang1 = gel(Uvargs, nU);
        else ang1 = gel(Uvargs, sind - 1);
        if (cmprr(ang1, ang2) > 0) ang2 = addrr(ang2, twopi);/*We loop past zero.*/
        GEN z = hdiscrandom_arc(R, ang1, ang2, prec);
        GEN found = afuchfindelts(X, gen_1, z, C, 1, NULL, NULL);
        if (!found) {/*Precision too low.*/
          *startingset = gerepilecopy(av, normbound_get_elts(U));
          return NULL;
        }
        if (lg(found) > 1) vectrunc_append(elts, gel(found, 1));/*Found an element!*/
      }
    }
    else {
      for (i = 1; i < N; i++){
        GEN z = hdiscrandom(R, prec);
        GEN found = afuchfindelts(X, gen_1, z, C, 1, NULL, NULL);
        if (!found) {/*Precision too low.*/
          if (!U) return gc_NULL(av);
          *startingset = gerepilecopy(av, normbound_get_elts(U));
          return NULL;
        }
        if (lg(found) > 1) vectrunc_append(elts, gel(found, 1));/*Found an element!*/
      }
    }
    if (DEBUGLEVEL > 0) pari_printf("%d elements found\n", lg(elts) - 1);
    if (lg(elts) > 1) U = normbasis(X, U, elts, &afuchtoklein, &afuchmul, &afuchconj, &afuchistriv, gdat);
    if (!U) {
      nU = 0;
      if (DEBUGLEVEL > 0) pari_printf("Current normalized basis has 0 sides\n\n");
      R = gadd(R, epsilon);/*Update R, maybe it was too small.*/
      gerepileall(av_mid, 1, &R);
      U = NULL;
      continue;
    }
    long newnU = lg(normbound_get_elts(U)) - 1;
    if (DEBUGLEVEL > 0) pari_printf("Current normalized basis has %d sides\n\n", newnU);
    GEN Uarea = normbound_get_area(U);
    if (gcmp(Uarea, areabound) < 0) break;/*Done*/
    if (pass > 1 && nU == newnU) R = gadd(R, epsilon);/*Updating R if we didn't change the number of sides.*/
    nU = newnU;
    if (gc_needed(av, 2)) {
      if (DEBUGLEVEL > 0) pari_printf("Garbage collection.\n");
      gerepileall(av_mid, 2, &U, &R);
    }
  }
  return gerepileupto(av, U);
}

/*Computes and stores the fundamental domain, recomputing X to more accuracy if the precision is too low.*/
void
afuchfdom(GEN X)
{
  pari_sp av = avma;
  GEN U = afuch_get_fdom(X);
  if (U) return;/*We already have U!*/
  GEN Gtype = afuch_get_type(X);
  if (typ(Gtype) != t_INT) pari_err_TYPE("Type should be 0, 1, 2, or 3", Gtype);
  long type = itos(Gtype);
  GEN allelts = obj_check(X, afuch_SAVEDELTS);
  if (allelts) {/*We already have a set of generators for everything, so just call normbasis on the appropriate thing.*/
    GEN S = gel(allelts, 1);/*Norm 1, incluced in everything*/
    switch (type) {
      case 3:
        S = shallowconcat(S, gel(allelts, 4));
      case 2:
        S = shallowconcat(S, gel(allelts, 3));
      case 1:
        S = shallowconcat(S, gel(allelts, 2));
    }
    U = normbasis(X, NULL, S, &afuchtoklein, &afuchmul, &afuchconj, &afuchistriv, afuch_get_gdat(X));
    obj_insert(X, afuch_FDOM, U);
    set_avma(av);
    return;
  }
  int precinc = 0;
  GEN startingset = NULL;
  for (;;) {
    U = afuchfdom_i(X, &startingset);
    if (U) break;/*Success!*/
    if (DEBUGLEVEL > 0) pari_warn(warner, "Increasing precision");
    precinc = 1;
    afuch_moreprec(X, 1);/*Increase the precision by 1.*/
    if (startingset) {/*We did find some elements, so we save them (after removing the 0's).*/
      long ls = lg(startingset), i;
      GEN newstart = vectrunc_init(ls);
      for (i = 1; i < ls; i++) {
        GEN elt = gel(startingset, i);
        if (!gequal0(elt)) vectrunc_append(newstart, elt);
      }
      startingset = newstart;
    }
  }
  if (precinc) {
    GEN tol = gdat_get_tol(afuch_get_gdat(X));
    if (DEBUGLEVEL > 0) pari_warn(warner, "Precision increased to %d, i.e. \\p%Pd", lg(tol), precision00(tol, NULL));
  }
  if (!type) {/*Looking for O^1 only.*/
    obj_insert(X, afuch_FDOM, U);
    set_avma(av);
    return;
  }
  GEN O1elts = normbound_get_elts(U), newelts, S;
  switch (type) {
    case 1:
      S = shallowconcat(O1elts, afuch_makeunitelts(X));
      break;
    case 2:
      newelts = afuch_makeALelts(X);
      S = shallowconcat(O1elts, shallowconcat1(newelts));
      break;
    case 3:/*case 3*/
      newelts = afuch_makenormelts(X);
      S = shallowconcat(O1elts, shallowconcat1(newelts));
      break;
    default:
      S = cgetg(1, t_VEC);/*In case we input a bad type value.*/
  }
  U = normbasis(X, NULL, S, &afuchtoklein, &afuchmul, &afuchconj, &afuchistriv, afuch_get_gdat(X));
  GEN elts = normbound_get_elts(U), kact = normbound_get_kact(U);
  long le = lg(elts), i;
  for (i = 1; i < le; i++) {/*Due to multiplication, we may have huge elements that can be cut down. We do this now.*/
    GEN den;
    gel(elts, i) = Q_primitive_part(gel(elts, i), &den);
    if (den) gel(kact, i) = afuchtoklein(X, gel(elts, i));/*Recompute this way in case of precision mishaps.*/
  }
  obj_insert(X, afuch_FDOM, U);
  if (type == 3) obj_insert(X, afuch_SAVEDELTS, vec_prepend(newelts, O1elts));
  set_avma(av);
}

/*Finds the image of the root geodesic of g in the fundamental domain. The return is a vector of [g, s1, s2, v1, v2, [a, b, c]], where each component runs from vertex v1 on side s1 to vertex v2 on side s2, which has equation ax+by=c=0 or 1. The components are listed in order.*/
GEN
afuchgeodesic(GEN X, GEN g)
{
  pari_sp av = avma;
  GEN gdat = afuch_get_gdat(X);
  GEN U = afuch_get_fdom(X);
  if (!U) { afuchfdom(X); U = afuch_get_fdom(X); }
  GEN O = afuch_get_O(X);
  int isO = !gequal1(O);
  if (isO) {/*Convert g to in O*/
    GEN Oinv = afuch_get_Oinv(X);
    g = QM_QC_mul(Oinv, g);
  }
  GEN geod = geodesic_fdom(X, U, g, afuchid(X), &afuchtoklein, &afuchmul, &afuchconj, gdat);
  if (!isO) return gerepileupto(av, geod);/*No conversion necessary*/
  long i, lgeo = lg(geod);
  for (i = 1; i < lgeo; i++) gmael(geod, i, 1) = QM_QC_mul(O, gmael(geod, i, 1));/*Convert back.*/
  return gerepilecopy(av, geod);
}

/*Given a totally real number field F (with variable not x), we return [pairs, areas, rprimes], where A=alginit(F, pairs[i]) gives an arithmetic Fuchsian group (with respect to the maximal order) whose area, areas[i], is between Amin and Amax, and the multiset of primes lying above the ramified ideals is rprimes[i]. In fact, we find all such algebras that are split only at the place "split". Can pass Amax as NULL to go from 0 to Amin. Currently, we do not treat Eichler orders here.*/
GEN
afuchlist(GEN F, GEN Amin, GEN Amax, long split)
{
  pari_sp av = avma, av2;
  long bits = bit_accuracy(3);
  GEN pol = nf_get_pol(F);
  GEN zetaval = lfun(pol, gen_2, bits);/*zeta_F(2)*/
  long n = nf_get_degree(F), twon = n << 1;
  GEN ar = gmul(gpow(nfdisc(pol), gdivsg(3, gen_2), 3), zetaval);/*zeta_F(2)*d_F^(3/2)*/
  ar = mpshift(ar, 3 - twon);
  ar = gmul(ar, gpowgs(mppi(3), 1 - twon));/*2^(3-2n)*Pi^(1-2n)*zeta_F(2)*d_F^(3/2). The area is ar*product of N(p)-1 across the prime ideals dividing D*/
  if (!Amax) { Amax = Amin; Amin = gen_0; }
  GEN phimin = gceil(gdiv(Amin, ar));
  GEN phimax = gfloor(gdiv(Amax, ar));/*The min and max possible values of the product of N(p)-1*/
  GEN plist = primes0(mkvec2(gen_2, addis(phimax, 1)));/*The possible primes lying above divisors of D in the given area range.*/
  long lp = lg(plist), i, j;
  GEN poss_ideal = vectrunc_init(n*lp), poss_norm = vectrunc_init(n*lp);
  for (i = 1; i < lp; i++) {
    GEN pdec = idealprimedec(F, gel(plist, i));
    for (j = 1; j < lg(pdec); j++) {
      GEN nm = subis(pr_norm(gel(pdec, j)), 1);
      if (gcmp(nm, phimax) <= 0) {/*It's possible!*/
        vectrunc_append(poss_ideal, gel(pdec, j));/*Add the ideal*/
        vectrunc_append(poss_norm, nm);/*Add the norm*/
      }
    }
  }
  GEN order = indexsort(poss_norm);/*Order smallest to largest.*/
  long par = (n - 1)%2, maxprimes = 0, leno = lg(order) - 1;/*The parity of the number of prime divisors we need.*/
  GEN pro = gen_1;
  while (cmpii(pro, phimax) <= 0 && maxprimes < leno){/*maxprimes keeps track of the maximal number of prime divisors we can take*/
    maxprimes++;
    pro = mulii(pro, gel(poss_norm, order[maxprimes]));
  }
  if (maxprimes%2 != par) maxprimes--;/*Make it line up with parity.*/
  GEN infram = const_vecsmall(n, 1);
  infram[split] = 0;/*The infinite ramification vector.*/
  GEN condition = cgetg(4, t_VEC);/*For specifying the algebra by Hasse invariants*/
  gel(condition, 1) = gen_2;
  gel(condition, 3) = infram;
  long maxalg = lp, foundalg = 0;
  GEN algdat = cgetg(maxalg + 1, t_VEC);/*Stores the data we find.*/
  long nprime;
  for (nprime = maxprimes; nprime >= 0; nprime = nprime - 2) {/*Look for all subsets with this number of primes that satisfies the condition.*/
    gel(condition, 2) = cgetg(3, t_VEC);
    gmael(condition, 2, 2) = const_vec(nprime, gen_1);/*All the finite places we specify should ramify.*/
    gmael(condition, 2, 1) = cgetg(nprime + 1, t_VEC);/*To store those primes.*/
    GEN S = cgetg(nprime + 1, t_VECSMALL);/*Look over nprime element subsets. I could do this more efficiently (in terms of failing a test => we fail the next set), but this is not the bottleneck and should be OK for now.*/
    for (i = 1; i <= nprime; i++) S[i] = i;
    for(;;) {
      av2 = avma;
      pro = gen_1;
      for (i = 1; i <= nprime; i++) pro = mulii(pro, gel(poss_norm, order[S[i]]));/*Find the total product.*/
      if (cmpii(pro, phimin) < 0) {
        set_avma(av2);
        if(nextsub(S, leno)) continue;/*To small, go on.*/
        break;/*We reached the end.*/
      }
      if (cmpii(pro, phimax) > 0) {/*To large! Let's skip.*/
        i = nprime;
        while (i > 1) {/*Finding the last block of consecutive numbers.*/
          if (S[i] - S[i - 1] != 1) break;
          i--;
        }
        if (i <= 1) break;/*S[1]...S[nprime] is consecutive and too big, so we will always lose from now on. Include i=0 in case there were never any primes (oo loop else).*/
        S[i - 1]++;
        for (j = i; j <= nprime; j++) S[j] = S[j - 1] + 1;
        set_avma(av2);
        continue;
      }
      for (i = 1; i <= nprime; i++) gmael3(condition, 2, 1, i) = gel(poss_ideal, order[S[i]]);/*Add these prime ideals to the ramification vector.*/
      GEN A = alginit(F, condition, -1, 1);
      GEN ab = algab(A);/*The pair [a, b]*/
      for (i = 1; i <= 2; i++) {
        GEN den = Q_denom(gel(ab, i));/*Currently, to use alginit(F, ab), ab needs to have no denominator.*/
        if (!equali1(den)) gel(ab, i) = gmul(gel(ab, i), sqri(den));/*Scale a/b*/
      }
      GEN curarea = gmul(ar, pro);/*The area*/
      GEN ramp = cgetg(nprime + 1, t_VEC);
      for (i = 1; i <= nprime; i++) gel(ramp, i) = pr_get_p(gel(poss_ideal, order[S[i]]));
      ramp = ZV_sort(ramp);
      GEN curdat = gerepilecopy(av2, mkvec3(ab, curarea, ramp));/*The data for this algebra.*/
      foundalg++;
      if (foundalg > maxalg) {/*Lengthen the vector*/
        maxalg <<= 1;/*Double length*/
        algdat = vec_lengthen(algdat, maxalg);/*Actually initalize it.*/
      }
      gel(algdat, foundalg) = curdat;
      if(nextsub(S, leno)) continue;/*To small, go on.*/
      break;/*We reached the end.*/
    }
  }
  setlg(algdat, foundalg + 1);/*Chop off the end.*/
  return gerepileupto(av, vecsort(algdat, mkvecsmall(2)));/*Sort by area*/
}

/*Essentially does forsubset_next to S, but by passing in the acutal Vecsmall we can skip large chunks by directly changing it. Returns 1 if we successfully changed it, 0 if we are done.*/
static int
nextsub(GEN S, long n)
{
  long i, cur = n, j;
  for (i = lg(S) - 1; i > 0; i--) {
    if (S[i] < cur) {
      S[i]++;
      for (j = i + 1; j < lg(S); j++) S[j] = S[j - 1] + 1;
      return 1;
    }
    cur--;
  }
  return 0;
}

/*Possible norms of normalizer elements.*/
GEN
afuchnormalizernorms(GEN X)
{
  pari_sp av = avma;
  GEN NN = afuch_get_normalizernorms(X);
  if (NN) return gcopy(NN);
  GEN A = afuch_get_alg(X);
  GEN F = alg_get_center(A);
  long prec = afuch_get_prec(X);
  GEN B = Buchall(F, 0, prec);
  GEN ram_disc = algramifiedplacesf(A), ramid;
  GEN O = afuch_get_O(X);
  if (gequal1(O)) ramid = ram_disc;
  else {/*Incorporate the level of O too*/
    GEN ram_level = algorderlevel(A, O, 1);
    long lr = lg(ram_disc), llev = lg(gel(ram_level, 1)), i;
    ramid = cgetg(lr + llev - 1, t_VEC);
    for (i = 1; i < lr; i++) gel(ramid, i) = gel(ram_disc, i);/*Copy these over*/
    for (i = 1; i < llev; i++) {
      GEN theid = idealpow(F, gcoeff(ram_level, i, 1), gcoeff(ram_level, i, 2));
      gel(ramid, i + lr - 1) = theid;
    }
  }
  GEN norms = normalizer_make_norms(B, algsplitoo(A), ramid, prec);
  obj_insert(X, afuch_NORMALIZERNORMS, norms);
  return gerepileupto(av, norms);
}

/*Stored order*/
GEN
afuchorder(GEN X)
{
  return gcopy(afuch_get_O(X));
}

/*Presentation*/
GEN
afuchpresentation(GEN X)
{
  pari_sp av = avma;
  GEN pres = afuch_get_pres(X);
  if (!pres) {
    GEN U = afuch_get_fdom(X);
    if (!U) { afuchfdom(X); U = afuch_get_fdom(X); }
    pres = presentation(X, U, afuchid(X), &afuchmul, &afuchisparabolic, &afuchistriv);
    obj_insert(X, afuch_PRES, pres);
  }
  GEN O = afuch_get_O(X);
  if (gequal1(O)) return gerepileupto(av, pres);/*O=1, no conversion necessary.*/
  GEN Opres = cgetg(4, t_VEC);
  long lgen, i;
  gel(Opres, 1) = cgetg_copy(gel(pres, 1), &lgen);
  for (i = 1; i < lgen; i++) gmael(Opres, 1, i) = QM_QC_mul(O, gmael(pres, 1, i));
  gel(Opres, 2) = gcopy(gel(pres, 2));
  gel(Opres, 3) = gcopy(gel(pres, 3));
  return gerepileupto(av, Opres);
}

/*Returns the sides of the fundamental domain.*/
GEN
afuchsides(GEN X)
{
  GEN U = afuch_get_fdom(X);
  if (!U) { afuchfdom(X); U = afuch_get_fdom(X); }
  GEN oldsides = normbound_get_sides(U);
  long lo, i, j;
  GEN sides = cgetg_copy(oldsides, &lo);
  for (i = 1; i < lo; i++) {
	gel(sides, i) = cgetg(4, t_VEC);
	for (j = 1; j <= 3; j++) gmael(sides, i, j) = gcopy(gmael(oldsides, i, j));
  }
  return sides;
}

/*Signature*/
GEN
afuchsignature(GEN X)
{
  pari_sp av = avma;
  GEN sig = afuch_get_sig(X);/*Already computed.*/
  if (sig) return sig;
  GEN U = afuch_get_fdom(X);
  if (!U) { afuchfdom(X); U = afuch_get_fdom(X); }
  sig = signature(X, U, afuchid(X), &afuchmul, &afuchisparabolic, &afuchistriv);
  obj_insert(X, afuch_SIG, sig);
  return gerepileupto(av, sig);
}

/*Returns the side pairing for X.*/
GEN
afuchspair(GEN X)
{
  GEN U = afuch_get_fdom(X);
  if (!U) { afuchfdom(X); U = afuch_get_fdom(X); }
  return gcopy(normbound_get_spair(U));
}

/*Writing an element as a word in the presentation. We do not reduce wrt the generators.*/
GEN
afuchword(GEN X, GEN g)
{
  pari_sp av = avma;
  GEN U = afuch_get_fdom(X);
  if (!U) { afuchfdom(X); U = afuch_get_fdom(X); }
  GEN P = afuch_get_pres(X);
  if (!P) {
    GEN pre = afuchpresentation(X);
    cgiv(pre);
    P = afuch_get_pres(X);/*The return value of afuchpresentation converts back to A from O, which we do not want.*/
  }
  GEN O = afuch_get_O(X);
  if (!gequal1(O)) {
    GEN Oinv = afuch_get_Oinv(X);
    g = QM_QC_mul(Oinv, g);/*Convert to in O.*/
  }
  return gerepileupto(av, word(X, U, P, g, &afuchtoklein, &afuchmul, &afuchconj, &afuchistriv, afuch_get_gdat(X)));
}

/*Returns the vertices of the fundamental domain in the Klein model (model=0), or unit disc model (model=1).*/
GEN
afuchvertices(GEN X, int model)
{
  GEN U = afuch_get_fdom(X);
  if (!U) { afuchfdom(X); U = afuch_get_fdom(X); }
  GEN verts = normbound_get_vcors(U);
  if (!model) return gcopy(verts);
  GEN tol = gdat_get_tol(afuch_get_gdat(X));
  long lv, i;
  GEN vdisc = cgetg_copy(verts, &lv);
  for (i = 1; i < lv; i++) gel(vdisc, i) = klein_to_disc(gel(verts, i), tol);
  return vdisc;
}


/*3: NON NORM 1 METHODS*/

/*Returns a vector of generators for F^{x}O^{unit norm}/F^{x}O^1, which is a (multiplicative) F_2 vector space. We can pass in unitnorms, which is the first entry of bnf_make_unitnorms.*/
static GEN
afuch_makeunitelts(GEN X)
{
  pari_sp av = avma;
  GEN A = afuch_get_alg(X);
  GEN F = alg_get_center(A);
  long prec = afuch_get_prec(X);
  GEN B = Buchall(F, 0, prec);
  GEN unitnorms = gel(bnf_make_unitnorms(B, algsplitoo(A), prec), 1);
  long lu, i;
  GEN elts = cgetg_copy(unitnorms, &lu);
  for (i = 1; i < lu; i++) gel(elts, i) = afuchfindoneelt_i(X, gel(unitnorms, i), NULL);
  return gerepileupto(av, elts);
}

/*Returns a vector of generators for AL(O)^+/F^{x}O^1, which is a (multiplicative) F_2 vector space.*/
static GEN
afuch_makeALelts(GEN X)
{
  pari_sp av = avma;
  GEN A = afuch_get_alg(X);
  GEN F = alg_get_center(A);
  long prec = afuch_get_prec(X);
  GEN B = Buchall(F, 0, prec);
  GEN ram_disc = algramifiedplacesf(A), ramid;
  GEN O = afuch_get_O(X);
  if (gequal1(O)) ramid = ram_disc;
  else {/*Incorporate the level of O too*/
    GEN ram_level = algorderlevel(A, O, 1);
    long lr = lg(ram_disc), llev = lg(gel(ram_level, 1)), i;
    ramid = cgetg(lr + llev - 1, t_VEC);
    for (i = 1; i < lr; i++) gel(ramid, i) = gel(ram_disc, i);/*Copy these over*/
    for (i = 1; i < llev; i++) {
      GEN theid = idealpow(F, gcoeff(ram_level, i, 1), gcoeff(ram_level, i, 2));
      gel(ramid, i + lr - 1) = theid;
    }
  }
  GEN ALnorms = AL_make_norms(B, algsplitoo(A), ramid, prec);
  GEN norms1 = gel(ALnorms, 1), norms2 = gel(ALnorms, 2);/*Unit, then AL norms*/
  long lu, i;
  GEN elts1 = cgetg_copy(norms1, &lu);
  for (i = 1; i < lu; i++) gel(elts1, i) = afuchfindoneelt_i(X, gel(norms1, i), NULL);
  GEN elts2 = cgetg_copy(norms2, &lu);
  for (i = 1; i < lu; i++) gel(elts2, i) = afuchfindoneelt_i(X, gel(norms2, i), NULL);
  return gerepilecopy(av, mkvec2(elts1, elts2));
}

/*Returns a vector of generators for N_{B^x}(O)^+/F^{x}O^1, which is a (multiplicative) F_2 vector space.*/
static GEN
afuch_makenormelts(GEN X)
{
  pari_sp av = avma;
  GEN norms = afuch_get_normalizernorms(X);
  if (!norms) {/*Must make the norms.*/
    GEN A = afuch_get_alg(X);
    GEN F = alg_get_center(A);
    long prec = afuch_get_prec(X);
    GEN B = Buchall(F, 0, prec);
    GEN ram_disc = algramifiedplacesf(A), ramid;
    GEN O = afuch_get_O(X);
    if (gequal1(O)) ramid = ram_disc;
    else {/*Incorporate the level of O too*/
      GEN ram_level = algorderlevel(A, O, 1);
      long lr = lg(ram_disc), llev = lg(gel(ram_level, 1)), i;
      ramid = cgetg(lr + llev - 1, t_VEC);
      for (i = 1; i < lr; i++) gel(ramid, i) = gel(ram_disc, i);/*Copy these over*/
      for (i = 1; i < llev; i++) {
        GEN theid = idealpow(F, gcoeff(ram_level, i, 1), gcoeff(ram_level, i, 2));
        gel(ramid, i + lr - 1) = theid;
      }
    }
    norms = normalizer_make_norms(B, algsplitoo(A), ramid, prec);
    obj_insert(X, afuch_NORMALIZERNORMS, norms);
  }
  GEN norms1 = gel(norms, 1), norms2 = gel(norms, 2), norms3 = gel(norms, 3);/*Unit, then AL norms, then normalizer norms*/
  long lu, i;
  GEN elts1 = cgetg_copy(norms1, &lu);
  for (i = 1; i < lu; i++) gel(elts1, i) = afuchfindoneelt_i(X, gel(norms1, i), NULL);
  GEN elts2 = cgetg_copy(norms2, &lu);
  for (i = 1; i < lu; i++) gel(elts2, i) = afuchfindoneelt_i(X, gel(norms2, i), NULL);
  GEN elts3 = cgetg_copy(norms3, &lu);
  for (i = 1; i < lu; i++) gel(elts3, i) = afuchfindoneelt_i(X, gel(norms3, i), NULL);
  return gerepilecopy(av, mkvec3(elts1, elts2, elts3));
}

/*
Let U be the unit group of B, which is a totally real number field. This method returns [v, uneg] where:
v is a basis for the totally positive units in U/U^2
uneg is a unit with sigma_i(nm)>0 if and only if i!=split, if it exists. If it does not exist, this entry is 0.
It is not possible for both uneg=0 and v to be empty, as v is empty <==> the generating basis of U is a basis of U/U^2 ==> uneg exists.
*/
GEN
bnf_make_unitnorms(GEN B, long split, long prec)
{
  pari_sp av = avma;
  GEN F = bnf_get_nf(B);
  GEN units = gel(bnfunits(B, NULL), 1);/*The n units: first n-1 are infinite order, last one is -1.*/
  long lu = lg(units), i, j;
  for (i = 1; i < lu; i++) gel(units, i) = lift(basistoalg(F, nffactorback(F, gel(units, i), NULL)));/*Now we have the units in algebraic form.*/
  GEN rts = nf_get_roots(F);
  GEN signs = cgetg(lu, t_MAT);/*Stores the signs of each unit at each place.*/
  for (i = 1; i < lu; i++) {
    GEN sgn = cgetg(lu, t_VECSMALL), un = gel(units, i);
    for (j = 1; j < lu; j++) sgn[j] = (1 - signe(poleval(un, gel(rts, j)))) >> 1;/*0 for sign 1, 1 for sign -1*/
    gel(signs, i) = sgn;
  }
  GEN ker = Flm_ker(signs, 2);/*Basis for the kernel, i.e. totally positive units.*/
  GEN vecei = vecsmall_ei(lu - 1, split);/*Positive outside of the split place, where negative.*/
  GEN left = Flm_Flc_invimage(signs, vecei, 2);/*Element giving this sign distribution, if it exists.*/
  GEN uneg;
  if (left) {
    uneg = gen_1;
    for (i = 1; i < lu; i++) if (left[i]) uneg = nfmul(F, uneg, gel(units, i));/*Multiply out the element.*/
    uneg = lift(basistoalg(F, uneg));
  }
  else uneg = gen_0;/*All units that are positive outside of split are positive at split too.*/
  long lv = lg(ker);
  GEN v = cgetg(lv, t_VEC);
  for (i = 1; i < lv; i++) {
    GEN elt = gen_1;
    GEN pat = gel(ker, i);
    for (j = 1; j < lu; j++) if (pat[j]) elt = nfmul(F, elt, gel(units, j));/*Multiply out the element.*/
    gel(v, i) = lift(basistoalg(F, elt));
  }
  return gerepilecopy(av, mkvec2(v, uneg));
}

/*
We make the possible norms for Atkin-Lehner elements, where B=bnfinit(F). Returns [vunit, vAL, uneg], as in bnf_make_unitnorms (if this did not have a uneg, then we might get one from here). ideals should be the vector of maximal prime power ideals dividing the reduced norm of the order (assuming Eichler order). We can also pass in B as the bnr initialized with respect to [1, [1,...,1,0,1,...,1]], where the 0 occurs at the unique split place.
*/
static GEN
AL_make_norms(GEN B, long split, GEN ideals, long prec)
{
  pari_sp av = avma;
  GEN F, C;
  long nF;
  if (nftyp(B) == typ_BNR) {/*We already made the BNR*/
    C = B;
    B = bnr_get_bnf(C);
    F = bnf_get_nf(B);
    nF = nf_get_degree(F);
  }
  else {
    F = bnf_get_nf(B);
    nF = nf_get_degree(F);
    GEN modulus = mkvec2(gen_1, const_vec(nF, gen_1));
    gmael(modulus, 2, split) = gen_0;/*Making the modulus. Norms of elements from our algebra must be positive at all non-split places, hence this choice.*/
    C = Buchray(B, modulus, nf_INIT);/*Compute the ray class field without generators.*/
  }
  GEN unitnorms = bnf_make_unitnorms(B, split, prec);/*We need this too, and this is why we saved a uneg if it existed.*/
  GEN orders = gel(bnr_get_clgp(C), 2);/*Odd orders we don't care about, only even, so let's modify this.*/
  long lo = lg(orders), i;
  GEN ordmod = cgetg(lo, t_VECSMALL);
  for (i = 1; i < lo; i++) {
    if (mod2(gel(orders, i))) ordmod[i] = 1;/*Everything is a square for this generator.*/
    else ordmod[i] = 2;
  }
  long li = lg(ideals), j;
  GEN alphas = cgetg_copy(ideals, &li);
  GEN mat = cgetg(li, t_MAT);/*Each of our ideals is alpha*prod(g_i^e_i) in the ray class group C. We store the alphas and the e_i's mod ordmod[i] since we care about whether we are in C^2 or not.*/
  for (i = 1; i < li; i++) {
    GEN prin = bnrisprincipal(C, gel(ideals, i), 1);/*Compute [ei, alpha]*/
    gel(alphas, i) = gel(prin, 2);/*Store alpha*/
    GEN col = cgetg(lo, t_VECSMALL);
    for (j = 1; j < lo; j++) col[j] = smodis(gmael(prin, 1, j), ordmod[j]);
    gel(mat, i) = col;/*We took the exponents modulo ordmod[j].*/
  }/*A product of the ideals is in C^2 if and only if the corresponding linear combination of the columns of mat is 0 in F_2.*/
  GEN ker = Flm_ker(mat, 2);/*A basis for the ideals.*/
  long lk = lg(ker);
  GEN newalphas = vectrunc_init(lk);
  for (i = 1; i < lk; i++) {
    GEN pattern = gel(ker, i);
    GEN alph = gen_1;
    for (j = 1; j < li; j++) {
      if (pattern[j] == 1) alph = nfmul(F, alph, gel(alphas, j));
    }
    GEN den;
    alph = Q_primitive_part(alph, &den);/*Scaling to be primitive and integral, which won't affect the signs in the embeddings.*/
    if (den) {/*We are actually only allowed to multiply by squares, so must fix this if we did not.*/
      GEN denfact = Q_factor(den);
      for (j = 1; j < lg(gel(denfact, 2)); j++) {
        if (mod2(gcoeff(denfact, j, 2))) alph = gmul(alph, gcoeff(denfact, j, 1));
      }
    }
    alph = lift(basistoalg(F, alph));
    if (gequal1(alph)) continue;/*It is possible to end up with 1, e.g. F=nfinit(y^3 - 104052*y - 12520924), split=1, ideals=idealprimedec(F, 5). Maybe there is a better way to handle this in general? Not sure if alph=1 is the only relation possible.*/
    vectrunc_append(newalphas, alph);
  }
  lk = lg(newalphas);
  /*We are almost there! We just have to modify the new found elements to be positive at the split place as well.*/
  GEN rt = gel(nf_get_roots(F), split);
  GEN swapper = gel(unitnorms, 2);/*This is -1 at the split place, if it exists.*/
  GEN ALnorms = vectrunc_init(lk);/*We either get this many, or one less if swapper=0.*/
  for (i = 1; i < lk; i++) {
    GEN nm = gel(newalphas, i);
    if (signe(poleval(nm, rt)) > 0) {
      vectrunc_append(ALnorms, nm);
      continue;/*All good!*/
    }
    if (gequal0(swapper)) {
      swapper = nm;
      continue;/*Use this element to swap the rest out.*/
    }
    nm = nfmul(F, nm, swapper);
    vectrunc_append(ALnorms, lift(basistoalg(F, nm)));/*Add it in.*/
  }
  return gerepilecopy(av, mkvec3(gel(unitnorms, 1), ALnorms, swapper));
}

/*
We make the possible norms for all normalizer elements. Returns [vunit, vAL, vclorder2, uneg], as in AL_make_unitnorms (if this did not have a uneg, then we might get one from here). ideals should be the maximal prime power ideals dividing the reduced norm of the order (assuming Eichler order).
*/
static GEN
normalizer_make_norms(GEN B, long split, GEN ideals, long prec)
{
  pari_sp av = avma;
  GEN F = bnf_get_nf(B);
  long nF = nf_get_degree(F);
  GEN modulus = mkvec2(gen_1, const_vec(nF, gen_1));
  gmael(modulus, 2, split) = gen_0;/*Making the modulus. Norms of elements from our algebra must be positive at all non-split places, hence this choice.*/
  GEN C = Buchray(B, modulus, nf_INIT | nf_GEN);/*Compute the ray class field with generators.*/
  GEN ALnorms = AL_make_norms(C, split, ideals, prec);/*We need this too, and this is why we saved a uneg if it existed.*/
  GEN orders = gel(bnr_get_clgp(C), 2);/*We want generators for the 2-order part of the class group (and do not care about 1).*/
  GEN gens = bnr_get_gen_nocheck(C);
  long lo = lg(orders), i;
  GEN ord2gens = vectrunc_init(lo);
  for (i = 1; i < lo; i++) {
    if (!mod2(gel(orders, i))) {
      GEN g = gel(gens, i);
      vectrunc_append(ord2gens, idealpows(F, g, itos(gel(orders, i)) >> 1));/*Save the order 2 generators.*/
    }
  }
  lo = lg(ord2gens);/*The actual generators.*/
  GEN Borders = gel(bnf_get_clgp(B), 2);/*Orders of the generators of B.*/
  long lBo = lg(Borders), j;
  GEN imCl = cgetg(lo, t_MAT);/*Find the image in Cl(R)[2], as we only keep a generating set for the survivors.*/
  for (i = 1; i < lo; i++) {
    GEN Clim = bnfisprincipal0(B, gel(ord2gens, i), 0);
    for (j = 1; j < lBo; j++) if (!gequal0(gel(Clim, j))) gel(Clim, j) = gen_1;/*We only care about 1 vs 0 for picking up part of the 2-torsion.*/
    gel(imCl, i) = Clim;
  }
  GEN left = gel(FpM_indexrank(imCl, gen_2), 2);/*Which columns we keep, i.e. do not get destroyed boosting up to CL(R)[2].*/
  long lgleft = lg(left);
  GEN alphs = cgetg(lgleft, t_VEC);
  for (i = 1; i < lgleft; i++) {
    GEN csqr = idealsqr(F, gel(ord2gens, left[i]));/*Principal*/
    GEN alph = gel(bnrisprincipal(C, csqr, 1), 2);/*Compute [e, alpha], where e=0.*/
    GEN den;
    alph = Q_primitive_part(alph, &den);/*Scaling to be primitive and integral, which won't affect the signs in the embeddings.*/
    if (den) {/*We are actually only allowed to multiply by squares, so must fix this if we did not.*/
      GEN denfact = Q_factor(den);
      for (j = 1; j < lg(gel(denfact, 2)); j++) {
        if (mod2(gcoeff(denfact, j, 2))) alph = gmul(alph, gcoeff(denfact, j, 1));
      }
    }
    gel(alphs, i) = lift(basistoalg(F, alph));/*Store the scaled alph*/
  }
  /*We are almost there! We just have to modify the new found elements to be positive at the split place as well.*/
  GEN rt = gel(nf_get_roots(F), split);
  GEN swapper = gel(ALnorms, 3);/*This is -1 at the split place, if it exists.*/
  GEN normalizernorms = vectrunc_init(lgleft);/*We either get this many, or one less if swapper=0.*/
  for (i = 1; i < lgleft; i++) {
    GEN nm = gel(alphs, i);
    if (signe(poleval(nm, rt)) > 0) {
      vectrunc_append(normalizernorms, nm);
      continue;/*All good!*/
    }
    if (gequal0(swapper)) {
      swapper = nm;
      continue;/*Use this element to swap the rest out.*/
    }
    nm = nfmul(F, nm, swapper);
    vectrunc_append(normalizernorms, lift(basistoalg(F, nm)));/*Add it in.*/
  }
  return gerepilecopy(av, mkvec3(gel(ALnorms, 1), gel(ALnorms, 2), normalizernorms));
}


/*3: ALGEBRA BASIC AUXILLARY METHODS*/

/*Conjugation formatted for the input of an afuch, for use in the geometry section. Since we work in O, entries are all in Z. We use this in place of alginv as well, to avoid a potential costly scaling.*/
static GEN
afuchconj(GEN X, GEN g)
{
  return ZM_ZC_mul(afuch_get_Oconj(X), g);
}

/*Returns the identity element*/
static GEN
afuchid(GEN X)
{
  return col_ei(lg(alg_get_tracebasis(afuch_get_alg(X))) - 1, 1);
}

/*Returns 1 if the element g (in O) is in the normalizer of O, 0 if not.*/
static int
afuchinnormalizer(GEN X, GEN g)
{
  pari_sp av = avma;
  long lO = lg(afuch_get_O(X)), i;
  GEN O1 = cgetg(lO, t_MAT), O2 = cgetg(lO, t_MAT);
  GEN v = col_ei(lO - 1, 1);
  gel(O1, 1) = g, gel(O2, 1) = g;/*O1=gO, O2=Og*/
  for (i = 2; i < lO; i++) {/*Form the conjugate matrix*/
    gel(v, i - 1) = gen_0;
    gel(v, i) = gen_1;
    gel(O1, i) = afuchmul(X, g, v);
    gel(O2, i) = afuchmul(X, v, g);
  }
  return gc_int(av, gequal(hnf(O1), hnf(O2)));
}

/*Returns 1 if g is parabolic, i.e. trd(g)^2=+/-2nrd(g), and 0 else.*/
static int
afuchisparabolic(GEN X, GEN g)
{
  pari_sp av = avma;
  GEN nm = afuchnorm_fast(X, g);
  GEN tr = afuchtrace(X, g);
  GEN F = alg_get_center(afuch_get_alg(X));
  GEN twonm = nfmul(F, nm, stoi(4));/*2*norm*/
  GEN trsq = nfsqr(F, tr);
  if (gequal(twonm, trsq) || gequal(twonm, gneg(trsq))) return gc_int(av, 1);
  return gc_int(av, 0);
}

/*Returns 1 if g is a scalar, i.e. g==conj(g).*/
static int
afuchistriv(GEN X, GEN g)
{
  pari_sp av = avma;
  GEN gconj = afuchconj(X, g);
  return gc_int(av, ZV_equal(g, gconj));
}

/*algmul for elements of O written in terms of the basis of O.*/
static GEN
afuchmul(GEN X, GEN g1, GEN g2)
{
  pari_sp av = avma;
  GEN T = afuch_get_Omultable(X);
  GEN g = ZC_Z_mul(ZM_ZC_mul(gel(T, 1), g2), gel(g1, 1));
  long lg1 = lg(g1), i;
  for (i = 2; i < lg1; i++) {
    g = ZC_add(g, ZC_Z_mul(ZM_ZC_mul(gel(T, i), g2), gel(g1, i)));
  }
  return gerepileupto(av, g);
}

/*Fast algnorm for elements written in terms of the basis of O.*/
static GEN
afuchnorm_fast(GEN X, GEN g)
{
  GEN F = alg_get_center(afuch_get_alg(X));
  GEN Onormdat = afuch_get_Onormdat(X);
  if (typ(Onormdat) == t_VEC) return afuchnorm_chol(F, Onormdat, g);
  return afuchnorm_mat(F, Onormdat, g);
}

/*If X has been initialized with Onorm_makechol, then this computes the norm of an element of O in basis form.*/
static GEN
afuchnorm_chol(GEN F, GEN chol, GEN g)
{
  pari_sp av = avma;
  GEN inds = gel(chol, 1);
  GEN coefs = gel(chol, 2);
  long n = lg(g), i, j;
  GEN s = gen_0;
  for (i = 1; i <= 4; i++) {
    GEN t = gel(g, inds[i]);
    for (j = inds[i] + 1; j < n; j++) t = nfadd(F, t, nfmul(F, gmael(coefs, i, j), gel(g, j)));
    t = nfsqr(F, t);
    s = nfadd(F, s, nfmul(F, gmael(coefs, i, inds[i]), t));
  }
  return gerepileupto(av, s);
}

/*If X has been initialized with Onorm_makemat, then this computes the norm of an element of O in basis form.*/
static GEN
afuchnorm_mat(GEN F, GEN Onorm, GEN g)
{
  pari_sp av = avma;
  long lg, i, j;
  GEN twog = cgetg_copy(g, &lg);
  for (i = 1; i < lg - 1; i++) gel(twog, i) = shifti(gel(g, i), 1);/*2g, except the last term which is uninitialized as we don't need it.*/
  GEN s = gen_0;
  for (i = 1; i < lg; i++) {
    GEN t = gen_0;
    for (j = 1; j < i; j++) t = nfadd(F, t, nfmul(F, gcoeff(Onorm, j, i), gel(twog, j)));
    t = nfadd(F, t, nfmul(F, gcoeff(Onorm, i, i), gel(g, i)));
    t = nfmul(F, t, gel(g, i));
    s = nfadd(F, s, t);
  }
  return gerepileupto(av, s);
}

/*Fast real approximation to algnorm.*/
static GEN
afuchnorm_real(GEN X, GEN g)
{
  return qfeval(afuch_get_Onormreal(X), g);
}

/*Given an element g of O (of non-zero norm) written in basis form, this returns the image of g in acting on the Klein model.*/
static GEN
afuchtoklein(GEN X, GEN g)
{
  pari_sp av = avma;
  GEN embs = afuch_get_kleinmats(X);
  GEN A = mulcri(gmael(embs, 1, 1), gel(g, 1));
  GEN B = mulcri(gmael(embs, 1, 2), gel(g, 2));
  long lg = lg(g), i;
  for (i = 2; i < lg; i++) {
    A = addcrcr(A, mulcri(gmael(embs, i, 1), gel(g, i)));
    B = addcrcr(B, mulcri(gmael(embs, i, 2), gel(g, i)));
  }
  return gerepilecopy(av, mkvec2(A, B));
}

/*Returns the trace of an element of X. We only use it in afuchisparabolic (which itself is only used in computing minimal cycles), so this is not as efficient as it could be (by saving the traces of a basis).*/
static GEN
afuchtrace(GEN X, GEN g)
{
  pari_sp av = avma;
  GEN Ag = QM_QC_mul(afuch_get_O(X), g);/*g in A*/
  return gerepileupto(av, algtrace(afuch_get_alg(X), Ag, 0));
}


/*3: FINDING ELEMENTS*/

/*Returns the quadratic form Q_{z, 0}^nm(g). If g has norm nm and is in the normalizer of O, then Q_{z, 0}(g)=cosh(d(gz, 0))+n-1 (n=deg(F), F is the centre of A). We output the symmetric matrix M, such that g~*M*g=Q_{z, 0}^nm(g). If nm!=1, then we can pass in the upper triangular part for Tr_{F/Q}(nrd(g)/nm) as tracepart, since this will be reused quite a bit (and we only save it for nm=1 in X). If passed as NULL, we will compute it here. We can also pass in the real value of nm (when nm!=1), or we can recompute it here if NULL.*/
static GEN
afuch_make_qf(GEN X, GEN nm, GEN z, GEN tracepart, GEN realnm)
{
  pari_sp av = avma;
  GEN qfmats = afuch_get_qfmats(X);
  GEN znorm = normcr(z);
  GEN rt = sqrtr(subsr(1, znorm));/*sqrt(1-|z|^2)*/
  GEN c = addrs(invr(rt), 1);/*1/sqrt(1-|z|^2)+1*/
  GEN zscale = invr(addrs(rt, 1));/*1/(1+sqrt(1-|z|^2))*/
  GEN w1 = mulrr(gel(z, 1), zscale), w2 = mulrr(gel(z, 2), zscale);/*w1+w2*i=w=z/(1+sqrt(1-|z|^2))*/
  GEN sumsq = addrr(sqrr(w1), sqrr(w2));/*w1^2+w2^2*/
  GEN N1 = rM_upper_r_mul(gel(qfmats, 1), sumsq);/*Msq*(w1^2+w2^2)*/
  GEN N2 = rM_upper_r_mul(gel(qfmats, 2), w1);/*M1*w1*/
  GEN N3 = rM_upper_r_mul(gel(qfmats, 3), w2);/*M2*w2*/
  GEN N = rM_upper_add(rM_upper_add(N1, N2), rM_upper_add(N3, gel(qfmats, 4)));/*N1+N2+N3+Mconst*/
  GEN qfup;
  if (gequal1(nm)) qfup = RgM_upper_add(rM_upper_r_mul(N, c), gel(qfmats, 5));/*The correct qf, in upper triangular form.*/
  else {/*Norm not 1, so we need to adjust qfmats[5].*/
    if (!tracepart) tracepart = afuch_make_traceqf(X, nm, NULL);/*Find Tr_{F/Q}(nrd(g)/nm).*/
    if (!realnm) {
      GEN A = afuch_get_alg(X), F = alg_get_center(A);
      long split = algsplitoo(A);
      long Fvar = nf_get_varn(F);
      GEN rt = gel(nf_get_roots(F), split);
      realnm = gsubst(nm, Fvar, rt);
    }
    c = mpdiv(c, realnm);
    qfup = RgM_upper_add(rM_upper_r_mul(N, c), tracepart);/*The correct qf, in upper triangular form.*/
  }
  long n = lg(qfup), i, j;
  for (i = 2; i < n; i++) for (j = 1; j < i; j++) gcoeff(qfup, i, j) = gcoeff(qfup, j, i);/*Make it symmetric, which is required for qfminim.*/
  return gerepilecopy(av, qfup);
}

/*Finds the norm nm elements of the normalizer such that Q_{z, 0}^nm(g)<=C. If maxelts is positive, this is the most number of returned elements. We will pick up elements g for which gz is close to 0. We use pruning, so there is no guarantee that we find elements that are there (but overall, it is quicker to find elements this way by repeating with different values of z). NOTE: we find element in the basis of O, NOT back to A. Can supply tracepart/realnm to use in afuch_make_qf (or pass as NULL).*/
static GEN
afuchfindelts(GEN X, GEN nm, GEN z, GEN C, long maxelts, GEN tracepart, GEN realnm)
{
  pari_sp av = avma;
  GEN tol = gdat_get_tol(afuch_get_gdat(X));
  long prec = lg(tol);
  GEN lowtol = deflowtol(prec);
  GEN M = afuch_make_qf(X, nm, z, tracepart, realnm);
  GEN v;
  if (gequal1(nm)) v = fincke_pohst_prune(M, C, 1, prec);
  else v = fincke_pohst_prune(M, C, 0, prec);/*ONLY prune if norm is 1.*/
  if(!v) return gc_NULL(av);/*Occurs when the precision is too low. Return NULL and maybe recompute above.*/
  int nomax, ind = 1, lv = lg(v), i;
  if (maxelts) nomax = 0;
  else { nomax = 1; maxelts = 10; }
  GEN found = cgetg(maxelts + 1, t_VEC);
  GEN A = afuch_get_alg(X);
  GEN F = alg_get_center(A);
  int checknormalizer = !gequal1(nm);/*If the norm is 1, no need to check the normalizer.*/
  if (nf_get_degree(F) == 1) {/*Over Q, so no real approximation required.*/
    for (i = 1; i < lv; i++) {
      GEN elt = gel(v, i);
      GEN norm = afuchnorm_fast(X, elt);
      if (!gequal(norm, nm)) continue;/*The norm was not nm.*/
      if (afuchistriv(X, elt)) continue;/*Ignore trivial elements.*/
      if (checknormalizer && !afuchinnormalizer(X, elt)) continue;/*Correct norm, but not in normalizer.*/
      gel(found, ind) = elt;
      ind++;
      if (ind <= maxelts) continue;/*We can find more*/
      if (!nomax) break;/*We have it the maximum return vector.*/
      maxelts <<= 1;/*Double and continue.*/
      found = vec_lengthen(found, maxelts);
    }  
  }
  else {
    if (!realnm) {/*Compute the real norm first.*/
      if (gequal1(nm)) realnm = gen_1;
      else {
        long split = algsplitoo(A);
        long Fvar = nf_get_varn(F);
        GEN rt = gel(nf_get_roots(F), split);
        realnm = gsubst(nm, Fvar, rt);
      }
    }
    for (i = 1; i < lv; i++) {
      GEN elt = gel(v, i);
      GEN realnorm = afuchnorm_real(X, elt);
      if (!toleq(realnorm, realnm, lowtol)) continue;/*Norm not nm up to tolerance.*/
      if (afuchistriv(X, elt)) continue;/*Ignore trivial elements.*/
      GEN norm = afuchnorm_fast(X, elt);
      norm = basistoalg(F, norm);
      if (!gequal(norm, nm)) continue;/*The norm was close to nm but not equal.*/
      if (checknormalizer && !afuchinnormalizer(X, elt)) continue;/*Correct norm, but not in normalizer.*/
      gel(found, ind) = elt;
      ind++;
      if (ind <= maxelts) continue;/*We can find more*/
      if (!nomax) break;/*We have it the maximum return vector.*/
      maxelts <<= 1;/*Double and continue.*/
      found = vec_lengthen(found, maxelts);
    }
  }
  found = vec_shorten(found, ind - 1);
  return gerepilecopy(av, found);
}

/*Finds and returns one non-trivial element of the positive normalizer of O with the given norm nm, by searching around random points.*/
static GEN
afuchfindoneelt_i(GEN X, GEN nm, GEN C)
{
  pari_sp av = avma, av2;
  long prec = afuch_get_prec(X);
  GEN R = afuch_get_R(X);
  GEN tracepart = NULL, realnm = NULL;
  if (!gequal1(nm)) {/*Set the trace part and real norm here, AND update C if not set.*/
    tracepart = afuch_make_traceqf(X, nm, NULL);/*Find Tr_{F/Q}(nrd(g)/nm).*/
    GEN A = afuch_get_alg(X), F = alg_get_center(A);
    long split = algsplitoo(A);
    long Fvar = nf_get_varn(F);
    GEN rt = gel(nf_get_roots(F), split);
    realnm = gsubst(nm, Fvar, rt);
	GEN nmnorm = nfnorm(F, nm);/*Norm to Q*/
	if (!C) {
	  C = afuch_get_bestC(X);
	  if (!equali1(nmnorm)) {/*We must scale it.*/
		/*
		Let L=nmnorm^(1/n). It turns out that the quadratic form we make, Q^nm(g)=Q_1(g)/sigma(nm)+Tr_{F/Q}(nrd(g)/nm), where sigma is the unique split oo place has the property that Q^nm<=C has the same number of solutions as Q^1<=CL. Assuming the time to compute this enumeration is A+B*(CL)^2n, we can do the same analysis as in my paper to get that
	      (2n-1)C^{2n}-2n^2C^{2n-1} = 1/L^{2n}(A/B)=1/nmnorm^2*A/B.
		Since we know this holds true with L=1 and we have the optimal value of C for this case, we can find A/B in this way, and solve for the new optimal C. Note that if L is smallish, this just changes C_opt to C_opt/L. If L gets larger, then we cannot ignore the second term and this no longer holds.
		*/
		long var = fetch_var(), n = nf_get_degree(F), twon = n << 1;
		long firstc = twon - 1, secondc = twon * n, i;
		GEN AB = gmul(gsubgs(gmulsg(firstc, C), secondc), gpowgs(C, firstc));/*(2n-1)C^(2n)-2n^2C^(2n-1)*/
		AB = gdiv(AB, sqri(nmnorm));/*Update the output, re-solve for C.*/
		GEN P = cgetg(twon + 3, t_POL);
		P[1] = evalvarn(var);
		gel(P, 2) = gneg(AB);/*Constant coefficient.*/
        for (i = 3; i <= twon; i++) gel(P, i) = gen_0;/*0's in the middle. P[i] corresp to x^(i-2)*/
        gel(P, twon + 1) = stoi(-secondc);
        gel(P, twon + 2) = stoi(firstc);
        P = normalizepol(P);/*Normalize it in place.*/
		GEN rts = realroots(P, mkvec2(stoi(n), mkoo()), prec);/*Real roots in [n, oo)*/
        delete_var();
        if (lg(rts) != 2) pari_err(e_MISC, "Could not find exactly one real root in [n, oo)!");
        C = gel(rts, 1);/*Finally, found it!*/
	  }
	}
  }
  else if (!C) C = afuch_get_bestC(X);
  long skip = 0;
  for (;;) {
    av2 = avma;
    GEN z = hdiscrandom(R, prec);
    GEN E = afuchfindelts(X, nm, z, C, 1, tracepart, realnm);
    if (!E) {/*Precision too low, we allow up to 20 exceptions before failing.*/
      skip++;
      if (skip > 20) pari_err_PREC("Precision too low for Fincke Pohst");
      continue;
    }
    if (lg(E) > 1) return gerepilecopy(av, gel(E, 1));
    set_avma(av2);
  }
}

/*afuchfindoneelt_i, except we convert the element back to being in A if it was in the order.*/
GEN
afuchfindoneelt(GEN X, GEN nm, GEN C)
{
  pari_sp av = avma;
  GEN O = afuch_get_O(X);
  GEN elt = afuchfindoneelt_i(X, nm, C);
  if (gequal1(O)) return elt;
  return gerepileupto(av, QM_QC_mul(O, elt));
}


/*3: ALGEBRA HELPER METHODS*/

/*Returns [a, b] where A=(a, b/F).*/
GEN
algab(GEN A)
{
  pari_sp av = avma;
  GEN L = alg_get_splittingfield(A), pol = rnf_get_pol(L);/*L=K(sqrt(a))*/
  long varn = rnf_get_varn(L);
  GEN a = gneg(gsubst(pol, varn, gen_0));/*Polynomial is of the form x^2-a, so plug in 0 and negate to get a*/
  GEN b = lift(alg_get_b(A));
  return gerepilecopy(av, mkvec2(a, b));
}

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

/*Given x=[e, f, g, h], returns x in the algebra representation.*/
GEN
alg1ijktoalg(GEN A, GEN x)
{
  pari_sp av = avma;
  GEN L = alg_get_splittingfield(A);/*L=F(i),*/
  long Lvar = rnf_get_varn(L);
  GEN Lx = pol_x(Lvar);/*x*/
  GEN e1 = gadd(gel(x, 1), gmul(gel(x, 2), Lx));/*e+fi*/
  GEN e2 = gsub(gel(x, 3), gmul(gel(x, 4), Lx));/*g-hi*/
  return gerepilecopy(av, mkcol2(e1, e2));
}

/*Given x=[e, f, g, h], returns x in the basis representation.*/
GEN
alg1ijktobasis(GEN A, GEN x)
{
  pari_sp av = avma;
  GEN xalg = alg1ijktoalg(A, x);
  return gerepileupto(av, algalgtobasis(A, xalg));
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

/*Returns the norm to Q of the discriminant of A*/
GEN
algdiscnorm(GEN A)
{
  pari_sp av=avma;
  GEN F = alg_get_center(A);/*Field*/
  GEN rams = algramifiedplacesf(A);
  GEN algdisc = gen_1;
  long lr = lg(rams), i;
  for (i = 1; i < lr; i++) algdisc = mulii(algdisc, idealnorm(F, gel(rams, i)));/*Norm to Q of the ramification*/
  return gerepileupto(av, algdisc);
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
voidalgmul(void *A, GEN x, GEN y)
{
  return algmul(*((GEN*)A), x, y);
}

/*Returns the vector of finite ramified places of the algebra A.*/
static GEN
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

/*Returns the reduced discriminant of A as an ideal in the centre.*/
GEN
algreduceddisc(GEN A)
{
  pari_sp av = avma;
  GEN F = alg_get_center(A);
  GEN rplaces = algramifiedplacesf(A);
  long l = lg(rplaces);
  GEN e = const_vecsmall(l - 1, 1);/*Vecsmall of 1's*/
  return gerepileupto(av, idealfactorback(F, rplaces, e, 0));
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


/*3: ALGEBRA ORDER METHODS*/

/*Given a1, ..., a4 in A, this returns d(a1,...,a4), i.e. det(trd(a_i*a_j))*/
static GEN
algd(GEN A, GEN a)
{
  pari_sp av = avma;
  long i, j;
  GEN M = cgetg(5, t_MAT);
  for (i = 1; i < 5; i++) gel(M, i) = cgetg(5, t_COL);
  for (i = 1; i < 5; i++) {/*Setting the coefficients of M*/
    gcoeff(M, i, i)=algtrace(A, algsqr(A, gel(a, i)), 0);
    for (j = i + 1; j < 5; j++) gcoeff(M, i, j) = gcoeff(M, j, i) = algtrace(A, algmul(A, gel(a, i), gel(a, j)), 0);
  }
  return gerepileupto(av, nfM_det(alg_get_center(A), M));
}

/*Checks if O is an order.*/
int
algisorder(GEN A, GEN O)
{
  pari_sp av = avma;
  GEN Oinv = QM_inv(O);
  long n = lg(O), i, j;
  for (i = 1; i < n; i++) {
    for (j = 1; j < n; j++) {
      GEN elt = algmul(A, gel(O, i), gel(O, j));
      if(!RgV_is_ZV(gmul(Oinv, elt))) return gc_int(av, 0);
    }
  }
  return gc_int(av, 1);
}

/*Given a vector of quaternionic elements Oalg in algebraic form given a basis for an order, we convert it to a matrix, whose columns express each basis element in terms of the natural order of A.*/
GEN
algorderalgtoorder(GEN A, GEN Oalg)
{
  long lO = lg(Oalg), i;
  GEN M = cgetg(lO, t_MAT);
  for (i = 1; i < lO; i++) gel(M, i) = algalgtobasis(A, gel(Oalg, i));
  return M;
}

/*Given an order O as a matrix of columns who generate the order, returns a vector of the columns in algebraic form. This is useful to save an order as being basis-independant (e.g. if you recompute the algebra A, you may change the chosen maximal order).*/
GEN
algordertoorderalg(GEN A, GEN O)
{
  pari_sp av = avma;
  long lO = lg(O), i;
  GEN v = cgetg(lO, t_VEC);
  for (i = 1; i < lO; i++) gel(v, i) = lift(algbasistoalg(A, gel(O, i)));
  return gerepilecopy(av, v);
}

/*Given an order O in A, returns the discriminant of the order, which is disc(algebra)*level*/
static GEN
algorderdisc(GEN A, GEN O, int reduced, int factored)
{
  pari_sp av = avma;
  GEN F = alg_get_center(A);/*The centre*/
  GEN zdisc = sqri(mulii(QM_det(O), algdiscnorm(A)));/*The norm to Q of the (non-reduced) discriminant*/
  GEN I = gen_0;/*Stores the ideal generated by d(a1, a2, a3, a4) as we run across the basis.*/
  long lbas=lg(O), i1, i2, i3, i4;
  for (i1 = 1; i1 < lbas; i1++) {
    for (i2 = i1 + 1; i2 < lbas; i2++) {
      for (i3 = i2 + 1; i3 < lbas; i3++) {
        for (i4 = i3 + 1; i4 < lbas; i4++) {
          GEN d = algd(A, mkvec4(gel(O, i1), gel(O, i2), gel(O, i3), gel(O, i4)));
          I = idealadd(F, I, d);
          if (equalii(idealnorm(F, I), zdisc)) {/*We are done!*/
            if (!reduced && !factored) return gerepilecopy(av, I);
            GEN fact = idealfactor(F, I);
            if (!reduced && factored) return gerepileupto(av, fact);/*From now on, we want it reduced*/
            long lp = lg(gel(fact, 1)), i;
            for (i = 1; i < lp; i++) gcoeff(fact, i, 2) = shifti(gcoeff(fact, i, 2), -1);/*Divide by 2*/
            if (factored) return gerepilecopy(av, fact);
            return gerepileupto(av, idealfactorback(F, fact, NULL, 0));
          }
        }
      }
    }
  }
  pari_warn(warner, "We did not succeed in finding the discriminant.");
  return gc_const(av, gen_0);
}

/*Returns the level of the order O.*/
GEN
algorderlevel(GEN A, GEN O, int factored)
{
  pari_sp av = avma;
  GEN F = alg_get_center(A);
  if (gequal1(O)) {/*Diagonal identity, so level 1.*/
    if (factored) return idealfactor(F, gen_1);
    return gen_1;
  }
  GEN Odisc = algorderdisc(A, O, 0, 0);
  GEN Adisc = idealpows(F, algreduceddisc(A), 2);
  GEN lsqr = idealdivexact(F, Odisc, Adisc);/*Divides exactly*/
  GEN fact = idealfactor(F, lsqr);
  long lp = lg(gel(fact, 1)), i;
  for (i = 1; i < lp; i++) gcoeff(fact, i, 2) = shifti(gcoeff(fact, i, 2), -1);/*Divide by 2*/
  if(factored) return gerepilecopy(av, fact);
  return gerepileupto(av, idealfactorback(F, fact, NULL, 0));
}



/*SECTION 4: FINCKE POHST FOR FLAG=2 WITH PRUNING*/


/*4: SUPPORTING METHODS TO FINCKE POHST*/

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

static GEN
clonefill(GEN S, long s, long t)
{ /* initialize to dummy values */
  GEN T = S, dummy = cgetg(1, t_STR);
  long i;
  for (i = s+1; i <= t; i++) gel(S,i) = dummy;
  S = gclone(S); if (isclone(T)) gunclone(T);
  return S;
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


/*4: MAIN FINCKE POHST METHODS*/

/*prunetype=0 means no pruning, 1 means linear pruning ala Schnorr-Euchner*/
GEN
qfminim_prune(GEN M, GEN C, int prunetype, long prec)
{
  GEN res = fincke_pohst_prune(M, C, prunetype, prec);
  if (!res) pari_err_PREC("qfminim");
  return res;
}

/*Solve q(x)=x~*M*x <= C, M is positive definite with real entries. We use pruning if prunetype=1. This is mostly a copy of fincke_pohst from bibli1.c.*/
GEN
fincke_pohst_prune(GEN M, GEN C, int prunetype, long PREC)
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
  GEN prune;
  if (!prunetype) prune = const_vec(lM - 1, gen_1);/*No funny business, just normal Fincke-Pohst.*/
  else {/*Linear pruning*/
    /*GEN con = dbltor(1.05);*/
    prune = cgetg(lM, t_VEC);
    long n = lM - 1;
    /*for (i = 1; i < lM; i++) gel(prune, i) = gmin_shallow(gen_1, divrs(mulrs(con, lM - i), n));*/
    for (i = 1; i < lM; i++) gel(prune, i) = rdivss(lM - i, n, prec);
  }
  GEN q = gaussred_from_QR(R, gprecision(Vnorm));
  if (q) res = smallvectors_prune(q, C, prune);
  if (!res) return gc_NULL(av);
  return gerepileupto(av, ZM_mul(U, res));
}

/*
The quadratic form is given by sum_{i=1}^n q_{i,i}(x_i+sum_{j=i+1}^n q_{i,j}*xj)^2
prune should be a vector of length n, entries real/integer in (0, 1], used to prune the search space.
*/
static GEN
smallvectors_prune(GEN q, GEN C, GEN prune)
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


/*SECTION 5: CHANGING ORDER METHODS, WHICH WERE DELETED FROM LIBPARI*/

/*Here because it was deleted from libpari*/
static GEN
my_alg_changeorder(GEN al, GEN ord)
{
  GEN al2, mt, iord, mtx;
  long i, n;
  pari_sp av = avma;

  if (!gequal0(gel(al, 10)))
    pari_err_DOMAIN("my_alg_changeorder","characteristic","!=", gen_0, gel(al, 10));
  n = alg_get_absdim(al);

  iord = QM_inv(ord);
  al2 = shallowcopy(al);

  gel(al2, 7) = RgM_mul(gel(al2, 7), ord);

  gel(al2, 8) = RgM_mul(iord, gel(al, 8));

  mt = cgetg(n + 1,t_VEC);
  gel(mt, 1) = matid(n);
  for (i = 2; i <= n; i++) {
    mtx = algbasismultable(al,gel(ord, i));
    gel(mt, i) = RgM_mul(iord, RgM_mul(mtx, ord));
  }
  gel(al2, 9) = mt;

  gel(al2, 11) = algtracebasis(al2);

  return gerepilecopy(av, al2);
}

static GEN
elementabsmultable(GEN mt, GEN x)
{
  GEN d, z = elementabsmultable_Z(mt, Q_remove_denom(x, &d));
  return (z && d)? ZM_Z_div(z, d): z;
}
static GEN
elementabsmultable_Fp(GEN mt, GEN x, GEN p)
{
  GEN z = elementabsmultable_Z(mt, x);
  return z? FpM_red(z, p): z;
}
static GEN
algbasismultable(GEN al, GEN x)
{
  pari_sp av = avma;
  GEN z, p = alg_get_char(al), mt = alg_get_multable(al);
  z = signe(p)? elementabsmultable_Fp(mt, x, p): elementabsmultable(mt, x);
  if (!z)
  {
    long D = lg(mt) - 1;
    set_avma(av); return zeromat(D, D);
  }
  return gerepileupto(av, z);
}

static GEN
algtracebasis(GEN al)
{
  pari_sp av = avma;
  GEN mt = alg_get_multable(al), p = alg_get_char(al);
  long i, l = lg(mt);
  GEN v = cgetg(l, t_VEC);
  if (signe(p)) for (i = 1; i < l; i++) gel(v, i) = FpM_trace(gel(mt, i), p);
  else          for (i = 1; i < l; i++) gel(v, i) = ZM_trace(gel(mt, i));
  return gerepileupto(av,v);
}

static GEN
elementabsmultable_Z(GEN mt, GEN x)
{
  long i, l = lg(x);
  GEN z = NULL;
  for (i = 1; i < l; i++)
  {
    GEN c = gel(x, i);
    if (signe(c))
    {
      GEN M = ZM_Z_mul(gel(mt, i), c);
      z = z? ZM_add(z, M): M;
    }
  }
  return z;
}

static GEN
FpM_trace(GEN x, GEN p)
{
  long i, lx = lg(x);
  GEN t;
  if (lx < 3) return lx == 1? gen_0: gcopy(gcoeff(x, 1, 1));
  t = gcoeff(x, 1, 1);
  for (i = 2; i < lx; i++) t = Fp_add(t, gcoeff(x, i, i), p);
  return t;
}

static GEN
ZM_trace(GEN x)
{
  long i, lx = lg(x);
  GEN t;
  if (lx < 3) return lx == 1? gen_0: gcopy(gcoeff(x, 1, 1));
  t = gcoeff(x, 1, 1);
  for (i = 2; i < lx; i++) t = addii(t, gcoeff(x, i, i));
  return t;
}



