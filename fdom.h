/*fdom.c methods*/

/*CONSTANTS*/

/*Normalized boundary*/
enum { normbound_ELTS = 1, normbound_SIDES, normbound_VCORS, normbound_VARGS, normbound_CROSS, normbound_KACT, normbound_AREA, normbound_SPAIR, normbound_INFINITE };

/*Sets the ordering of the fundamental domain vector.*/
enum { afuch_A = 1, afuch_O, afuch_OINV, afuch_OCONJ, afuch_OMULTABLE, afuch_ONORMDAT, afuch_TYPE, afuch_ONORMREAL, afuch_KLEINMATS, afuch_QFMATS, afuch_GDAT, afuch_FDOMDAT, afuch_FDOM, afuch_SIGNATURE, afuch_PRES, afuch_ELLIPTIC, afuch_SAVEDELTS, afuch_NORMALIZERNORMS };

/*Constants for fundamental domain computation*/
enum { fdomdat_O1AREA = 1, fdomdat_BESTC, fdomdat_R, fdomdat_EPSILON, fdomdat_PASSES };


/*INLINE SHALLOW RETRIEVAL METHODS*/

/*1: GEOMETRIC DATA*/
INLINE GEN
gdat_get_tol(GEN gd) { return gel(gd, 1); }
INLINE GEN
gdat_get_p(GEN gd) { return gel(gd, 2); }

/*2: NORMALIZED BOUNDARY*/
INLINE GEN
normbound_get_elts(GEN U) { return gel(U, normbound_ELTS); }
INLINE GEN
normbound_get_sides(GEN U) { return gel(U, normbound_SIDES); }
INLINE GEN
normbound_get_vcors(GEN U) { return gel(U, normbound_VCORS); }
INLINE GEN
normbound_get_vargs(GEN U) { return gel(U, normbound_VARGS); }
INLINE long
normbound_get_cross(GEN U) { return itos(gel(U, normbound_CROSS)); }
INLINE GEN
normbound_get_kact(GEN U) { return gel(U, normbound_KACT); }
INLINE GEN
normbound_get_area(GEN U) { return gel(U, normbound_AREA); }
INLINE GEN
normbound_get_spair(GEN U) { return gel(U, normbound_SPAIR); }
INLINE GEN
normbound_get_infinite(GEN U) { return gel(U, normbound_INFINITE); }

/*3: ARITHMETIC FUCHSIAN GROUPS*/
INLINE GEN
afuch_get_alg(GEN X) { return gel(X, afuch_A); }
INLINE GEN
afuch_get_O(GEN X) { return gel(X, afuch_O); }
INLINE GEN
afuch_get_Oinv(GEN X) { return gel(X, afuch_OINV); }
INLINE GEN
afuch_get_Oconj(GEN X) { return gel(X, afuch_OCONJ); }
INLINE GEN
afuch_get_Omultable(GEN X) { return gel(X, afuch_OMULTABLE); }
INLINE GEN
afuch_get_Onormdat(GEN X) { return gel(X, afuch_ONORMDAT); }
INLINE GEN
afuch_get_type(GEN X) { return gel(X, afuch_TYPE); }
INLINE GEN
afuch_get_Onormreal(GEN X) { return gel(X, afuch_ONORMREAL); }
INLINE GEN
afuch_get_kleinmats(GEN X) { return gel(X, afuch_KLEINMATS); }
INLINE GEN
afuch_get_qfmats(GEN X) { return gel(X, afuch_QFMATS); }
INLINE GEN
afuch_get_gdat(GEN X) { return gel(X, afuch_GDAT); }
INLINE GEN
afuch_get_O1area(GEN X) { return gmael(X, afuch_FDOMDAT, fdomdat_O1AREA); }
INLINE GEN
afuch_get_bestC(GEN X) { return gmael(X, afuch_FDOMDAT, fdomdat_BESTC); }
INLINE GEN
afuch_get_R(GEN X) { return gmael(X, afuch_FDOMDAT, fdomdat_R); }
INLINE GEN
afuch_get_epsilon(GEN X) { return gmael(X, afuch_FDOMDAT, fdomdat_EPSILON); }
INLINE GEN
afuch_get_passes(GEN X) { return gmael(X, afuch_FDOMDAT, fdomdat_PASSES); }
INLINE GEN
afuch_get_fdom(GEN X) { return gel(X, afuch_FDOM); }
INLINE GEN
afuch_get_signature(GEN X) { return gel(X, afuch_SIGNATURE); }
INLINE GEN
afuch_get_presentation(GEN X) { return gel(X, afuch_PRES); }
INLINE GEN
afuch_get_elliptic(GEN X) { return gel(X, afuch_ELLIPTIC); }
INLINE GEN
afuch_get_savedelts(GEN X) { return gel(X, afuch_SAVEDELTS); }
INLINE GEN
afuch_get_normalizernorms(GEN X) { return gel(X, afuch_NORMALIZERNORMS); }
INLINE long
afuch_get_prec(GEN X) { return realprec(gdat_get_tol(afuch_get_gdat(X))); }

/*SECTION 1: GEOMETRIC METHODS*/

/*1: MATRIX ACTION ON GEOMETRY*/
GEN klein_act_i(GEN M, GEN z);
GEN klein_act(GEN M, GEN z, long prec);
GEN pgl_act(GEN M, GEN z);

/*1: TRANSFER BETWEEN MODELS*/
GEN disc_to_klein(GEN z);
GEN disc_to_plane(GEN z, GEN p);
GEN klein_to_disc(GEN z, GEN tol);
GEN klein_to_plane(GEN z, GEN p, GEN tol);
GEN plane_to_disc(GEN z, GEN p);
GEN plane_to_klein(GEN z, GEN p);

/*1: DISTANCES/AREAS*/
GEN hdiscrandom(GEN R, long prec);

/*SECTION 2: FUNDAMENTAL DOMAIN GEOMETRY*/

/*SECTION 3: QUATERNION ALGEBRA METHODS*/

/*3: INITIALIZE SYMMETRIC SPACE*/
GEN afuchinit(GEN A, GEN O, GEN type, int flag, long prec);
GEN afuchmoreprec(GEN X, long inc);
GEN afuchnewp(GEN X, GEN p);
GEN afuchnewtype(GEN X, GEN type);

/*3: ALGEBRA FUNDAMENTAL DOMAIN METHODS*/
GEN afuchalg(GEN X);
GEN afucharea(GEN X);
GEN afuchelliptic(GEN X);
GEN afuchelts(GEN X);
int afuchelttype(GEN X, GEN g);
GEN afuchmakefdom(GEN X);
GEN afuchgeodesic(GEN X, GEN g);
GEN afuchlist(GEN F, GEN bounds, long split);
GEN afuchminimalcycles(GEN X);
GEN afuchnormalizernorms(GEN X);
GEN afuchorder(GEN X);
GEN afuchpresentation(GEN X);
GEN afuchsides(GEN X);
GEN afuchsignature(GEN X);
GEN afuchspair(GEN X);
GEN afuchword(GEN X, GEN g);
GEN afuchvertices(GEN X, int model);

/*3: NON NORM 1 METHODS*/
GEN bnf_make_unitnorms(GEN B, long split, long prec);
GEN normalizer_make_norms(GEN B, long split, GEN ideals, long prec);

/*3: ALGEBRA BASIC AUXILLARY METHODS*/
GEN afuchconj(GEN X, GEN g);
int afuchistriv(GEN X, GEN g);
GEN afuchmul(GEN X, GEN g1, GEN g2);

/*3: FINDING ELEMENTS*/
GEN afuchfindoneelt(GEN X, GEN nm, GEN C);

/*3: ALGEBRA HELPER METHODS*/
GEN algab(GEN A);
GEN algalgto1ijk(GEN A, GEN x);
GEN algbasisto1ijk(GEN A, GEN x);
GEN alg1ijktoalg(GEN A, GEN x);
GEN alg1ijktobasis(GEN A, GEN x);
GEN algdiscnorm(GEN A);
GEN algmulvec(GEN A, GEN G, GEN L);
GEN algramifiedplacesf(GEN A);
GEN algreduceddisc(GEN A);
long algsplitoo(GEN A);

/*3: ALGEBRA ORDER METHODS*/
int algisorder(GEN A, GEN O);
GEN algorderalgtoorder(GEN A, GEN Oalg);
GEN algordertoorderalg(GEN A, GEN O);
GEN algorderlevel(GEN A, GEN O, int factored);

/*SECTION 4: FINCKE POHST FOR FLAG=2 WITH PRUNING*/

/*4: MAIN FINCKE POHST METHODS*/
GEN qfminim_prune(GEN M, GEN C, int prunetype, long prec);
GEN fincke_pohst_prune(GEN M, GEN C, int prunetype, long PREC);


/*fdom_extra.c methods*/

/*SECTION 1: VISUALIZATION*/
void afuchfdom_latex(GEN X, char *filename, int model, int boundcircle, int compile, int open);
void afuchfdom_python(GEN X, char *filename);
void afuchgeodesic_python(GEN X, GEN g, char *filename);
void fdomviewer(char *input);

/*SECTION 2: EICHLER ORDERS*/
GEN algeichlerorder(GEN A, GEN I);

/*SECTION 3: TESTING AND TUNING*/
long afuchcheck(GEN X);
GEN tune_Cn(long n, GEN Cmin, GEN Cmax, long testsperalg, long tests, long prec);
