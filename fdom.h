/*fdom.c methods*/

enum {afuch_A = 1, afuch_ONORMREAL, afuch_KLEINMATS, afuch_QFMATS, afuch_GDAT, afuch_FDOMDAT, afuch_FDOM, afuch_SIG, afuch_PRES, afuch_SAVEDELTS, afuch_NORMALIZERNORMS};

/*INLINE SHALLOW RETRIEVAL METHODS*/

/*1: GEOMETRIC DATA*/
INLINE GEN
gdat_get_tol(GEN gd) { return gel(gd, 1); }
INLINE GEN
gdat_get_p(GEN gd) { return gel(gd, 2); }

/*2: NORMALIZED BOUNDARY*/
INLINE GEN
normbound_get_elts(GEN U) { return gel(U, 1); }
INLINE GEN
normbound_get_sides(GEN U) { return gel(U, 2); }
INLINE GEN
normbound_get_vcors(GEN U) { return gel(U, 3); }
INLINE GEN
normbound_get_vargs(GEN U) { return gel(U, 4); }
INLINE long
normbound_get_cross(GEN U) { return itos(gel(U, 5)); }
INLINE GEN
normbound_get_kact(GEN U) { return gel(U, 6); }
INLINE GEN
normbound_get_area(GEN U) { return gel(U, 7); }
INLINE GEN
normbound_get_spair(GEN U) { return gel(U, 8); }
INLINE GEN
normbound_get_infinite(GEN U) { return gel(U, 9); }

/*3: ARITHMETIC FUCHSIAN GROUPS*/
INLINE GEN
afuch_get_O(GEN X) { return gel(X, 1); }
INLINE GEN
afuch_get_Oinv(GEN X) { return gel(X, 2); }
INLINE GEN
afuch_get_Oconj(GEN X) { return gel(X, 3); }
INLINE GEN
afuch_get_Omultable(GEN X) { return gel(X, 4); }
INLINE GEN
afuch_get_Onormdat(GEN X) { return gel(X, 5); }
INLINE GEN
afuch_get_type(GEN X) { return gel(X, 6); }
INLINE GEN
afuch_get_alg(GEN X) { return obj_check(X, afuch_A); }
INLINE GEN
afuch_get_Onormreal(GEN X) { return obj_check(X, afuch_ONORMREAL); }
INLINE GEN
afuch_get_kleinmats(GEN X) { return obj_check(X, afuch_KLEINMATS); }
INLINE GEN
afuch_get_qfmats(GEN X) { return obj_check(X, afuch_QFMATS); }
INLINE GEN
afuch_get_gdat(GEN X) { return obj_check(X, afuch_GDAT); }
INLINE GEN
afuch_get_O1area(GEN X) { GEN fdomdat = obj_check(X, afuch_FDOMDAT); return fdomdat? gel(fdomdat, 1) : NULL; }
INLINE GEN
afuch_get_bestC(GEN X) { GEN fdomdat = obj_check(X, afuch_FDOMDAT); return fdomdat? gel(fdomdat, 2) : NULL; }
INLINE GEN
afuch_get_R(GEN X) { GEN fdomdat = obj_check(X, afuch_FDOMDAT); return fdomdat? gel(fdomdat, 3) : NULL; }
INLINE GEN
afuch_get_epsilon(GEN X) { GEN fdomdat = obj_check(X, afuch_FDOMDAT); return fdomdat? gel(fdomdat, 4) : NULL; }
INLINE GEN
afuch_get_passes(GEN X) { GEN fdomdat = obj_check(X, afuch_FDOMDAT); return fdomdat? gel(fdomdat, 5) : NULL; }
INLINE GEN
afuch_get_fdom(GEN X) { return obj_check(X, afuch_FDOM); }
INLINE GEN
afuch_get_sig(GEN X) { return obj_check(X, afuch_SIG); }
INLINE GEN
afuch_get_pres(GEN X) { return obj_check(X, afuch_PRES); }
INLINE GEN
afuch_get_savedelts(GEN X) { return obj_check(X, afuch_SAVEDELTS); }
INLINE GEN
afuch_get_normalizernorms(GEN X) { return obj_check(X, afuch_NORMALIZERNORMS); }
INLINE long
afuch_get_prec(GEN X) { return lg(gdat_get_tol(afuch_get_gdat(X))); }

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
GEN afuch_newtype(GEN X, GEN type);
void afuch_changep(GEN X, GEN p);

/*3: ALGEBRA FUNDAMENTAL DOMAIN METHODS*/
GEN afucharea(GEN X);
GEN afuchelts(GEN X);
void afuchfdom(GEN X);
GEN afuchgeodesic(GEN X, GEN g);
GEN afuchlist(GEN F, GEN Amin, GEN Amax, long split);
GEN afuchnormalizernorms(GEN X);
GEN afuchpresentation(GEN X);
GEN afuchsignature(GEN X);
GEN afuchspair(GEN X);
GEN afuchword(GEN X, GEN g);

/*3: NON NORM 1 METHODS*/
GEN bnf_make_unitnorms(GEN B, long split, long prec);

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

/*SECTION 2: TUNING*/
GEN tune_Cn(long n, GEN Cmin, GEN Cmax, long testsperalg, long tests, long prec);
