/*fdom.c methods*/

enum {afuch_A = 1, afuch_ONORMREAL, afuch_KLEINMATS, afuch_QFMATS, afuch_GDAT, afuch_FDOMDAT, afuch_FDOM, afuch_PRES};

/*INLINE SHALLOW RETRIEVAL METHODS*/

/*1: GEOMETRIC DATA*/
INLINE GEN
gdat_get_tol(GEN gd){return gel(gd, 1);}
INLINE GEN
gdat_get_lowtol(GEN gd){return gel(gd, 2);}
INLINE GEN
gdat_get_p(GEN gd){return gel(gd, 3);}

/*2: NORMALIZED BOUNDARY*/
INLINE GEN
normbound_get_elts(GEN U){return gel(U, 1);}
INLINE GEN
normbound_get_sides(GEN U){return gel(U, 2);}
INLINE GEN
normbound_get_vcors(GEN U){return gel(U, 3);}
INLINE GEN
normbound_get_vargs(GEN U){return gel(U, 4);}
INLINE long
normbound_get_cross(GEN U){return itos(gel(U, 5));}
INLINE GEN
normbound_get_kact(GEN U){return gel(U, 6);}
INLINE GEN
normbound_get_area(GEN U){return gel(U, 7);}
INLINE GEN
normbound_get_spair(GEN U){return gel(U, 8);}
INLINE GEN
normbound_get_infinite(GEN U){return gel(U, 9);}

/*3: ARITHMETIC FUCHSIAN GROUPS*/
INLINE GEN
afuch_get_O(GEN X){return gel(X, 1);}
INLINE GEN
afuch_get_Oinv(GEN X){return gel(X, 2);}
INLINE GEN
afuch_get_Oconj(GEN X){return gel(X, 3);}
INLINE GEN
afuch_get_Omultable(GEN X){return gel(X, 4);}
INLINE GEN
afuch_get_Onormdat(GEN X){return gel(X, 5);}
INLINE GEN
afuch_get_type(GEN X){return gel(X, 6);}
INLINE GEN
afuch_get_alg(GEN X){return gmael(X, 7, afuch_A);}
INLINE GEN
afuch_get_Onormreal(GEN X){return gmael(X, 7, afuch_ONORMREAL);}
INLINE GEN
afuch_get_kleinmats(GEN X){return gmael(X, 7, afuch_KLEINMATS);}
INLINE GEN
afuch_get_qfmats(GEN X){return gmael(X, 7, afuch_QFMATS);}
INLINE GEN
afuch_get_gdat(GEN X){return gmael(X, 7, afuch_GDAT);}
INLINE GEN
afuch_get_area(GEN X){return gmael3(X, 7, afuch_FDOMDAT, 1);}
INLINE GEN
afuch_get_bestC(GEN X){return gmael3(X, 7, afuch_FDOMDAT, 2);}
INLINE GEN
afuch_get_R(GEN X){return gmael3(X, 7, afuch_FDOMDAT, 3);}
INLINE GEN
afuch_get_epsilon(GEN X){return gmael3(X, 7, afuch_FDOMDAT, 4);}
INLINE GEN
afuch_get_passes(GEN X){return gmael3(X, 7, afuch_FDOMDAT, 5);}
INLINE GEN
afuch_get_fdom(GEN X){return gmael(X, 7, afuch_FDOM);}
INLINE GEN
afuch_get_pres(GEN X){return gmael(X, 7, afuch_PRES);}

/*SECTION 1: GEOMETRIC METHODS*/

/*1: LINES AND ARCS*/
/*1: MATRIX ACTION ON GEOMETRY*/
GEN klein_act_i(GEN M, GEN z);
GEN klein_act(GEN M, GEN z, long prec);
GEN pgl_act(GEN M, GEN z);

/*1: TRANSFER BETWEEN MODELS*/
/*1: DISTANCES/AREAS*/
/*1: OPERATIONS ON COMPLEX REALS*/
/*1: TOLERANCE*/
GEN deftol(long prec);
GEN deflowtol(long prec);


/*SECTION 2: FUNDAMENTAL DOMAIN GEOMETRY*/

/*2: NORMALIZED BOUNDARY*/
/*2: NORMALIZED BOUNDARY APPENDING*/
/*2: NORMALIZED BOUNDARY ANGLES*/
/*2: NORMALIZED BASIS*/
/*2: NORMALIZED BOUNDARY REDUCTION*/
/*2: PRESENTATION*/

/*SECTION 3: QUATERNION ALGEBRA METHODS*/

/*3: INITIALIZE SYMMETRIC SPACE*/
GEN afuchinit(GEN A, GEN O, GEN type, GEN p, long prec);

/*3: ALGEBRA FUNDAMENTAL DOMAIN METHODS*/
GEN afuchfdom(GEN X);
GEN afuchicirc(GEN X, GEN g);
GEN afuchklein(GEN X, GEN g);
GEN afuchnormbasis(GEN X, GEN G);
GEN afuchnormbound(GEN X, GEN G);
GEN afuchnormbound_append(GEN X, GEN U, GEN G);
GEN afuchredelt(GEN X, GEN U, GEN g, GEN z);

/*3: ALGEBRA BASIC AUXILLARY METHODS*/

/*3: FINDING ELEMENTS*/
GEN afuchfindelts(GEN X, GEN z, GEN C, long maxelts);
GEN afuchqf(GEN X, GEN z);

/*3: ALGEBRA HELPER METHODS*/
GEN algalgto1ijk(GEN A, GEN x);
GEN algbasisto1ijk(GEN A, GEN x);
GEN algmulvec(GEN A, GEN G, GEN L);

/*3: ALGEBRA ORDER METHODS*/

