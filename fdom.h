/*fdom.c methods*/

/*INLINE SHALLOW RETRIEVAL METHODS*/

/*1: GEOMETRIC DATA*/
INLINE GEN
gdat_get_tol(GEN gd){return gel(gd, 1);}
INLINE GEN
gdat_get_p(GEN gd){return gel(gd, 2);}
INLINE GEN
gdat_get_pscale(GEN gd){return gel(gd, 3);}

/*2: NORMALIZED BOUNDARY*/
INLINE GEN
normbound_get_elts(GEN U){return gel(U, 1);}
INLINE GEN
normbound_get_sides(GEN U){return gel(U, 2);}
INLINE GEN
normbound_get_vcors(GEN U){return gel(U, 3);}
INLINE GEN
normbound_get_vargs(GEN U){return gel(U, 4);}
INLINE GEN
normbound_get_kact(GEN U){return gel(U, 5);}
INLINE GEN
normbound_get_area(GEN U){return gel(U, 6);}
INLINE GEN
normbound_get_spair(GEN U){return gel(U, 7);}
INLINE GEN
normbound_get_gdat(GEN U){return gel(U, 8);}


/*SECTION 1: GEOMETRIC METHODS*/

/*1: LINES AND ARCS*/
/*1: MATRIX ACTION ON GEOMETRY*/
GEN klein_act(GEN M, GEN z);
GEN pgl_act(GEN M, GEN z);

/*1: TRANSFER BETWEEN MODELS*/
/*1: DISTANCES/AREAS*/
GEN hdiscradius(GEN area, long prec);
GEN hdiscrandom(GEN R, long prec);
GEN hdiscrandom_arc(GEN R, GEN ang1, GEN ang2, long prec);
GEN hdist(GEN z1, GEN z2, long flag, long prec);

/*1: OPERATIONS ON COMPLEX REALS*/
/*1: TOLERANCE*/
GEN deftol(long prec);


/*SECTION 2: FUNDAMENTAL DOMAIN GEOMETRY*/

/*2: NORMALIZED BOUNDARY*/

/*SECTION 3: QUATERNION ALGEBRA METHODS*/

/*3: INITIALIZE SYMMETRIC SPACE*/
GEN algsymminit(GEN A, GEN O, GEN type, GEN p, long prec);

/*3: ALGEBRA METHODS FOR GEOMETRY*/

/*3: ALGEBRA HELPER METHODS*/
GEN algalgto1ijk(GEN A, GEN x);
GEN algbasisto1ijk(GEN A, GEN x);
GEN algramifiedplacesf(GEN A);


