/*fdom.c methods*/

/*SECTION 0: MISCELLANEOUS METHODS*/

/*0: LISTS*/
GEN veclist_append(GEN v, long *vind, long *vlen, GEN x);
GEN vecsmalllist_append(GEN v, long *vind, long *vlen, long x);

/*0: SHORT VECTORS IN LATTICES*/


/*SECTION 1: GEOMETRIC METHODS*/

/*1: SEGMENT RETRIEVAL METHODS*/
INLINE GEN
seg_get_a(GEN L){return gel(L, 1);}
INLINE GEN
seg_get_b(GEN L){return gel(L, 2);}
INLINE GEN
seg_get_c(GEN L){return gel(L, 3);}
INLINE GEN
seg_get_start(GEN L){return gel(L, 4);}
INLINE GEN
seg_get_end(GEN L){return gel(L, 5);}

/*1: MATRIX ACTION ON GEOMETRY*/
GEN klein_act(GEN M, GEN z);
GEN pgl_act(GEN M, GEN z);

/*1: DISTANCES/AREAS*/
GEN hdiscradius(GEN area, long prec);
GEN hdiscrandom(GEN R, long prec);
GEN hdiscrandom_arc(GEN R, GEN ang1, GEN ang2, long prec);
GEN hdist(GEN z1, GEN z2, long flag, long prec);

/*1: TOLERANCE*/
GEN deftol(long prec);


/*SECTION 2: FUNDAMENTAL DOMAIN GEOMETRY*/

/*2: GEOMDATA RETRIEVAL METHODS*/
INLINE long
gdat_get_prec(GEN gd){return gel(gd, 1)[2];}
INLINE GEN
gdat_get_tol(GEN gd){return gel(gd, 2);}
INLINE GEN
gdat_get_p(GEN gd){return gel(gd, 3);}
INLINE GEN
gdat_get_pc(GEN gd){return gel(gd, 4);}
INLINE GEN
gdat_get_pscale(GEN gd){return gel(gd, 5);}


