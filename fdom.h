/*FDOM METHODS*/


/*SECTION 1: BASE METHODS*/

/*LISTS*/
GEN veclist_append(GEN v, long *vind, long *vlen, GEN x);
GEN vecsmalllist_append(GEN v, long *vind, long *vlen, long x);

/*SHORT VECTORS IN LATTICES*/


/*SECTION 2: GEOMETRY METHODS*/

/*SHALLOW RETRIEVAL METHODS*/
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

/*BASIC LINE, CIRCLE, AND POINT OPERATIONS*/

/*DISTANCES/AREAS*/

/*FUNDAMENTAL DOMAIN COMPUTATION*/

/*FUNDAMENTAL DOMAIN OTHER COMPUTATIONS*/

/*HELPER METHODS*/
GEN deftol(long prec);
