//STRUCTURES

typedef struct listtype1{//A circular list of GENs, stores data, next, and previous terms
  GEN data; 
  struct listtype1 *next;
  struct listtype1 *prev;
}clist;

typedef struct listtype2{//A generic linked list of GENs, stores data and next term
  GEN data; 
  struct listtype2 *next;
}glist;

typedef struct listtype3{//A generic linked list of longs, stores data and next term
  long data; 
  struct listtype3 *next;
}llist;


//METHODS


//SECTION 1: BASE METHODS

//INFINITY 
GEN addoo(GEN a, GEN b);
GEN divoo(GEN a, GEN b);

//LISTS
void clist_free(clist *l, long length);
void clist_putbefore(clist **head_ref, GEN new_data);
void clist_putafter(clist **head_ref, GEN new_data);
GEN clist_togvec(clist *l, long length, int dir);
void glist_free(glist *l);
GEN glist_pop(glist **head_ref);
void glist_putstart(glist **head_ref, GEN new_data);
GEN glist_togvec(glist *l, long length, int dir);
GEN glist_togvec_append(glist *l, GEN v, long length, int dir);
void llist_free(llist *l);
long llist_pop(llist **head_ref);
void llist_putstart(llist **head_ref, long new_data);
GEN llist_togvec(llist *l, long length, int dir);
GEN llist_tovecsmall(llist *l, long length, int dir);

//SHORT VECTORS IN LATTICES
GEN lat_smallvectors(GEN A, GEN C1, GEN C2, GEN condition, int onesign, int isintegral, int rdataonly, long prec);
GEN lat_smallvectors_givendata(GEN chol, GEN U, GEN perminv, GEN C1, GEN C2, GEN condition, int onesign, long prec);
GEN lat_smallvectors_tc(GEN A, GEN C1, GEN C2, int onesign, int isintegral, long prec);
GEN lat_smallvectors_cholesky(GEN Q, GEN C1, GEN C2, GEN condition, int onesign, long prec);
GEN mat_choleskydecomp(GEN A, int rcoefs, long prec);
GEN mat_choleskydecomp_tc(GEN A, int rcoefs, long prec);


//SECTION 2: GEOMETRY METHODS


//BASIC LINE, CIRCLE, AND POINT OPERATIONS
GEN arc_init(GEN c, GEN p1, GEN p2, int dir, long prec);
GEN arc_init_tc(GEN c, GEN p1, GEN p2, int dir, long prec);
GEN arc_midpoint(GEN c, GEN p1, GEN p2, GEN tol, long prec);
GEN arc_midpoint_tc(GEN c, GEN p1, GEN p2, long prec);
GEN circle_angle(GEN c1, GEN c2, GEN p, GEN tol, long prec);
GEN circle_angle_tc(GEN c1, GEN c2, GEN p, long prec);
GEN circle_fromcp(GEN cent, GEN p, long prec);
GEN circle_fromppp(GEN p1, GEN p2, GEN p3, GEN tol, long prec);
GEN circle_fromppp_tc(GEN p1, GEN p2, GEN p3, long prec);
GEN circle_tangentslope(GEN c, GEN p, long prec);
GEN circle_tangentslope_tc(GEN c, GEN p, long prec);
GEN crossratio(GEN a, GEN b, GEN c, GEN d);
GEN line_angle(GEN l1, GEN l2, long prec);
GEN line_fromsp(GEN s, GEN p);
GEN line_frompp(GEN p1, GEN p2);
GEN mat_eval(GEN M, GEN x);
GEN mat_eval_tc(GEN M, GEN x);
GEN midpoint(GEN p1, GEN p2);
GEN mobius(GEN M, GEN c, GEN tol, long prec);
GEN mobius_tc(GEN M, GEN c, long prec);
GEN perpbis(GEN p1, GEN p2);
GEN radialangle(GEN c, GEN p, GEN tol, long prec);
GEN radialangle_tc(GEN c, GEN p, long prec);
GEN slope(GEN p1, GEN p2);

//INTERSECTION OF LINES/CIRCLES
GEN arc_int(GEN c1, GEN c2, GEN tol, long prec);
GEN arc_int_tc(GEN c1, GEN c2, long prec);
GEN arcseg_int(GEN c, GEN l, GEN tol, long prec);
GEN arcseg_int_tc(GEN c, GEN l, long prec);
GEN circle_int(GEN c1, GEN c2, GEN tol, long prec);
GEN circle_int_tc(GEN c1, GEN c2, long prec);
GEN circleline_int(GEN c, GEN l, GEN tol, long prec);
GEN circleline_int_tc(GEN c, GEN l, long prec);
GEN genseg_int(GEN s1, GEN s2, GEN tol, long prec);
GEN genseg_int_tc(GEN s1, GEN s2, long prec);
GEN line_int(GEN l1, GEN l2, GEN tol, long prec);
GEN line_int_tc(GEN l1, GEN l2, long prec);
int onarc(GEN c, GEN p, GEN tol, long prec);
int onarc_tc(GEN c, GEN p, long prec);
int onseg(GEN l, GEN p, GEN tol, long prec);
int onseg_tc(GEN l, GEN p, long prec);
GEN seg_int(GEN l1, GEN l2, GEN tol, long prec);
GEN seg_int_tc(GEN l1, GEN l2, long prec);

//DISTANCES
GEN hdist(GEN z1, GEN z2, long prec);
GEN hdist_tc(GEN z1, GEN z2, long prec);
GEN hdist_ud(GEN z1, GEN z2, long prec);
GEN hpolygon_area(GEN circles, GEN vertices, GEN tol, long prec);
GEN hpolygon_area_tc(GEN circles, GEN vertices, long prec);

//FUNDAMENTAL DOMAIN COMPUTATION
GEN edgepairing(GEN U, GEN tol, int rboth, long prec);
GEN edgepairing_tc(GEN U, long prec);
GEN isometriccircle_mats(GEN g, GEN mats, GEN *data, GEN (*gamtopsl)(GEN *, GEN, long), GEN tol, long prec);
GEN isometriccircle_psl(GEN g, GEN p, GEN tol, long prec);
GEN isometriccircle_psl_mats(GEN g, GEN mats, GEN tol, long prec);
GEN isometriccircle_psu(GEN g, GEN tol, long prec);
GEN normalizedbasis(GEN G, GEN U, GEN mats, GEN gamid, GEN *data, GEN (*gamtopsl)(GEN *, GEN, long), GEN (*eltmul)(GEN *, GEN, GEN), GEN (*eltinv)(GEN *, GEN), int (*istriv)(GEN *, GEN), GEN tol, long prec);
GEN normalizedboundary_append(GEN Ubase, GEN G, GEN mats, GEN id, GEN tol, long prec);
GEN normalizedboundary_givenU(GEN Ubase, GEN G, GEN mats, GEN id, GEN *data, GEN (*gamtopsl)(GEN *, GEN, long), GEN tol, long prec);
GEN normalizedboundary_givencircles(GEN G, GEN mats, GEN id, GEN tol, long prec);
GEN normalizedboundary(GEN G, GEN mats, GEN id, GEN *data, GEN (*gamtopsl)(GEN *, GEN, long), GEN tol, long prec);
long normalizedboundary_outside(GEN U, GEN z, GEN tol, long prec);
long normalizedboundary_outside_tc(GEN U, GEN z, long prec);
GEN normalizedboundary_sideint(GEN U, GEN c, int start, GEN tol, long prec);
GEN psltopsu(GEN g, GEN p);
GEN psltopsu_mats(GEN g, GEN M);
GEN psltopsu_transmats(GEN p);
GEN psl_roots(GEN M, GEN tol, long prec);
GEN randompoint_ud(GEN R, long prec);
GEN reduceelt_givennormbound(GEN U, GEN g, GEN z, GEN gamid, GEN *data, GEN (*gamtopsl)(GEN *, GEN, long), GEN (*eltmul)(GEN *, GEN, GEN), GEN tol, long prec);
GEN reduceelt_givenpsu(GEN G, GEN Gmats, GEN g, GEN gmat, GEN z, GEN gamid, GEN *data, GEN (*eltmul)(GEN *, GEN, GEN), GEN tol, long prec);
GEN reducepoint(GEN U, GEN z, GEN gamid, GEN *data, GEN (*eltmul)(GEN *, GEN, GEN), GEN tol, long prec);
GEN rootgeodesic_fd(GEN U, GEN g, GEN gamid, GEN *data, GEN (*gamtopsl)(GEN *, GEN, long), GEN (*eltmul)(GEN *, GEN, GEN), GEN (*eltinv)(GEN *, GEN), GEN tol, long prec);
GEN rootgeodesic_ud(GEN M, GEN mats, GEN tol, long prec);
GEN rootgeodesic_uhp(GEN M, GEN tol, long prec);
GEN rootgeodesic_uhp_tc(GEN M, long prec);

//PRINTING TO PLOTVIEWER
void python_printarcs(GEN arcs, char *filename, int view, char *extrainput, long prec);
void python_plotviewer(char *input);
void python_printfdom(GEN U, char *filename, long prec);

//HELPER METHODS
GEN anglediff(GEN ang, GEN bot, GEN tol, long prec);
GEN atanoo(GEN x, long prec);
GEN deftol(long prec);
int gcmp_strict(void *data, GEN x, GEN y);
int geom_check(GEN c);
GEN shiftangle(GEN ang, GEN bot, GEN tol, long prec);
long tolcmp(GEN x, GEN y, GEN tol, long prec);
int tolcmp_sort(void *data, GEN x, GEN y);
int toleq(GEN x, GEN y, GEN tol, long prec);


//SECTION 3: QUATERNIONIC METHODS


//QUATERNION METHODS
GEN algfdarea(GEN A, long prec);
GEN algfd(GEN A, GEN p, int dispprogress, GEN ANRdata, GEN area, long prec);
GEN algramifiedplacesf(GEN A);
GEN algsmallnorm1elts(GEN A, GEN C, GEN p, GEN z, long prec);