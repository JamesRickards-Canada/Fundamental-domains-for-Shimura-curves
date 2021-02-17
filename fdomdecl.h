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
GEN divoo(GEN a, GEN b);

//LISTS
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
GEN mat_choleskydecomp(GEN A, int rcoefs, long prec);
GEN mat_nfcholesky(GEN nf, GEN A);

//SECTION 2: GEOMETRY METHODS


//DISTANCES
GEN hdist(GEN z1, GEN z2, long prec);
GEN hdist_ud(GEN z1, GEN z2, long prec);

//FUNDAMENTAL DOMAIN COMPUTATION
GEN isometriccircle_mats(GEN g, GEN mats, GEN *data, GEN (*gamtopsl)(GEN *, GEN, long), GEN tol, long prec);
GEN isometriccircle_psu(GEN g, GEN tol, long prec);
GEN normalizedbasis(GEN G, GEN U, GEN mats, GEN gamid, GEN *data, GEN (*gamtopsl)(GEN *, GEN, long), GEN (*eltmul)(GEN *, GEN, GEN), GEN (*eltinv)(GEN *, GEN), int (*istriv)(GEN *, GEN), GEN tol, long prec);
GEN normalizedboundary(GEN G, GEN mats, GEN id, GEN *data, GEN (*gamtopsl)(GEN *, GEN, long), GEN tol, long prec);
GEN psltopsu(GEN g, GEN p);
GEN psltopsu_mats(GEN g, GEN M);
GEN psltopsu_transmats(GEN p);
GEN randompoint_ud(GEN R, long prec);
GEN randompoint_udarc(GEN R, GEN ang1, GEN ang2, long prec);
GEN reduceelt_givennormbound(GEN U, GEN g, GEN z, GEN gamid, GEN *data, GEN (*gamtopsl)(GEN *, GEN, long), GEN (*eltmul)(GEN *, GEN, GEN), GEN tol, long prec);
GEN reduceelt_givenpsu(GEN G, GEN Gmats, GEN g, GEN gmat, GEN z, GEN gamid, GEN *data, GEN (*eltmul)(GEN *, GEN, GEN), GEN tol, long prec);
GEN reducepoint(GEN U, GEN z, GEN gamid, GEN *data, GEN (*eltmul)(GEN *, GEN, GEN), GEN tol, long prec);
GEN rootgeodesic_fd(GEN U, GEN g, GEN gamid, GEN *data, GEN (*gamtopsl)(GEN *, GEN, long), GEN (*eltmul)(GEN *, GEN, GEN), GEN (*eltinv)(GEN *, GEN), GEN tol, long prec);
GEN rootgeodesic_ud(GEN M, GEN mats, GEN tol, long prec);
GEN rootgeodesic_uhp(GEN M, GEN tol, long prec);

//PRINTING TO PLOTVIEWER
void python_printarcs(GEN arcs, char *filename, int view, char *extrainput, long prec);
void python_plotviewer(char *input);
void python_printfdom(GEN U, char *filename, long prec);

//HELPER METHODS
GEN deftol(long prec);


//SECTION 3: QUATERNIONIC METHODS


//QUATERNION METHODS
GEN algfdomarea(GEN A, long prec);
GEN algfdom(GEN A, GEN p, int dispprogress, GEN area, GEN ANRdata, long prec);
GEN algfdomreduce(GEN A, GEN U, GEN g, GEN z, long prec);
GEN algfdomrootgeodesic(GEN A, GEN U, GEN g, long prec);
GEN algmulvec(GEN A, GEN G, GEN L);
GEN algramifiedplacesf(GEN A);
GEN algsmallnorm1elts(GEN A, GEN C, GEN p, GEN z, long prec);



//TEMPORARY
GEN algnormform(GEN A, long prec);
GEN algabsrednorm(GEN A, GEN p, GEN z, long prec);
GEN algfdom_test(GEN A, GEN p, int dispprogress, GEN area, GEN ANRdata, long prec);
