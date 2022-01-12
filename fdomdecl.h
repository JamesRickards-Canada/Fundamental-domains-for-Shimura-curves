//STRUCTURES

typedef struct listtype1{//A generic linked list of GENs, stores data and next term
  GEN data; 
  struct listtype1 *next;
}glist;

typedef struct listtype2{//A generic linked list of longs, stores data and next term
  long data; 
  struct listtype2 *next;
}llist;


//FDOM METHODS


//SECTION 1: BASE METHODS

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
GEN mat_nfcholesky(GEN nf, GEN A);

//SECTION 2: GEOMETRY METHODS

//BASIC LINE, CIRCLE, AND POINT OPERATIONS
GEN mat_eval(GEN M, GEN x);
GEN mobius_gp(GEN M, GEN c, long prec);

//DISTANCES/AREAS
GEN hdiscrandom(GEN R, long prec);
GEN hdist(GEN z1, GEN z2, long flag, long prec);

//FUNDAMENTAL DOMAIN COMPUTATION
GEN isometriccircle_mats(GEN g, GEN mats, GEN data, GEN (*gamtopsl)(GEN, GEN, long), GEN tol, long prec);
GEN isometriccircle_psu(GEN g, GEN tol, long prec);
GEN normalizedbasis(GEN G, GEN U, GEN mats, GEN gamid, GEN data, GEN (*gamtopsl)(GEN, GEN, long), GEN (*eltmul)(GEN, GEN, GEN), GEN (*eltinv)(GEN, GEN), int (*istriv)(GEN, GEN), GEN tol, long prec);
GEN normalizedboundary(GEN G, GEN mats, GEN id, GEN data, GEN (*gamtopsl)(GEN, GEN, long), GEN tol, long prec);
GEN psltopsu(GEN g, GEN p);
GEN psltopsu_mats(GEN g, GEN M);
GEN psltopsu_transmats(GEN p);
GEN reduceelt_givennormbound(GEN U, GEN g, GEN z, GEN data, GEN (*gamtopsl)(GEN, GEN, long), GEN (*eltmul)(GEN, GEN, GEN), GEN (*eltinv)(GEN, GEN), GEN tol, long prec);
GEN reducepoint(GEN U, GEN z, GEN gamid, GEN data, GEN (*gamtopsl)(GEN, GEN, long), GEN (*eltmul)(GEN, GEN, GEN), GEN tol, long prec);
GEN rootgeodesic_ud(GEN M, GEN mats, GEN tol, long prec);
GEN rootgeodesic_uhp(GEN M, GEN tol, long prec);

//FUNDAMENTAL DOMAIN OTHER COMPUTATIONS
GEN minimalcycles_bytype(GEN U, GEN gamid, GEN data, GEN (*eltmul)(GEN, GEN, GEN), GEN (*elttrace)(GEN, GEN), int (*istriv)(GEN, GEN));
GEN normalizedboundary_oosides(GEN U);
GEN rootgeodesic_fd(GEN U, GEN g, GEN gamid, GEN data, GEN (*gamtopsl)(GEN, GEN, long), GEN (*eltmul)(GEN, GEN, GEN), GEN (*eltinv)(GEN, GEN), GEN tol, long prec);
GEN presentation(GEN U, GEN gamid, GEN data, GEN (*eltmul)(GEN, GEN, GEN), GEN (*elttrace)(GEN, GEN), int (*istriv)(GEN, GEN));
GEN signature(GEN U, GEN gamid, GEN data, GEN (*eltmul)(GEN, GEN, GEN), GEN (*elttrace)(GEN, GEN), int (*istriv)(GEN, GEN));
GEN word(GEN U, GEN P, GEN g, GEN data, GEN (*gamtopsl)(GEN, GEN, long), GEN (*eltmul)(GEN, GEN, GEN), GEN (*eltinv)(GEN, GEN), int (*istriv)(GEN, GEN), GEN tol, long prec);

//HELPER METHODS
GEN deftol(long prec);


//SECTION 3: QUATERNIONIC METHODS

//QUATERNION METHODS
GEN algabsrednorm(GEN A, GEN p, GEN z1, GEN z2, long prec);
GEN algfdom(GEN A, GEN O, GEN p, int dispprogress, int dumppartial, GEN partialset, GEN constants, long prec);
GEN algfdomarea(GEN A, GEN O, int lessprec, long prec);
GEN algfdom_bestC(GEN A, GEN O, long prec);
GEN algfdomminimalcycles(GEN U, long prec);
GEN algfdompresentation(GEN U, long prec);
GEN algfdomreduce(GEN g, GEN U, GEN z, long prec);
GEN algfdomrootgeodesic(GEN g, GEN U, long prec);
GEN algfdomsignature(GEN U, long prec);
GEN algfdomword(GEN g, GEN P, GEN U, long prec);
GEN algmoreprec(GEN A, long increment, long prec);
GEN algmulvec(GEN A, GEN G, GEN L);
GEN algnormalizedbasis(GEN A, GEN O, GEN G, GEN p, long prec);
GEN algnormalizedboundary(GEN A, GEN O, GEN G, GEN p, long prec);
GEN algnormdisc(GEN A);
GEN algramifiedplacesf(GEN A);
GEN algsmallnorm1elts(GEN A, GEN O, GEN p, GEN C, GEN z1, GEN z2, int type, long prec);

//FUNDAMENTAL DOMAIN RETRIEVAL METHODS
GEN algfdomalg(GEN U);
GEN algfdom_get_alg(GEN U);
GEN algfdomorder(GEN U);
GEN algfdom_get_order(GEN U);

//HELPER METHODS
GEN algnorm_chol(GEN nf, GEN decomp, GEN x);
GEN qalg_absrednormqf(GEN Q, GEN mats, GEN z1, GEN z2, GEN normformpart, long prec);
GEN qalg_fdomarea(GEN Q, long computeprec, long prec);
GEN qalg_fdombestC(GEN Q, long prec);
GEN qalg_fdominitialize(GEN A, GEN O, GEN level, long prec);
GEN qalg_normform(GEN Q);
GEN qalg_smallnorm1elts_qfminim(GEN Q, GEN p, GEN C, GEN z1, GEN z2, long maxelts, GEN nfdecomp, GEN nformpart, long prec);
GEN qalg_smallnorm1elts_condition(GEN Q, GEN p, GEN C, GEN z1, GEN z2, long maxelts, GEN nform, GEN nformpart, long prec);

//BASIC OPERATIONS FOR NORMALIZED BASIS ET AL
GEN qalg_fdominv(GEN data, GEN x);
GEN qalg_fdomm2rembed(GEN data, GEN x, long prec);
GEN qalg_fdommul(GEN data, GEN x, GEN y);
GEN qalg_fdomtrace(GEN data, GEN x);
int qalg_istriv(GEN data, GEN x);

//3: SHALLOW RETRIEVAL METHODS
GEN qalg_get_alg(GEN Q);
GEN qalg_get_rams(GEN Q);
GEN qalg_get_varnos(GEN Q);
GEN qalg_get_roots(GEN Q);
GEN qalg_get_order(GEN Q);
GEN qalg_get_level(GEN Q);


//FDOM_EXTRA METHODS


//SECTION 1: GEOMETRY METHODS
void python_printarcs(GEN arcs, char *filename, int view, char *extrainput, long prec);
void python_plotviewer(char *input);
void python_printfdom(GEN U, char *filename, long prec);

//SECTION 2: QUATERNIONIC METHODS
GEN alg1ijktoalg(GEN A, GEN x);
GEN alg1ijktobasis(GEN A, GEN x);
GEN algalgto1ijk(GEN A, GEN x);
GEN algbasisto1ijk(GEN A, GEN x);
GEN algeichler_conj(GEN A, GEN x);
GEN algorderconj(GEN A, GEN x, GEN O);
GEN algorderdisc(GEN A, GEN O, int reduced, int factored);
GEN algorderlevel(GEN A, GEN O, int factored);
GEN algreduceddisc(GEN A);
GEN algshimura(GEN F, GEN D, long place, long maxcomptime, int allowswap);
GEN algshimura_ab(GEN F, GEN D, long place, long maxcomptime, int allowswap);
GEN algswapab(GEN A);
GEN smallalgebras(GEN F, long nwant, GEN Dmin, GEN Dmax, long maxcomptime, int allowswap);
GEN smallalgebras_area(GEN nf, GEN Amin, GEN Amax, int retD, int maxcomptime, int allowswap);

GEN algsmallelts(GEN A, GEN O, GEN nm, GEN p, GEN C, GEN z1, GEN z2, long prec);
GEN qalg_smallelts_qfminim(GEN Q, GEN nm, GEN p, GEN C, GEN z1, GEN z2, long maxelts, GEN nfdecomp, GEN nformpart, long prec);


//SECTION 3: PRODUCING DATA FOR MY PAPER

//OPTIMIZING THE VALUE OF C FOR ENUMERATION
GEN enum_bestC(GEN A, GEN p, GEN scale, long ntrials, long mintesttime, long prec);
GEN enum_bestC_range(GEN Aset, GEN p, GEN scale, long ntrials, long mintesttime, char *fname, int isArange, int compile, int WSL, long prec);
GEN enum_successrate(GEN A, GEN p, GEN C, long Ntests, GEN R, long prec);
GEN enum_successrate_range(GEN A, GEN p, GEN Cmin, GEN Cmax, long ntrials, long Ntests, GEN R, char *fname, int compile, int WSL, long prec);
GEN enum_time(GEN A, GEN p, GEN Cset, long mintesttime, long prec);
GEN enum_time_range(GEN A, GEN p, GEN Cmin, GEN Cmax, long ntrials, long mintesttime, char *fdata, int compile, int WSL, long prec);
long enum_timeforNelts(GEN A, GEN p, GEN C, long nelts, GEN R, int type, long prec);
void enum_timeforNelts_range(GEN A, GEN p, GEN Cmin, GEN Cmax, long ntrials, long nelts, char *fname, int type, int compile, int WSL, long prec);

//3: NUMBER OF ELEMENTS REQUIRED TO GENERATE ALGEBRA
GEN algfdom_nelts(GEN A, GEN p, GEN CNRdata, int type, long prec);

//REGRESSIONS & PLOTS
GEN OLS(GEN X, GEN y, int retrsqr);
GEN OLS_nointercept(GEN X, GEN y, int retrsqr);
GEN OLS_single(GEN x, GEN y, int retrsqr);
GEN rsquared(GEN X, GEN y, GEN fit);
void plot_compile(char *fname, int WSL);
