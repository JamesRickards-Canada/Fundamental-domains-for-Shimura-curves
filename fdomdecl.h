/*FDOM METHODS*/


/*SECTION 1: BASE METHODS*/

/*LISTS*/
GEN veclist_append(GEN v, long *vind, long *vlen, GEN x);
GEN vecsmalllist_append(GEN v, long *vind, long *vlen, long x);

/*SHORT VECTORS IN LATTICES*/
GEN mat_nfcholesky(GEN nf, GEN A);

/*SECTION 2: GEOMETRY METHODS*/

/*BASIC LINE, CIRCLE, AND POINT OPERATIONS*/
GEN mat_eval(GEN M, GEN x);
GEN mobius_gp(GEN M, GEN c, long prec);

/*DISTANCES/AREAS*/
GEN hdiscradius(GEN area, long prec);
GEN hdiscrandom(GEN R, long prec);
GEN hdiscrandom_arc(GEN R, GEN ang1, GEN ang2, long prec);
GEN hdist(GEN z1, GEN z2, long flag, long prec);

/*FUNDAMENTAL DOMAIN COMPUTATION*/
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

/*FUNDAMENTAL DOMAIN OTHER COMPUTATIONS*/
GEN minimalcycles_bytype(GEN U, GEN gamid, GEN data, GEN (*eltmul)(GEN, GEN, GEN), GEN (*elttrace)(GEN, GEN), int (*istriv)(GEN, GEN));
GEN normalizedboundary_oosides(GEN U);
GEN rootgeodesic_fd(GEN U, GEN g, GEN gamid, GEN data, GEN (*gamtopsl)(GEN, GEN, long), GEN (*eltmul)(GEN, GEN, GEN), GEN (*eltinv)(GEN, GEN), GEN tol, long prec);
GEN presentation(GEN U, GEN gamid, GEN data, GEN (*eltmul)(GEN, GEN, GEN), GEN (*elttrace)(GEN, GEN), int (*istriv)(GEN, GEN));
GEN signature(GEN U, GEN gamid, GEN data, GEN (*eltmul)(GEN, GEN, GEN), GEN (*elttrace)(GEN, GEN), int (*istriv)(GEN, GEN));
GEN word(GEN U, GEN P, GEN g, GEN data, GEN (*gamtopsl)(GEN, GEN, long), GEN (*eltmul)(GEN, GEN, GEN), GEN (*eltinv)(GEN, GEN), int (*istriv)(GEN, GEN), GEN tol, long prec);

/*HELPER METHODS*/
GEN deftol(long prec);
INLINE GEN
gc_0vec(pari_sp av){set_avma(av);return cgetg(1, t_VEC);}/*Resets avma and returns the vector []*/


/*SECTION 3: QUATERNIONIC METHODS*/

/*QUATERNION ALGEBRA NON-FDOM METHODS*/
int algelttype(GEN A, GEN g);
int algisorder(GEN A, GEN O);
GEN alg_get_ab(GEN A);
GEN algmoreprec(GEN A, long increment, long prec);
GEN algmulvec(GEN A, GEN G, GEN L);
GEN algnormdisc(GEN A);
GEN algorderdisc(GEN A, GEN O, int reduced, int factored);
GEN algorderlevel(GEN A, GEN O, int factored);
GEN algramifiedplacesf(GEN A);
GEN algreduceddisc(GEN A);

/*QUATERNION METHODS*/
GEN algabsrednorm(GEN A, GEN p, GEN z1, GEN z2, long prec);
GEN algfdom(GEN A, GEN O, int dispprogress, GEN constants, long prec);
GEN algfdomarea(GEN A, GEN O, int lessprec, long prec);
GEN algfdom_bestC(GEN A, GEN O, long prec);
GEN algfdomminimalcycles(GEN U, long prec);
GEN algfdompresentation(GEN U, long prec);
GEN algfdomreduce(GEN g, GEN U, GEN z, long prec);
GEN algfdomrootgeodesic(GEN g, GEN U, long prec);
GEN algfdomsignature(GEN U, long prec);
GEN algfdomword(GEN g, GEN P, GEN U, long prec);
GEN algnormalizedbasis(GEN A, GEN O, GEN G, GEN p, long prec);
GEN algnormalizedboundary(GEN A, GEN O, GEN G, GEN p, long prec);
GEN algsmallnorm1elts(GEN A, GEN O, GEN p, GEN C, GEN z1, GEN z2, int type, long prec);

/*FUNDAMENTAL DOMAIN RETRIEVAL METHODS*/
GEN algfdomalg(GEN U);
GEN algfdomorder(GEN U);

/*HELPER METHODS*/
GEN algnorm_chol(GEN nf, GEN decomp, GEN x);
GEN qalg_absrednormqf(GEN Q, GEN mats, GEN z1, GEN z2, GEN normformpart, long prec);
GEN qalg_fdomarea(GEN Q, long computeprec, long prec);
GEN qalg_fdombestC(GEN Q, long prec);
GEN qalg_fdominitialize(GEN A, GEN O, long prec);
GEN qalg_normform(GEN Q);
GEN qalg_smallnorm1elts_qfminim(GEN Q, GEN p, GEN C, GEN z1, GEN z2, long maxelts, GEN nfdecomp, GEN nformpart, long prec);
GEN qalg_smallnorm1elts_condition(GEN Q, GEN p, GEN C, GEN z1, GEN z2, long maxelts, GEN nform, GEN nformpart, long prec);

/*BASIC OPERATIONS FOR NORMALIZED BASIS ET AL*/
GEN qalg_fdominv(GEN data, GEN x);
GEN qalg_fdomm2rembed(GEN data, GEN x, long prec);
GEN qalg_fdommul(GEN data, GEN x, GEN y);
GEN qalg_fdomtrace(GEN data, GEN x);
int qalg_istriv(GEN data, GEN x);

/*3: INLINE SHALLOW RETRIEVAL METHODS*/
/*Shallow method to return the algebra*/
INLINE GEN
qalg_get_alg(GEN Q){return gel(Q, 1);}
/*Shallow method to get the ramification*/
INLINE GEN
qalg_get_rams(GEN Q){return gel(Q, 2);}
/*Shallow method to return the variable numbers of K and L*/
INLINE GEN
qalg_get_varnos(GEN Q){return gel(Q, 3);}
/*Shallow method to return the roots of K and L.*/
INLINE GEN
qalg_get_roots(GEN Q){return gel(Q, 4);}
/*Shallow method to return the order*/
INLINE GEN
qalg_get_order(GEN Q){return gel(Q, 5);}
/*Shallow method to return the level of the order*/
INLINE GEN
qalg_get_level(GEN Q){return gel(Q, 6);}
/*Shallow method to return the tolerance*/
INLINE GEN
qalg_get_tol(GEN Q){return gel(Q, 7);}
/*Returns the qalg*/
INLINE GEN algfdom_get_qalg(GEN U){return gel(U, 9);}
/*Shallow version of algfdomalg*/
INLINE GEN algfdom_get_alg(GEN U){return qalg_get_alg(algfdom_get_qalg(U));}
/*Shallow version of algfdomorder*/
INLINE GEN
algfdom_get_order(GEN U){return qalg_get_order(algfdom_get_qalg(U));}

/*EXTRA*/
GEN alg_changeorder(GEN al, GEN ord);

/*FDOM_EXTRA METHODS*/


/*SECTION 1: GEOMETRY METHODS*/
void fdom_latex(GEN U, char *filename, int boundcircle, int compile, int open, long prec);
void python_printarcs(GEN arcs, char *filename, int view, char *extrainput, long prec);
void python_plotviewer(char *input);
void python_printfdom(GEN U, char *filename, long prec);

/*SECTION 2: QUATERNIONIC METHODS*/
GEN alg1ijktoalg(GEN A, GEN x);
GEN alg1ijktobasis(GEN A, GEN x);
GEN algalgto1ijk(GEN A, GEN x);
GEN algbasisto1ijk(GEN A, GEN x);
void algeichler1(GEN A, GEN I);
GEN algeichler2(GEN A, GEN O);
GEN algeichler_conj(GEN A, GEN x);
GEN algorderconj(GEN A, GEN x, GEN O);
GEN algshimura(GEN F, GEN D, long place, long maxcomptime, int allowswap);
GEN algshimura_ab(GEN F, GEN D, long place, long maxcomptime, int allowswap);
GEN algswapab(GEN A);
GEN smallalgebras(GEN F, long nwant, GEN Dmin, GEN Dmax, long maxcomptime, int allowswap);
GEN smallalgebras_area(GEN nf, GEN Amin, GEN Amax, int retD, int maxcomptime, int allowswap);

GEN algsmallelts(GEN A, GEN O, GEN nm, GEN p, GEN C, GEN z1, GEN z2, long prec);
GEN qalg_smallelts_qfminim(GEN Q, GEN nm, GEN p, GEN C, GEN z1, GEN z2, long maxelts, GEN nfdecomp, GEN nformpart, long prec);


/*SECTION 3: PRODUCING DATA FOR MY PAPER*/

/*OPTIMIZING THE VALUE OF C FOR ENUMERATION*/
GEN enum_bestC(GEN A, GEN p, GEN scale, long ntrials, long mintesttime, long prec);
GEN enum_bestC_range(GEN Aset, GEN p, GEN scale, long ntrials, long mintesttime, char *fname, int isArange, int compile, int WSL, long prec);
GEN enum_successrate(GEN A, GEN p, GEN C, long Ntests, GEN R, long prec);
GEN enum_successrate_range(GEN A, GEN p, GEN Cmin, GEN Cmax, long ntrials, long Ntests, GEN R, char *fname, int compile, int WSL, long prec);
GEN enum_time(GEN A, GEN p, GEN Cset, long mintesttime, long prec);
GEN enum_time_range(GEN A, GEN p, GEN Cmin, GEN Cmax, long ntrials, long mintesttime, char *fdata, int compile, int WSL, long prec);
long enum_timeforNelts(GEN A, GEN p, GEN C, long nelts, GEN R, int type, long prec);
void enum_timeforNelts_range(GEN A, GEN p, GEN Cmin, GEN Cmax, long ntrials, long nelts, char *fname, int type, int compile, int WSL, long prec);

/*3: NUMBER OF ELEMENTS REQUIRED TO GENERATE ALGEBRA*/
GEN algfdom_nelts(GEN A, GEN p, GEN CNRdata, int type, long prec);

/*REGRESSIONS & PLOTS*/
GEN OLS(GEN X, GEN y, int retrsqr);
GEN OLS_nointercept(GEN X, GEN y, int retrsqr);
GEN OLS_single(GEN x, GEN y, int retrsqr);
GEN rsquared(GEN X, GEN y, GEN fit);
void plot_compile(char *fname, int WSL);
