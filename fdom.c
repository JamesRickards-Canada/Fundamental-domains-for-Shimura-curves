//Fundamental domain computation

//INCLUSIONS

#ifndef PARILIB
#define PARILIB
#include <pari/pari.h>
#endif

#ifndef METHDECL
#define METHDECL
#include "fdomdecl.h"
#endif

//DEFINITIONS

//The length (lg, so technically length+1) of a circle/line and arc/segment, and a normalized boundary
#define CIRCLEN 4
#define ARCLEN 9
#define NORMBOUND 9

//STATIC DECLARATIONS

//1: SHORT VECTORS IN LATTICES
static GEN quadraticintegernf(GEN nf, GEN A, GEN B, GEN C, long prec);
static GEN smallvectors_cholesky(GEN Q, GEN C, long maxelts, GEN condition, long prec);
static GEN smallvectors_nfcondition(GEN A, GEN C, long maxelts, GEN condition, long prec);

//2: BASIC LINE, CIRCLE, AND POINT OPERATIONS
static GEN arc_init(GEN c, GEN p1, GEN p2, int dir, long prec);
static GEN arc_midpoint(GEN c, GEN p1, GEN p2, GEN tol, long prec);
static GEN circle_angle(GEN c1, GEN c2, GEN p, GEN tol, long prec);
static GEN circle_fromcp(GEN cent, GEN p, long prec);
static GEN circle_fromppp(GEN p1, GEN p2, GEN p3, GEN tol, long prec);
static GEN circle_tangentslope(GEN c, GEN p, long prec);
static GEN line_fromsp(GEN s, GEN p);
static GEN line_frompp(GEN p1, GEN p2, GEN tol, long prec);
static GEN midpoint(GEN p1, GEN p2);
static GEN mobius(GEN M, GEN c, GEN tol, long prec);
static GEN mobius_arcseg(GEN M, GEN c, int isarc, GEN tol, long prec);
static GEN mobius_circle(GEN M, GEN c, GEN tol, long prec);
static GEN mobius_line(GEN M, GEN l, GEN tol, long prec);
static GEN perpbis(GEN p1, GEN p2, GEN tol, long prec);
static GEN radialangle(GEN c, GEN p, GEN tol, long prec);
static GEN slope(GEN p1, GEN p2, GEN tol, long prec);

//2: INTERSECTION OF LINES/CIRCLES
static GEN arc_int(GEN c1, GEN c2, GEN tol, long prec);
static GEN arcseg_int(GEN c, GEN l, GEN tol, long prec);
static GEN circle_int(GEN c1, GEN c2, GEN tol, long prec);
static GEN circleline_int(GEN c, GEN l, GEN tol, long prec);
static GEN line_int(GEN l1, GEN l2, GEN tol, long prec);
static int onarc(GEN c, GEN p, GEN tol, long prec);
static int onseg(GEN l, GEN p, GEN tol, long prec);

//2: DISTANCES
static GEN hpolygon_area(GEN circles, GEN vertices, GEN tol, long prec);

//2: FUNDAMENTAL DOMAIN COMPUTATION
static GEN edgepairing(GEN U, GEN tol, int rboth, long prec);
static GEN normalizedbasis_shiftpoint(GEN c, GEN r, int initial, long prec);
static GEN normalizedboundary_append(GEN Ubase, GEN G, GEN mats, GEN id, GEN tol, long prec);
static GEN normalizedboundary_givenU(GEN Ubase, GEN G, GEN mats, GEN id, GEN *data, GEN (*gamtopsl)(GEN *, GEN, long), GEN tol, long prec);
static GEN normalizedboundary_givencircles(GEN G, GEN mats, GEN id, GEN tol, long prec);
static long normalizedboundary_outside(GEN U, GEN z, GEN tol, long prec);
static GEN normalizedboundary_sideint(GEN U, GEN c, int start, GEN tol, long prec);
static GEN psl_roots(GEN M, GEN tol, long prec);

//2: FUNDAMENTAL DOMAIN OTHER COMPUTATIONS
static GEN minimalcycles(GEN pair);

//2: GEOMETRIC HELPER METHODS
static GEN anglediff(GEN ang, GEN bot, GEN tol, long prec);
static GEN atanoo(GEN x, long prec);
static int gcmp_strict(void *data, GEN x, GEN y);
static int geom_check(GEN c);
static GEN shiftangle(GEN ang, GEN bot, GEN tol, long prec);
static long tolcmp(GEN x, GEN y, GEN tol, long prec);
static int tolcmp_sort(void *data, GEN x, GEN y);
static int toleq(GEN x, GEN y, GEN tol, long prec);

//3: FUNDAMENTAL DOMAIN COMPUTATION
static GEN qalg_fdom(GEN Q, GEN p, int dispprogress, int dumppartial, GEN partialset, GEN C, GEN R, GEN passes, int type, GEN tol, long prec);

//3: HELPER METHODS
static long algsplitoo(GEN A);
static GEN qalg_normform_givenbasis(GEN Q, GEN basis);
static GEN qalg_basis_conj(GEN Q, GEN x);



//SECTION 1: BASE METHODS



//INFINITY 


//Divides a and b, and allows for oo and division by 0. Returns oo for 0/0.
GEN divoo(GEN a, GEN b){//No garbage collection necessary
  if(gequal0(b)){//b=0
    if(gcmpgs(a,0)>=0) return mkoo();
    return mkmoo();
  }
  if(typ(a)==t_INFINITY){//a=+/-oo
	if(gsigne(a)==gsigne(b)) return mkoo();
	return mkmoo();
  }
  if(typ(b)==t_INFINITY) return gen_0;
  return gdiv(a,b);
}



//LISTS


//List of GENs

//Frees the memory pari_malloc'ed by glist
void glist_free(glist *l){
  glist *temp=l;
  while(l!=NULL){
	temp=l;
	l=l->next;
	pari_free(temp);
  }
}

//Removes the last element of the glist and returns it without copying
GEN glist_pop(glist **head_ref){
  glist *temp=*head_ref;
  GEN x=temp->data;
  *head_ref=temp->next;
  pari_free(temp);
  return x;
}

//Put an element at the start of the glist
void glist_putstart(glist **head_ref, GEN new_data){
  glist *new_elt = (glist*)pari_malloc(sizeof(glist)); 
  new_elt->data = new_data; 
  new_elt->next = *head_ref; 
  *head_ref = new_elt; 

}

//dir=1 means forward, dir=-1 means backwards. Returns the list as a vector, makes a clean copy. This also frees the list, but we also need to clean up the list data at the list creation location. The passed in pointer to l should NOT be used as it no longer points to a valid address.
GEN glist_togvec(glist *l, long length, int dir){
	glist *lcopy=l;
	GEN rvec=cgetg(length+1, t_VEC);
	if(dir==1){
	  long lind=1;
	  while(l!=NULL && lind<=length){
		  gel(rvec,lind)=gcopy(l->data);
		  l=l->next;
		  lind++;
	  }
	  if(lind<=length){//Couldn't finish.
	    pari_err(e_MISC,"List length is too long");
	  }
	}
	else{
      long lind=length;
	  while(l!=NULL && lind>0){
		gel(rvec,lind)=gcopy(l->data);
		l=l->next;
		lind--;
	  }
	  if(lind>0){//Couldn't finish.
	    pari_err(e_MISC,"List length is too long");
	  }
	}
	glist_free(lcopy);
	return rvec;
}

//Appends l to the end of v, returning a clean copy. dir=-1 means forward, dir=-1 backward. This also frees the list, but we also need to clean up the list data at the list creation location. The passed in pointer to l should NOT be used as it no longer points to a valid address.
GEN glist_togvec_append(glist *l, GEN v, long length, int dir){
	glist *lcopy=l;
	long vlen=lg(v), rveclen=length+vlen;
	GEN rvec=cgetg(rveclen, t_VEC);
	for(long i=1;i<vlen;i++) gel(rvec, i)=gcopy(gel(v, i));//Copying v
	if(dir==1){
	  long lind=vlen;
	  while(l!=NULL && lind<rveclen){
		  gel(rvec,lind)=gcopy(l->data);
		  l=l->next;
		  lind++;
	  }
	  if(lind<rveclen){//Couldn't finish.
	    pari_err(e_MISC,"List length is too long");
	  }
	}
	else{
      long lind=rveclen-1;
	  while(l!=NULL && lind>=vlen){
		gel(rvec,lind)=gcopy(l->data);
		l=l->next;
		lind--;
	  }
	  if(lind>=vlen){//Couldn't finish.
	    pari_err(e_MISC,"List length is too long");
	  }
	}
	glist_free(lcopy);
	return rvec;
}


//List of longs

//Frees the memory pari_malloc'ed by llist
void llist_free(llist *l){
  llist *temp=l;
  while(l!=NULL){
	temp=l;
	l=l->next;
	pari_free(temp);
  }
}

//Removes the last element of the llist and returns it
long llist_pop(llist **head_ref){
  llist *temp=*head_ref;
  long x=temp->data;
  *head_ref=temp->next;
  pari_free(temp);
  return x;
}

//Put an element at the start of the llist
void llist_putstart(llist **head_ref, long new_data){
  llist *new_elt = (llist*)pari_malloc(sizeof(llist)); 
  new_elt->data = new_data; 
  new_elt->next = *head_ref; 
  *head_ref = new_elt; 
}

//dir=1 means forward, dir=-1 means backwards. Returns the list as a vector. This also frees the list. The passed in pointer to l should NOT be used as it no longer points to a valid address.
GEN llist_togvec(llist *l, long length, int dir){//No garbage collection necessary with longs!
  llist *lcopy=l;
  GEN rvec=cgetg(length+1, t_VEC);
  if(dir==1){
    long lind=1;
    while(l!=NULL && lind<=length){
	  gel(rvec,lind)=stoi(l->data);
	  l=l->next;
	  lind++;
    }
    if(lind<=length){//Couldn't finish.
      pari_err(e_MISC,"List length is too long");
    }
  }
  else{
    long lind=length;
    while(l!=NULL && lind>0){
      gel(rvec,lind)=stoi(l->data);
	  l=l->next;
	  lind--;
    }
    if(lind>0){//Couldn't finish.
      pari_err(e_MISC,"List length is too long");
    }
  }
  llist_free(lcopy);
  return(rvec);
}

//dir=1 means forward, dir=-1 means backwards. Returns the list as a VECSMALL. This also frees the list. The passed in pointer to l should NOT be used as it no longer points to a valid address.
GEN llist_tovecsmall(llist *l, long length, int dir){//No garbage collection necessary with longs!
  llist *lcopy=l;
  GEN rvec=cgetg(length+1, t_VECSMALL);
  if(dir==1){
    long lind=1;
    while(l!=NULL && lind<=length){
	  rvec[lind]=l->data;
	  l=l->next;
	  lind++;
    }
    if(lind<=length){//Couldn't finish.
      pari_err(e_MISC,"List length is too long");
    }
  }
  else{
    long lind=length;
    while(l!=NULL && lind>0){
      rvec[lind]=l->data;
	  l=l->next;
	  lind--;
    }
    if(lind>0){//Couldn't finish.
      pari_err(e_MISC,"List length is too long");
    }
  }
  llist_free(lcopy);
  return(rvec);
}


//SHORT VECTORS IN LATTICES


//Computes the Cholesky decomposition of A (nxn symmetric matrix) whose coefficients are in nf. Returns the nxn matrix B so that x^TAx is expressible as sum(i=1..n)b_ii(x_i+sum(j=i+1..n)b_ijxj)^2. Note that A does not need to be positive definite, but it CANNOT represent 0 (as then we cannot write the quadratic form as a sum/difference of squares). The result in this case will be meaningless.
GEN mat_nfcholesky(GEN nf, GEN A){
  pari_sp top=avma;
  long n=lg(A)-1;//A is nxn
  GEN M=gcopy(A);//Will be manipulating the entries, so need to copy A.
  for(long i=1;i<n;i++){
	if(gequal0(gcoeff(M, i, i))){
	  for(long j=i+1;j<=n;j++){
	    gcoeff(M, j, i)=gcopy(gcoeff(M, i, j));//M[j,i]=M[i,j]
	  }
	}
	else{
	  for(long j=i+1;j<=n;j++){
	    gcoeff(M, j, i)=gcopy(gcoeff(M, i, j));//M[j,i]=M[i,j]
		pari_CATCH(CATCH_ALL){
		  gcoeff(M, i, j)=gen_0;
		}
		pari_TRY{//Only issue is a current bug with pari when M[i,j]=0
	      gcoeff(M, i, j)=nfdiv(nf, gcoeff(M, i, j), gcoeff(M, i, i));//M[i,j]=M[i,j]/M[i,i]
		}
		pari_ENDCATCH
	  }
	}
	for(long j=i+1;j<=n;j++){
	  for(long k=j;k<=n;k++) gcoeff(M, j, k)=nfsub(nf, gcoeff(M, j, k), nfmul(nf, gcoeff(M, j, i), gcoeff(M, i, k)));//M[j,k]=M[j,k]-M[j,i]*M[i,k];
	}
  }
  GEN ret=cgetg_copy(M, &n);//M stores the coeff, but we should delete the lower diagonal
  for(long i=1;i<n;i++){//Column i
    gel(ret, i)=cgetg(n, t_COL);
	for(long j=1;j<=i;j++) gcoeff(ret, j, i)=gcopy(gcoeff(M, j, i));
	for(long j=i+1;j<n;j++) gcoeff(ret, j, i)=gen_0;
  }
  return gerepileupto(top, ret);
}

//Solves Ax^2+Bx+C=0 in the integers (A, B, C belong to nf) and returns the vector of solutions.
static GEN quadraticintegernf(GEN nf, GEN A, GEN B, GEN C, long prec){
  pari_sp top=avma;
  if(gequal0(A)){//Actually a linear. This case occurs when dealing with small vectors in the quaternion algebra ramified nowhere.
	if(gequal0(B)) return cgetg(1, t_VEC);//We say there are no soln's when A=B=0
	GEN x=gneg(nfdiv(nf, C, B));//Solution
	x=lift(basistoalg(nf, x));
	if(typ(x)==t_INT) return gerepilecopy(top, mkvec(x));//Solution!
	avma=top;
	return cgetg(1, t_VEC);
  }
  GEN Ap=nfeltembed(nf, A, gen_1, prec);//Doesn't matter which place, since if we have an integer solution iff it works for all places
  GEN Bp=nfeltembed(nf, B, gen_1, prec);
  GEN Cp=nfeltembed(nf, C, gen_1, prec);
  
  GEN disc=gsub(gsqr(Bp), gmulsg(4, gmul(Ap, Cp)));//B^2-4AC
  GEN rt=real_i(gsqrt(disc, prec));//We use greal in case disc==0 but it is approximated as -E^(-large).
  GEN mBp=gneg(Bp), twoAp=gmulsg(2, Ap);//-B, 2A

  GEN r1=ground(gdiv(gadd(mBp, rt), twoAp));//We have enough precision that the rounded answer will be correct if there is an integral solution
  GEN r2=ground(gdiv(gsub(mBp, rt), twoAp));
  GEN rposs;//The possible roots
  if(equalii(r1, r2)) rposs=mkvec(r1);
  else rposs=mkvec2(r1, r2);
  //Now we plug back in and check.
  GEN rts=vectrunc_init(3), res, r;//At most 2 roots
  for(long i=1;i<lg(rposs);i++){
	r=gel(rposs, i);
	res=nfadd(nf, nfmul(nf, nfadd(nf, nfmul(nf, A, r), B), r), C);//Plug it in
	if(gequal0(res)) vectrunc_append(rts, r);//gequal0(res)=1 no matter what representation it is in.
  }
  return gerepilecopy(top, rts);
}

//We follow the article "Improved Methods for Calculating Vectors of Short Length in a Lattice, Including Complexity Analysis" by Fincke and Pohst (Mathematics of Computation, Vol. 44, No. 170 (Apr., 1985), pp 463-471

//Follows Algorithm 2.12 in Fincke-Pohst, where we pass in a condition, which is [nf, M, n] where x^T*M*x=n must also be satisfied, with M and n being elements of the number field nf. Note that a good portion of the code in this method is copied from fincke_pohst in bibli1.c.
static GEN smallvectors_nfcondition(GEN A, GEN C, long maxelts, GEN condition, long prec){  
  pari_sp top=avma;
  long l=lg(A), newprec=prec;//n+1
  GEN U=lllfp(A, 0.75, LLL_GRAM | LLL_IM);
  if(lg(U) != lg(A)) return NULL;
  GEN R=qf_apply_RgM(A, U);
  long rprec=gprecision(R);
  if(rprec) newprec=rprec;
  else{
    newprec = DEFAULTPREC + nbits2extraprec(gexpo(R));
    if (newprec<prec) newprec=prec;
  }
  R=qfgaussred_positive(R);
  if (!R) return NULL;//In case there was an issue with R.
  for(long i=1; i<l; i++){
    GEN s = gsqrt(gcoeff(R,i,i), newprec);
    gcoeff(R,i,i) = s;
    for(long j=i+1;j<l;j++) gcoeff(R,i,j) = gmul(s, gcoeff(R,i,j));
  }
  /* now R~*R = A in LLL basis */
  GEN Rinv = RgM_inv_upper(R);
  if(!Rinv) return NULL;
  GEN Rinvtrans = shallowtrans(Rinv);
  GEN V = lll(Rinvtrans);
  if(lg(V)!=lg(Rinvtrans)) return NULL;
  Rinvtrans = RgM_mul(Rinvtrans, V);
  V = ZM_inv(shallowtrans(V), NULL);
  R = RgM_mul(R, V);
  U = U? ZM_mul(U, V): V;

  l = lg(R);
  GEN vnorm = cgetg(l, t_VEC);
  for(long j=1; j<l; j++) gel(vnorm, j) = gnorml2(gel(Rinvtrans, j));
  GEN rperm = cgetg(l,t_MAT);
  GEN uperm = cgetg(l,t_MAT);
  GEN perm = indexsort(vnorm);
  for(long i=1; i<l; i++) {uperm[l-i] = U[perm[i]]; rperm[l-i] = R[perm[i]]; }
  U = uperm;
  R = rperm;
  GEN res=cgetg(1, t_VEC);
  pari_CATCH(e_PREC) { }
  pari_TRY {
    GEN q = gaussred_from_QR(R, gprecision(vnorm));
    if (!q) pari_err_PREC("smallvectors_nfcondition");
    GEN nf=gel(condition, 1);
    GEN newcond=cgetg(4, t_VEC);
    gel(newcond, 1)=nf;
    gel(newcond, 2)=nfM_mul(nf, shallowtrans(U), nfM_mul(nf, gel(condition, 2), U));//M->U^T*M*U
    gel(newcond, 3)=gel(condition, 3);
    res = smallvectors_cholesky(q, C, maxelts, newcond, prec);//The small entries
  } pari_ENDCATCH;
  GEN ret=cgetg_copy(res, &l);
  for(long i=1;i<l;i++) gel(ret, i)=ZM_ZC_mul(U, gel(res, i));
  return gerepilecopy(top, ret);
}

//Q is the Cholesky decomposition of a matrix A, this computes all vectors x such that x^T*A*x<=C (i.e. Q(x)<=C_2 where Q(x)=sum(i=1..n)q_ii(x_i+sum(j=i+1..n)q_ijxj)^2. Algorithm 2.8 of Fincke Pohst. We also pass in a condition, which is [nf, M, n] where x^T*M*x=n must also be satisfied (M and n live in the number field nf). (the point: M gives an indefinite norm condition on the vector, and A combines this norm with other info to make a positive definite form. We use the condition when finding small norm 1 elements of a quaternion algebra.) This does NOT work if the condition can represent 0, as we rely on writing conditon as the sum/difference of squares.
static GEN smallvectors_cholesky(GEN Q, GEN C, long maxelts, GEN condition, long prec){
  pari_sp tiptop=avma, top, mid;
  GEN nf=gel(condition, 1);
  GEN condchol=mat_nfcholesky(nf, gel(condition, 2));//Cholesky of the condition.

  top=avma;
  long np1=lg(Q), n=np1-1;//Number of variables+1 and n
  GEN T=zerovec(n);//Stores the ''tail'' of x, as we work from the back to the front (x_m to x_1), i.e. C-sum(j=i+1)^n of qii*(x_i+U_i)^2
  GEN U=zerovec(n);//U[i] will store the sum of j=i+1 to n of q_{ij}x_j
  GEN UB=zerovec(n);//UB represents the upper bound for x_i
  
  GEN Tcond=zerovec(n), Aco, Bco, Cco;//T, but for the condition
  GEN Ucond=zerovec(n);//U, but for the condition
  
  GEN x=zerocol(n), Z;//x represents the solution
  long i=np1-1, count=0, totcount=0;//i represents the current index, initially set to n. count=the number of solutions
  gel(T, n)=gcopy(C);//initialize T[n]=C
  gel(U, n)=gen_0;//Clearly U[n]=0. U[i]=sum(j=i+1,n,q[i,j]*x[j]);
  
  gel(Tcond, n)=nfsub(nf, gen_0, gel(condition, 3));
  gel(Ucond, n)=gen_0;//Clearly 0

  int step=2;//Represents the current step of algorithm 2.8
  int xpass0=0;
  GEN x1sols;
  glist *S=NULL;//Pointer to the list start
  GEN v=cgetg(1, t_VEC);//The list, is used for garbage collection partway through
  while(step>0){
	if(gc_needed(top, 1)){
	  mid=avma;
	  v=glist_togvec_append(S, v, count, 1);
	  count=0;
	  S=NULL;
	  T=gcopy(T);
	  U=gcopy(U);
	  UB=gcopy(UB);
	  x=gcopy(x);
	  Tcond=gcopy(Tcond);
	  Ucond=gcopy(Ucond);
	  gerepileallsp(top, mid, 7, &v, &T, &U, &UB, &x, &Tcond, &Ucond);
	}
	if(step==2){
	  Z=gsqrt(gabs(gdiv(gel(T, i), gcoeff(Q, i, i)), prec), prec);//The inner square root should be positive always, but could run into issue if T=0 and rounding puts it <0. Z=sqrt(T[i]/Q[i,i])
	  gel(UB, i)=gfloor(gsub(Z, gel(U, i)));//UB[i]=floor(Z-U[i]);
	  gel(x, i)=gsubgs(gceil(gneg(gadd(Z, gel(U, i)))), 1);//x[i]=ceil(-Z-U[i])-1;
	  step=3;
	}
    if(step==3){
	  gel(x, i)=gaddgs(gel(x, i), 1);//x[i]=x[i]+1
	  if(gcmp(gel(x, i), gel(UB, i))<=0) step=5;//If x[i]<=UB[i], goto step 5
	  else step=4; //If x[i]>UB[i], goto step 4
	}
	if(step==4){
	  i++;
	  step=3;
	  continue;//May as well go back to start
	}
	if(step==5){
	  i--;
	  gel(U, i)=gen_0;
	  for(long j=i+1;j<np1;j++) gel(U, i)=gadd(gel(U, i), gmul(gcoeff(Q, i, j), gel(x, j)));//U[i]=sum(j=i+1,n,q[i,j]*x[j]);
	  gel(Ucond, i)=gen_0;
	  if(!gequal0(gcoeff(condchol, i, i))){//ith row of condchol is non-zero, so something to add to Ucond
	    for(long j=i+1;j<np1;j++) gel(Ucond, i)=nfadd(nf, gel(Ucond, i), nfmul(nf, gcoeff(condchol, i, j), gel(x, j)));
	  }
	  if(!gequal0(gcoeff(condchol, i+1, i+1))){//i+1th row of condchol is non-zero, so something to add to Tcond
	    gel(Tcond, i)=nfadd(nf, gel(Tcond, i+1), nfmul(nf, gcoeff(condchol, i+1, i+1), nfsqr(nf, nfadd(nf, gel(x, i+1), gel(Ucond, i+1)))));
	  }
	  else gel(Tcond, i)=gcopy(gel(Tcond, i+1));//Same
	  gel(T, i)=gsub(gel(T, i+1), gmul(gcoeff(Q, i+1, i+1), gsqr(gadd(gel(x, i+1), gel(U, i+1)))));//T[i]=T[i+1]-q[i+1,i+1]*(x[i+1]+U[i+1])^2;
	  if(i==1){step=6;}//We have a condition to deal with!
	  else{//Go back now
		step=2;
		continue;
	  }
	}
	if(step==6){//Dealing with extra condtions
	  step=3;//The next step, if we make it out.
	  i=2;//We go back to i=2, if we make it out
	  
	  GEN q11U1=nfmul(nf, gcoeff(condchol, 1, 1), gel(Ucond, 1));
	  Aco=gcoeff(condchol, 1, 1);
	  Bco=nfmul(nf, q11U1, gen_2);
	  Cco=nfadd(nf, gel(Tcond, 1), nfmul(nf, q11U1, gel(Ucond, 1)));//Ax_1^2+Bx_1+C=0 is necessary

	  gel(x, 1)=gen_0;
	  x1sols=quadraticintegernf(nf, Aco, Bco, Cco, prec);//Tcond_1+q_11(x_1+Ucond_1)^2
	  if(gequal0(x)) xpass0=1;//This is the last check
	  for(long j=1;j<lg(x1sols);j++){//We don't actually check that Q(x)<=C, as what we really care about are norm 1 vectors, and if we happen to discover one slightly outside of the range, there is no issue.
		if(xpass0 && signe(gel(x1sols, j))!=-1) continue;//x is 0 (except the first coefficient), so the first coefficent has to be negative.
		gel(x, 1)=gel(x1sols, j);//Now we are good, all checks out.
		glist_putstart(&S, gcopy(x));
	    count++;
		if(maxelts!=0){
		  totcount++;//We can't use count, since this resets if we garbage collect
		  if(totcount<maxelts) continue;
		  return gerepileupto(tiptop, glist_togvec_append(S, v, count, -1));//We hit the maximal number of return elements
	    }
	  }
	  if(xpass0){step=0;continue;}//Game over, we are done!
	  continue;
	}
  }
  return gerepileupto(tiptop, glist_togvec_append(S, v, count, -1));
}



//SECTION 2: GEOMETRY METHODS



//Circle is stored as [centre, radius, 0], where 0 means a bona fide cicle (and not a line)
//Circle arc is [centre, radius, start pt, end pt, start angle, end angle, dir, 0]. It is the arc counterclockwise from startpt to endpt, and dir=1 means oriented counterclockwise, and dir=-1 means oriented clockwise. This can also be left uninitialized if arc is not oriented. The final 0 represents that we have an arc and not a segment.
//Line is stored as [slope, intercept, 1], where the line is y=slope*x+intercept unless slope=oo, where it is x=intercept instead, and the 1 just means a bona fide line (not circle).
//Line segment is stored as [slope, intercept, startpt, endpt, 0, ooendptor, dir, 1] (the extra 0 is to format it the same as a circle arc). The final 1 is to signal a segment. dir=1 means we go from startpt to endpt in the upper half plane, and -1 means we go through infinity (only when neither endpoint is infinite). If one endpoint is oo, ooendptor stores which way we get to it. If ooendptor=1, this means the segment travels either vertically up or right, and -1 means the arc is vertically down or left.

//GEN tol -> The tolerance.


//BASIC LINE, CIRCLE, AND POINT OPERATIONS


//Creates the arc on circle c going from p1 to p2 counterclockwise. if dir=1, the arc is oriented counterclockwise, else it is clockwise (i.e. clockwise from p2 to p1). If dir=0, we take it to be unoriented
static GEN arc_init(GEN c, GEN p1, GEN p2, int dir, long prec){
  pari_sp top=avma;
  GEN ang2=radialangle(c, p2, gen_0, prec);//No tolerance need
  GEN arc=cgetg(ARCLEN, t_VEC);
  gel(arc, 1)=gcopy(gel(c, 1));
  gel(arc, 2)=gcopy(gel(c, 2));
  gel(arc, 3)=gcopy(p1);
  gel(arc, 4)=gcopy(p2);
  gel(arc, 5)=radialangle(c, p1, gen_0, prec);//start angle
  gel(arc, 6)=shiftangle(ang2, gel(arc, 5), gen_0, prec);//end angle; no need for tolerance.
  if(dir==1) gel(arc, 7)=gen_1;
  else if(dir==0) gel(arc, 7)=gen_0;
  else gel(arc, 7)=gen_m1;
  gel(arc, 8)=gen_0;
  return gerepileupto(top, arc);
}

//Returns the midpoint of the arc between p1 and p2 (counterclockwise) on c.
static GEN arc_midpoint(GEN c, GEN p1, GEN p2, GEN tol, long prec){
  pari_sp top=avma;
  GEN pts=circleline_int(c, perpbis(p1, p2, tol, prec), tol, prec);
  GEN a1=radialangle(c, p1, gen_0, prec);
  GEN a2=shiftangle(radialangle(c, p2, gen_0, prec), a1, gen_0, prec);//No tolerance concerns
  GEN angint1=shiftangle(radialangle(c, gel(pts, 1), gen_0, prec), a1, gen_0, prec);//the angle formed by pts[1] to c with base a1. Again, no tolerance need.
  if(gcmp(a2, angint1)==1) return gerepilecopy(top, gel(pts, 1));//No need for tolerance, as if this is an issue our points p1 and p2 would be equal up to tolerance.
  return gerepilecopy(top, gel(pts, 2));
}

//Returns the angle of intersection between circles c1 and c2 which intersect at p, with the angle formed by rotating the tangent to c1 at p counterclockwise to the tangent to c2 at p.
static GEN circle_angle(GEN c1, GEN c2, GEN p, GEN tol, long prec){
  pari_sp top=avma;
  GEN s1=circle_tangentslope(c1, p, prec);
  GEN s2=circle_tangentslope(c2, p, prec);
  GEN ang=anglediff(atanoo(s2, prec), atanoo(s1, prec), tol, prec), pi=mppi(prec);//Difference in angles in [0, 2*Pi]
  int topi=tolcmp(ang, pi, tol, prec);//We want to be in the range [0, Pi), so potentially subtract off Pi
  if(topi==-1) return gerepileupto(top, ang);
  else if(topi==0){avma=top;return gen_0;}//Same angle
  return gerepileupto(top, gsub(ang, pi));//Subtract off pi
}

//Circle with centre cent passing through p
static GEN circle_fromcp(GEN cent, GEN p, long prec){
  pari_sp top=avma;
  GEN pmcent=gsub(p, cent);
  GEN ret=cgetg(4, t_VEC);
  gel(ret, 1)=gcopy(cent);
  gel(ret, 2)=gabs(pmcent, prec);
  gel(ret, 3)=gen_0;
  return gerepileupto(top, ret);
}

//Circle through 3 points (with one allowed to be oo, making a line instead)
static GEN circle_fromppp(GEN p1, GEN p2, GEN p3, GEN tol, long prec){
  if(typ(p1)==t_INFINITY) return line_frompp(p2, p3, tol, prec);
  if(typ(p2)==t_INFINITY) return line_frompp(p1, p3, tol, prec);
  if(typ(p3)==t_INFINITY) return line_frompp(p1, p2, tol, prec);//Lines
  pari_sp top=avma;
  GEN l1=perpbis(p1, p2, tol, prec), l2=perpbis(p1, p3, tol, prec);
  GEN centre=line_int(l1, l2, tol, prec);//centre is intersection of perp bisectors.
  if(typ(centre)==t_INFINITY) return gerepileupto(top, line_frompp(p1, p2, tol, prec));//p1, p2, p3 collinear
  return gerepileupto(top, circle_fromcp(centre, p1, prec));//The circle!
}

//Returns the slope of the tangent to c at p
static GEN circle_tangentslope(GEN c, GEN p, long prec){
  pari_sp top=avma;
  GEN c1mp=gsub(gel(c, 1), p);//c[1]-p
  GEN c1mpr=real_i(c1mp);
  GEN c1mpi=imag_i(c1mp);
  return gerepileupto(top, divoo(c1mpr, gneg(c1mpi)));//divoo(real(c[1])-real(p),imag(p)-imag(c[1])));
}

//The line through p with slope s
static GEN line_fromsp(GEN s, GEN p){
  if(typ(s)==t_INFINITY){//oo slope
    GEN ret=cgetg(4, t_VEC);
    gel(ret, 1)=mkoo();
    gel(ret, 2)=greal(p);//x-intercept
    gel(ret, 3)=gen_1;
    return ret;
  }
  pari_sp top=avma;
  GEN strealp=gmul(s, real_i(p)), imagp=imag_i(p);
  GEN ret=cgetg(4, t_VEC);
  gel(ret, 1)=gcopy(s);
  gel(ret, 2)=gsub(imagp, strealp);//y=intercept
  gel(ret, 3)=gen_1;
  return gerepileupto(top, ret);
}

//Line through two points
static GEN line_frompp(GEN p1, GEN p2, GEN tol, long prec){
  pari_sp top=avma;
  return gerepileupto(top, line_fromsp(slope(p1, p2, tol, prec), p1));
}

//Mx, where M is a 2x2 matrix and x is complex or infinite.
GEN mat_eval(GEN M, GEN x){
  pari_sp top = avma;
  if(typ(x)==t_INFINITY) return gerepileupto(top,divoo(gcoeff(M, 1, 1), gcoeff(M, 2, 1)));
  return gerepileupto(top,divoo(gadd(gmul(gcoeff(M, 1, 1), x), gcoeff(M, 1, 2)), gadd(gmul(gcoeff(M, 2, 1), x), gcoeff(M, 2, 2))));
}

//Midpoint of p1 and p2
static GEN midpoint(GEN p1, GEN p2){
  pari_sp top=avma;
  return gerepileupto(top, gdivgs(gadd(p1, p2), 2));
}

//Mobius transform accessible in GP
GEN mobius_gp(GEN M, GEN c, long prec){
  pari_sp top=avma;
  GEN tol=deftol(prec);
  return gerepileupto(top, mobius(M, c, tol, prec));
}

//Returns M(c), for c a circle/line/arc/segment
static GEN mobius(GEN M, GEN c, GEN tol, long prec){
  int i=geom_check(c);
  switch(i){
    case 0: return mobius_circle(M, c, tol, prec);
    case 1: return mobius_line(M, c, tol, prec);
    case 2: return mobius_arcseg(M, c, 1, tol, prec);
    case 3: return mobius_arcseg(M, c, 0, tol, prec);
  }
  pari_err_TYPE("Please input a circle/line/arc/segment", c);
  return gen_0;//Never reach this far
}

//Mobius map acting on an arc or segment
static GEN mobius_arcseg(GEN M, GEN c, int isarc, GEN tol, long prec){
  pari_sp top=avma;//We start by finding the new circle/line(ignoring the arc/segment bits)
  GEN endpt1, endpt2, extpt;//Endpoints and an extra point (used to have 3 points to translate)
  endpt1=mat_eval(M, gel(c, 3));
  endpt2=mat_eval(M, gel(c, 4));
  //Now we must find an extra point extpt
  if(isarc==1) extpt=mat_eval(M, arc_midpoint(c, gel(c, 3), gel(c, 4), tol, prec));//arc
  else{//segment
    if(typ(gel(c, 3))==t_INFINITY){//Start point infinity
      GEN u;
      if(gequal(gel(c, 6), gen_1)) u=gen_m2;//segment goes vertically up or right
      else u=gen_2;//segment goes vertically down or left
      if(typ(gel(c, 1))==t_INFINITY) extpt=mat_eval(M, gadd(gel(c, 4), gmul(u, gen_I())));//Vertical line, M(c[4]+U*I)
      else extpt=mat_eval(M, gadd(gel(c, 4), gmul(u, gaddsg(1, gmul(gel(c, 1), gen_I())))));//non-vertical line, M(c[4]+u+u*c[1]*I)
    }
    else if(typ(gel(c, 4))==t_INFINITY){//End point infinity
      GEN u;
      if(gequal(gel(c, 6), gen_1)) u=gen_2;//segment goes vertically up or right
      else u=gen_m2;//segment goes vertically down or left
      if(typ(gel(c, 1))==t_INFINITY) extpt=mat_eval(M, gadd(gel(c, 3), gmul(u, gen_I())));//Vertical line, M(c[3]+U*I)
      else extpt=mat_eval(M, gadd(gel(c, 3), gmul(u, gaddsg(1, gmul(gel(c, 1), gen_I())))));//non-vertical line, M(c[3]+u+u*c[1]*I)
    }
    else{//Start/end points in the plane
      if(gequal(gel(c, 7), gen_1)) extpt=mat_eval(M, midpoint(gel(c, 3), gel(c, 4)));//Does not go through oo, can take midpoint
      else extpt=mat_eval(M, mkoo());//Use oo, since the line goes through there
    }
  }
  //Now we have the 3 new points used to define the new arc/segment. Let's finish the process.
  GEN ret=cgetg(ARCLEN, t_VEC);//The returned arc/segment
  GEN newcirc=circle_fromppp(endpt1, endpt2, extpt, tol, prec);//The new circle/line
  gel(ret, 1)=gel(newcirc, 1);//Slope/centre
  gel(ret, 2)=gel(newcirc, 2);//x or y intercept/radius
  gel(ret, 3)=endpt1;//Start point
  gel(ret, 4)=endpt2;//end point. These may be in the wrong order, so will fix later if so.
  if(gequal0(gel(newcirc, 3))){//Circle
    gel(ret, 8)=gen_0;//arc
    gel(ret, 5)=radialangle(newcirc, endpt1, tol, prec);//angle 1
    gel(ret, 6)=shiftangle(radialangle(newcirc, endpt2, tol, prec), gel(ret, 5), tol, prec);//angle 2
    if(isarc) gel(ret, 7)=gel(c, 7);//Temporary; match the old to the new direction
    else gel(ret, 7)=gen_1;//Temporary; since for segments we have a bona fide start and end, we start with the direction being 1.
    if(!onarc(ret, extpt, tol, prec)){//Must swap start/endpoints, angles and the direction
      gel(ret, 3)=endpt2;
      gel(ret, 4)=endpt1;
      GEN tempang=gel(ret, 5);//angle to endpt1
      gel(ret, 5)=shiftangle(gel(ret, 6), gen_0, tol, prec);//The angle to endpt2 shifted to [0, 2*Pi)
      gel(ret, 6)=shiftangle(tempang, gel(ret, 5), tol, prec);//Angle to endpt1 shifted with a base of the angle to endpt2
      gel(ret, 7)=gneg(gel(ret, 7));//We now run backwards
    }
  }
  else{//Line
    if(isarc==1 && gequal(gel(c, 7), gen_m1)){//We need to reverse the order of the points, because the arc ran backwards.
      gel(ret, 3)=endpt2;
      gel(ret, 4)=endpt1;
    }
    gel(ret, 8)=gen_1; //segment
    gel(ret, 5)=gen_0;//Unused
    if(typ(endpt1)==t_INFINITY || typ(endpt2)==t_INFINITY){//oo endpoint
      gel(ret, 6)=gen_1;//Temporary, assume we go up/right
      gel(ret, 7)=gen_0;//Not used
      if(!onseg(ret, extpt, tol, prec)) gel(ret, 6)=gen_m1;//We were wrong, and go down/left.
    }
    else{//Both endpoints finite
      gel(ret, 6)=gen_0;
      gel(ret, 7)=gen_1;//Temporary, assume we go through the plane only
      if(!onseg(ret, extpt, tol, prec)) gel(ret, 7)=gen_m1;//We were wrong, and go through oo.
    }
  }
  return gerepilecopy(top, ret);
}

//Mobius map acting on circle
static GEN mobius_circle(GEN M, GEN c, GEN tol, long prec){
  pari_sp top=avma;
  GEN p1=mat_eval(M, gadd(gel(c, 1), gel(c, 2)));//M(c[1]+c[2])
  GEN p2=mat_eval(M, gadd(gel(c, 1), gmul(gel(c, 2), gen_I())));//M(c[1]+c[2]*I)
  GEN p3=mat_eval(M, gsub(gel(c, 1), gel(c, 2)));//M(c[1]-c[2])
  return gerepileupto(top, circle_fromppp(p1, p2, p3, tol, prec));
}

//Mobius map acting on line
static GEN mobius_line(GEN M, GEN l, GEN tol, long prec){
  pari_sp top=avma;
  GEN p1, p2, p3, I=gen_I();
  if(typ(gel(l, 1))==t_INFINITY){//Vertical line
    p1=mat_eval(M, gel(l, 2));//M(x-intercept)
    p2=mat_eval(M, gadd(gel(l, 2), I));//M(x-intercept+I)
    p3=mat_eval(M, gsub(gel(l, 2), I));//M(x-intercept-I)
  }
  else{//Non-vertical line
    GEN slopeIp1=gaddgs(gmul(gel(l, 1), I), 1);//1+Slope*I
    GEN p1base=gmul(gel(l, 2), I);//y-intercept
    GEN p2base=gadd(p1base, slopeIp1);//y-intercept+1+slope*I
    GEN p3base=gadd(p2base, slopeIp1);//y-intercept+2+2*slope*I
    p1=mat_eval(M, p1base);p2=mat_eval(M, p2base);p3=mat_eval(M, p3base);
  }
  return gerepileupto(top, circle_fromppp(p1, p2, p3, tol, prec));
}

//Perpendicular bisector of distinct points
static GEN perpbis(GEN p1, GEN p2, GEN tol, long prec){
  pari_sp top=avma;
  return gerepileupto(top, line_fromsp(divoo(gen_m1, slope(p1, p2, tol, prec)), midpoint(p1, p2)));
}

//Angle between p and the centre of c, in the range [0, 2*Pi)
static GEN radialangle(GEN c, GEN p, GEN tol, long prec){
  pari_sp top=avma;
  return gerepileupto(top, shiftangle(garg(gsub(p, gel(c, 1)), prec), gen_0, tol, prec));
}

//The slope of the line through p1, p2
static GEN slope(GEN p1, GEN p2, GEN tol, long prec){
  pari_sp top=avma;
  GEN p2mp1=gsub(p2, p1);
  GEN ftop=imag_i(p2mp1);
  GEN fbot=real_i(p2mp1);
  if(toleq(fbot, gen_0, tol, prec)) fbot=gen_0;
  if(toleq(ftop, gen_0, tol, prec)) ftop=gen_0;
  return gerepileupto(top, divoo(ftop, fbot));
}


//INTERSECTION OF LINES/CIRCLES


//Returns the intersection points of two arcs
static GEN arc_int(GEN c1, GEN c2, GEN tol, long prec){
  pari_sp top=avma;
  GEN ipts=circle_int(c1, c2, tol, prec);
  if(lg(ipts)==1){avma=top;return cgetg(1, t_VEC);}//No intersection
  if(lg(ipts)==2){//One intersection point (tangent circles)
    if(!onarc(c1, gel(ipts, 1), tol, prec)){avma=top;return cgetg(1, t_VEC);}//Not on arc 1
    if(!onarc(c2, gel(ipts, 1), tol, prec)){avma=top;return cgetg(1, t_VEC);}//Not on arc 2
    return gerepilecopy(top, ipts);//On arc
  }
  //Two intersections
  int i1=onarc(c1, gel(ipts, 1), tol, prec);
  if(i1==1) i1=onarc(c2, gel(ipts, 1), tol, prec);//Now i1==1 iff the ipts[1] is on both c1 and c2
  int i2=onarc(c1, gel(ipts, 2), tol, prec);
  if(i2==1) i2=onarc(c2, gel(ipts, 2), tol, prec);//Now i2==1 iff the ipts[2] is on both c1 and c2
  if(i1==1){
    if(i2==1) return gerepilecopy(top, ipts);//Both pts on the arcs
    GEN ret=cgetg(2, t_VEC);//Just point 1
    gel(ret, 1)=gcopy(gel(ipts, 1));
    return gerepileupto(top, ret);
  }
  //Now i1=0
  if(i2==0){avma=top;return cgetg(1, t_VEC);}//Not on either arc
  GEN ret=cgetg(2, t_VEC);//Just point 2
  gel(ret, 1)=gcopy(gel(ipts, 2));
  return gerepileupto(top, ret);
}

//Returns the intersection points of an arc and a segment
static GEN arcseg_int(GEN c, GEN l, GEN tol, long prec){
  pari_sp top=avma;
  GEN ipts=circleline_int(c, l, tol, prec);
  if(lg(ipts)==1){avma=top;return cgetg(1, t_VEC);}//No intersection
  if(lg(ipts)==2){//One intersection point (tangent circles)
    if(!onarc(c, gel(ipts, 1), tol, prec)){avma=top;return cgetg(1, t_VEC);}//Not on arc
    if(!onseg(l, gel(ipts, 1), tol, prec)){avma=top;return cgetg(1, t_VEC);}//Not on segment
    return gerepilecopy(top, ipts);//On both
  }
  //Two intersections
  int i1=onarc(c, gel(ipts, 1), tol, prec);
  if(i1==1) i1=onseg(l, gel(ipts, 1), tol, prec);//Now i1==1 iff the ipts[1] is on both c and l
  int i2=onarc(c, gel(ipts, 2), tol, prec);
  if(i2==1) i2=onseg(l, gel(ipts, 2), tol, prec);//Now i2==1 iff the ipts[2] is on both c and l
  if(i1==1){
    if(i2==1) return gerepilecopy(top, ipts);//Both pts on both
    GEN ret=cgetg(2, t_VEC);//Just point 1
    gel(ret, 1)=gcopy(gel(ipts, 1));
    return gerepileupto(top, ret);
  }
  //Now i1=0
  if(i2==0){avma=top;return cgetg(1, t_VEC);}//Not on either
  GEN ret=cgetg(2, t_VEC);//Just point 2
  gel(ret, 1)=gcopy(gel(ipts, 2));
  return gerepileupto(top, ret);
}

//Returns the set of points in the intersection of circles c1, c2
static GEN circle_int(GEN c1, GEN c2, GEN tol, long prec){
  pari_sp top=avma;
  GEN a1=real_i(gel(c1, 1)), b1=imag_i(gel(c1, 1)), r1=gel(c1, 2);//x, y coords and radius of c1
  GEN a2=real_i(gel(c2, 1)), b2=imag_i(gel(c2, 1)), r2=gel(c2, 2);//x, y coords and radius of c2
  GEN a1ma2=gsub(a1, a2), b1mb2=gsub(b1, b2), x1, x2, y1, y2;
  int oneint=0;
  if(gcmp(gabs(a1ma2, prec), gabs(b1mb2, prec))>=0){//We want to divide by the larger of the two quantities to maximize precision and avoid errors when the centres are on the same line.
    if(toleq(a1ma2, gen_0, tol, prec)==1){avma=top;return cgetg(1, t_VEC);}//Same centre, cannot intersect.
    //u=(r1^2-r2^2+b2^2-b1^2+a2^2-a1^2)/(2*a2-2*a1)-a1;
    GEN u=gsub(gdiv(gsub(gadd(gsub(gadd(gsub(gsqr(r1), gsqr(r2)), gsqr(b2)), gsqr(b1)), gsqr(a2)), gsqr(a1)), gmulgs(a1ma2, -2)), a1);
    GEN v=gneg(gdiv(b1mb2, a1ma2));//v=(b1-b2)/(a2-a1), and x=a1+u+vy
    GEN uvmb1=gsub(gmul(u, v), b1);//uv-b1
    GEN vsqrp1=gaddgs(gsqr(v), 1);//v^2+1
    GEN rtpart=gsub(gsqr(uvmb1), gmul(vsqrp1, gadd(gsqr(b1), gsub(gsqr(u), gsqr(r1)))));//(u*v-b1)^2-(v^2+1)*(b1^2+u^2-r1^2)
    oneint=tolcmp(rtpart, gen_0, tol, prec);//Comparing rtpart to 0
    if(oneint==-1){avma=top;return cgetg(1, t_VEC);}//rtpart must be square rooted, so if it's negative the circles do not intersect
    if(oneint==0){//One intersection, so we take rtpart=0. This is CRUCIAL, as taking the square root kills our precision if we don't do this here.
      y1=gdiv(gneg(uvmb1), vsqrp1);//y1=(b1-u*v)/(1*v^2+1)
      x1=gadd(gadd(a1, u), gmul(v, y1));//x1=a1+u+v*y1
    }
    else{
      y1=gdiv(gadd(gneg(uvmb1), gsqrt(rtpart, prec)), vsqrp1);//y1=(b1-u*v+sqrt((u*v-b1)^2-*(v^2+1)*(b1^2+u^2-r1^2)))/(1*v^2+1)
      y2=gadd(gneg(y1), gdiv(gmulgs(uvmb1, -2), vsqrp1));//y2=-y1+(2*b1-2*u*v)/(v^2+1)
      GEN a1pu=gadd(a1, u);
      x1=gadd(a1pu, gmul(v, y1));//x1=a1+u+v*y1
      x2=gadd(a1pu, gmul(v, y2));//x1=a1+u+v*y2
    }
  }
  else{
    if(toleq(b1mb2, gen_0, tol, prec)==1){avma=top;return cgetg(1, t_VEC);}//Same centre, cannot intersect.
    //u=(r1^2-r2^2+b2^2-b1^2+a2^2-a1^2)/(2*b2-2*b1)-b1;
    GEN u=gsub(gdiv(gsub(gadd(gsub(gadd(gsub(gsqr(r1), gsqr(r2)), gsqr(b2)), gsqr(b1)), gsqr(a2)), gsqr(a1)), gmulgs(b1mb2, -2)), b1);
    GEN v=gneg(gdiv(a1ma2, b1mb2));//v=(a1-a2)/(b2-b1), and y=b1+u+vx
    GEN uvma1=gsub(gmul(u, v), a1);//uv-a1
    GEN vsqrp1=gaddgs(gsqr(v), 1);//v^2+1
    GEN rtpart=gsub(gsqr(uvma1), gmul(vsqrp1, gadd(gsqr(a1), gsub(gsqr(u), gsqr(r1)))));//(u*v-a1)^2-(v^2+1)*(a1^2+u^2-r1^2))
    oneint=tolcmp(rtpart, gen_0, tol, prec);//Comparing rtpart to 0
    if(oneint==-1){avma=top;return cgetg(1, t_VEC);}//rtpart must be square rooted, so if it's negative the circles do not intersect
    if(oneint==0){//One intersection, so we take rtpart=0. This is CRUCIAL, as taking the square root kills our precision if we don't do this here.
      x1=gdiv(gneg(uvma1), vsqrp1);//x1=(a1-u*v)/(v^2+1);
      y1=gadd(gadd(b1, u), gmul(v, x1));//y1=b1+u+v*x1;
    }
    else{
      x1=gdiv(gadd(gneg(uvma1), gsqrt(rtpart, prec)), vsqrp1);//x1=(a1-u*v+sqrt((u*v-a1)^2-(v^2+1)*(a1^2+u^2-r1^2)))/(v^2+1);
      x2=gadd(gneg(x1), gdiv(gmulgs(uvma1, -2), vsqrp1));//x2=-x1+(2*a1-2*u*v)/(v^2+1)
      GEN b1pu=gadd(b1, u);
      y1=gadd(b1pu, gmul(v, x1));//y1=b1+u+v*x1;
      y2=gadd(b1pu, gmul(v, x2));//y2=b1+u+v*x2;
    }
  }
  if(oneint==0){//One point of intersection (0 pts of intersection was already dealt with and returned)
    GEN y1I=gmul(gen_I(), y1);
    GEN ret=cgetg(2, t_VEC);
    gel(ret, 1)=gadd(x1, y1I);
    return gerepileupto(top, ret);
  }
  GEN y1I=gmul(gen_I(), y1), y2I=gmul(gen_I(), y2);
  GEN ret=cgetg(3, t_VEC);
  gel(ret, 1)=gadd(x1, y1I);
  gel(ret, 2)=gadd(x2, y2I);
  return gerepileupto(top, ret);
}

//Returns the intersection points of c and l
static GEN circleline_int(GEN c, GEN l, GEN tol, long prec){
  pari_sp top=avma;
  if(typ(gel(l, 1))==t_INFINITY){
    GEN x1=gel(l, 2);
    GEN rtpart=gsub(gsqr(gel(c, 2)), gsqr(gsub(x1, real_i(gel(c, 1)))));//c[2]^2-(x1-real(c[1]))^2
    if(gsigne(rtpart)==-1){avma=top;return cgetg(1, t_VEC);}//No intersections.
    GEN y1=gadd(imag_i(gel(c, 1)), gsqrt(rtpart, prec));//y1=imag(c[1])+sqrt(c[2]^2-(x1-real(c[1]))^2)
    if(toleq(rtpart, gen_0, tol, prec)){//Only one intersection point
      GEN ret=cgetg(2, t_VEC);
      gel(ret, 1)=cgetg(3, t_COMPLEX);
      gel(gel(ret, 1), 1)=gcopy(x1);
      gel(gel(ret, 1), 2)=gcopy(y1);
      return gerepileupto(top, ret);
    }
    //Two intersection points
    GEN y1py2=gmulgs(imag_i(gel(c, 1)), 2);//2*imag(c[1])
    GEN ret=cgetg(3, t_VEC);
    gel(ret, 1)=cgetg(3, t_COMPLEX);gel(ret, 2)=cgetg(3, t_COMPLEX);
    gel(gel(ret, 1), 1)=gcopy(x1);
    gel(gel(ret, 1), 2)=gcopy(y1);
    gel(gel(ret, 2), 1)=gcopy(x1);
    gel(gel(ret, 2), 2)=gsub(y1py2, y1);
    return gerepileupto(top, ret);
  }
  //Now y=mx+b with m finite
  GEN A=gaddgs(gsqr(gel(l, 1)), 1);//l[1]^2+1
  GEN l2mic1=gsub(gel(l, 2), imag_i(gel(c, 1)));//l[2]-imag(c[1])
  GEN B=gadd(gmulgs(real_i(gel(c, 1)), -2), gmulsg(2, gmul(gel(l, 1), l2mic1)));//-2*real(c[1])+2*l[1]*(l[2]-imag(c[1]))
  GEN C=gadd(gsqr(real_i(gel(c, 1))), gsub(gsqr(l2mic1), gsqr(gel(c, 2))));//real(c[1])^2+(l[2]-imag(c[1]))^2-c[2]^2
  GEN rtpart=gsub(gsqr(B), gmulsg(4, gmul(A, C)));
  int rtpartsig=tolcmp(rtpart, gen_0, tol, prec);
  if(rtpartsig==-1){avma=top;return cgetg(1, t_VEC);}//No intersection
  if(rtpartsig==0){//One root, and rtpart=0
    GEN x1=gdiv(B, gmulgs(A, -2));//-B/(2A)
    GEN y1part=gmul(gel(l, 1), x1);//l[1]*x1
    GEN ret=cgetg(2, t_VEC);
    gel(ret, 1)=cgetg(3, t_COMPLEX);
    gel(gel(ret, 1), 1)=gcopy(x1);
    gel(gel(ret, 1), 2)=gadd(y1part, gel(l, 2));//y1=l[1]*x1+l[2];
    return gerepileupto(top, ret);
  }
  //Two roots
  GEN x1=gdiv(gsub(gsqrt(rtpart, prec), B), gmulgs(A, 2));//x1=(-B+sqrt(B^2-4*A*C))/(2*A);
  GEN x2=gsub(gneg(gdiv(B, A)), x1);//-B/A-x1
  GEN y1part=gmul(gel(l, 1), x1);//l[1]*x1
  GEN y2part=gmul(gel(l, 1), x2);//l[1]*x2
  GEN ret=cgetg(3, t_VEC);
  gel(ret, 1)=cgetg(3, t_COMPLEX);gel(ret, 2)=cgetg(3, t_COMPLEX);
  gel(gel(ret, 1), 1)=gcopy(x1);
  gel(gel(ret, 1), 2)=gadd(y1part, gel(l, 2));//l[1]*x1+l[2];
  gel(gel(ret, 2), 1)=gcopy(x2);
  gel(gel(ret, 2), 2)=gadd(y2part, gel(l, 2));//l[1]*x2+l[2];
  return gerepileupto(top, ret);
}

//The intersection of two lines
static GEN line_int(GEN l1, GEN l2, GEN tol, long prec){
  GEN s1=gel(l1, 1), s2=gel(l2, 1);//Slopes
  if(toleq(s1, s2, tol, prec)) return mkoo();//Parallel or equal
  pari_sp top=avma;
  if(typ(s1)==t_INFINITY){//l1 vertical
    GEN ypart=gmul(s2, gel(l1, 2));//s2*l1[2]
    GEN ipt=cgetg(3, t_COMPLEX);
    gel(ipt, 1)=gcopy(gel(l1, 2));
    gel(ipt, 2)=gadd(ypart, gel(l2, 2));//s2*l1[2]+l2[2]
    return gerepileupto(top, ipt);
  }
  if(typ(s2)==t_INFINITY){//l2 vertical
    GEN ypart=gmul(s1, gel(l2, 2));//s1*l2[2]
    GEN ipt=cgetg(3, t_COMPLEX);
    gel(ipt, 1)=gcopy(gel(l2, 2));
    gel(ipt, 2)=gadd(ypart, gel(l1, 2));//s1*l2[2]+l1[2]
    return gerepileupto(top, ipt);
  }
  GEN x=gdiv(gsub(gel(l2, 2), gel(l1, 2)), gsub(s1, s2));//(l2[2]-l1[2])/(s1-s2)
  GEN ypart=gmul(s1, x);//s1*x
  GEN ipt=cgetg(3, t_COMPLEX);
  gel(ipt, 1)=gcopy(x);
  gel(ipt, 2)=gadd(ypart, gel(l1, 2));//s1*x+l1[2]
  return gerepileupto(top, ipt);
}

//p is assumed to be on the circle defined by c; this checks if it is actually on the arc (running counterclockwise from c[3] to c[4]).
static int onarc(GEN c, GEN p, GEN tol, long prec){
  if(lg(c)==CIRCLEN) return 1;//Allow input of just a circle, so the return is trivially 1
  if(toleq(gel(c, 3), p, tol, prec)) return 1;//p=the start point. We have this done seperately in case rounding errors take the angle to <c[5], as this will cause issues with the shifting angle.
  pari_sp top=avma;
  GEN angle=shiftangle(radialangle(c, p, tol, prec), gel(c, 5), tol, prec);//Getting the angle in the range [c[5], c[5]+2*Pi)
  if(tolcmp(angle, gel(c, 6), tol, prec)<=0){avma=top;return 1;}//On the arc
  avma=top;
  return 0;//Beyond the arc.
}

//p is assumed to be on the line defined by l; this checks if it is actually on the segment l
static int onseg(GEN l, GEN p, GEN tol, long prec){
  if(lg(l)==CIRCLEN) return 1;//Allow input of a line, so return is trivially 1
  if(typ(p)==t_INFINITY){//p is the point at oo
    if(!gequal0(gel(l, 6)) || gequal(gel(l, 7), gen_m1)) return 1;//oo is an endpoint OR the seg passes through oo
    return 0;//If not, does not pass through oo
  }
  pari_sp top=avma;
  //Okay, now p is not oo and l is a line segment
  if(typ(gel(l, 1))==t_INFINITY){//Vertical line
    if(typ(gel(l, 3))==t_INFINITY){//Start point in oo
      if(equali1(gel(l, 6))){//p must lie BELOW l[4]
          if(tolcmp(imag_i(p), imag_i(gel(l, 4)), tol, prec)<=0){avma=top;return 1;}//Lies below l[4]
          avma=top;return 0;//lies above l[4]
      }
      //p must lie ABOVE l[4]
      if(tolcmp(imag_i(p), imag_i(gel(l, 4)), tol, prec)>=0){avma=top;return 1;}//Lies above l[4]
      avma=top;return 0;//lies below l[4]
    }
    if(typ(gel(l, 4))==t_INFINITY){//End point is oo
      if(equali1(gel(l, 6))){//p must lie ABOVE l[3]
          if(tolcmp(imag_i(p), imag_i(gel(l, 3)), tol, prec)>=0){avma=top;return 1;}//Lies above l[3]
          avma=top;return 0;//lies below l[3]
      }
      //p must lie BELOW l[3]
      if(tolcmp(imag_i(p), imag_i(gel(l, 3)), tol, prec)<=0){avma=top;return 1;}//Lies below l[3]
      avma=top;return 0;//lies above l[3]
    }
    //Start and end points are finite
    int i1=tolcmp(imag_i(gsub(p, gel(l, 3))), gen_0, tol, prec);//sign of imag(p)-imag(l[3])
    int i2=tolcmp(imag_i(gsub(p, gel(l, 4))), gen_0, tol, prec);//sign of imag(p)-imag(l[4])
    avma=top;
    if(i1==0 || i2==0) return 1;//endpoint
    if(i1==i2){//p on the same side of l[3] and l[4], so return 1 iff l passes through oo
        if(gequal(gel(l, 7), gen_1)) return 0;//Not through oo
        return 1;//through oo
    }
    //p is between l[3] and l[4], so return 1 iff l does not pass through oo
    if(gequal(gel(l, 7), gen_1)) return 1;//not through oo
    return 0;//through oo
  }
  //Non-vertical line
  if(typ(gel(l, 3))==t_INFINITY){//Start point in oo
    if(equali1(gel(l, 6))){//p must lie LEFT OF l[4]
      if(tolcmp(real_i(p), real_i(gel(l, 4)), tol, prec)<=0){avma=top;return 1;}//Lies left of l[4]
      avma=top;return 0;//lies right of  l[4]
    }
    //p must lie RIGHT OF l[4]
    if(tolcmp(real_i(p), real_i(gel(l, 4)), tol, prec)>=0){avma=top;return 1;}//Lies right of l[4]
    avma=top;return 0;//lies left of l[4]
  }
  if(typ(gel(l, 4))==t_INFINITY){//End point is oo
    if(equali1(gel(l, 6))){//p must lie RIGHT OF l[3]
      if(tolcmp(real_i(p), real_i(gel(l, 3)), tol, prec)>=0){avma=top;return 1;}//Lies right of l[3]
      avma=top;return 0;//lies below l[3]
    }
    //p must lie LEFT OF l[3]
    if(tolcmp(real_i(p), real_i(gel(l, 3)), tol, prec)<=0){avma=top;return 1;}//Lies left of l[3]
    avma=top;return 0;//lies above l[3]
  }
  //Start and end points are finite
  int i1=tolcmp(real_i(gsub(p, gel(l, 3))), gen_0, tol, prec);//sign of real(p)-real(l[3])
  int i2=tolcmp(real_i(gsub(p, gel(l, 4))), gen_0, tol, prec);//sign of real(p)-real(l[4])
  avma=top;
  if(i1==0 || i2==0) return 1;//endpoint
  if(i1==i2){//p on the same side of l[3] and l[4], so return 1 iff l passes through oo
    if(gequal(gel(l, 7), gen_1)) return 0;//Not through oo
    return 1;//through oo
  }
  //p is between l[3] and l[4], so return 1 iff l does not pass through oo
  if(gequal(gel(l, 7), gen_1)) return 1;//not through oo
  return 0;//through oo
}


//DISTANCES/AREAS


//Given a radius R>0, this returns the area of the hyperbolic disc of radius R. The formula is 4*Pi*sinh(R/2)^2
GEN hdiscarea(GEN R, long prec){
  pari_sp top=avma;
  return gerepileupto(top, gtofp(gmul(Pi2n(2, prec), gsqr(gsinh(gdivgs(R, 2), prec))), prec));
}

//Given the area of a hyperbolic disc, this returns the radius.
GEN hdiscradius(GEN area, long prec){
  pari_sp top=avma;
  return gerepileupto(top, gtofp(gmulgs(gasinh(gsqrt(gdiv(area, Pi2n(2, prec)), prec), prec), 2), prec));
}

//z1 and z2 are complex numbers, this computes the hyperbolic distance between them.
GEN hdist(GEN z1, GEN z2, long prec){
  pari_sp top=avma;
  GEN expd;
  pari_CATCH(CATCH_ALL){
	avma=top;
	pari_CATCH_reset();
	pari_err_TYPE("Please enter two complex numbers in the upper half plane", mkvec2(z1, z2));
	return gen_0;
  }
  pari_TRY{
    GEN x1=gel(z1,1);
    GEN y1=gel(z1,2);
    GEN x2=gel(z2,1);
    GEN y2=gel(z2,2);
    GEN x=gaddsg(1,gdiv(gadd(gsqr(gsub(x2,x1)),gsqr(gsub(y2,y1))),gmul(gmulsg(2,y1),y2)));
    expd=gadd(x,gsqrt(gsubgs(gsqr(x), 1), prec));
  }
  pari_ENDCATCH
  return gerepileupto(top,glog(expd, prec));
}

//The hyperbolic distance between z1 and z2 in the unit disc model
GEN hdist_ud(GEN z1, GEN z2, long prec){
  pari_sp top=avma;
  GEN a = gabs(gsubsg(1, gmul(z1, conj_i(z2))), prec);//|1-z1*conj(z2)|
  GEN b = gabs(gsub(z1, z2), prec);//|z1-z2|
  GEN num=gadd(a, b);
  GEN denom=gsub(a, b);
  GEN ret;
  pari_CATCH(e_INV){
    avma=top;
	pari_CATCH_reset();
    return mkoo();
  }
  pari_TRY{
    ret=gerepileupto(top, glog(gdiv(num, denom), prec));//log((a+b)/(a-b))
  }
  pari_ENDCATCH
  return ret;
}

//Given the isometric circles (in order) and the vertices (so that circles[i] intersect circles[i+1] is vertices[i]), returns the area of the convex hyperbolic polygon. If one of the sides is oo (an entry of circles is 0), the answer will be oo. If there are n vertices with angles a1,...,an, then the area is is (n-2)Pi-sum(ai)
static GEN hpolygon_area(GEN circles, GEN vertices, GEN tol, long prec){
  pari_sp top=avma;
  long blen=lg(circles);
  if(blen==1 || gequal0(gel(circles, 1))) return mkoo();//No cicles or the first is 0, i.e. an infinite side
  GEN ang, area=gmulsg(blen-3, mppi(prec));//We subtract off from area.
  for(long i=1;i<blen-1;i++){
    if(gequal0(gel(circles, i+1))){avma=top;return mkoo();}//The next side is infinite
    ang=circle_angle(gel(circles, i+1), gel(circles, i), gel(vertices, i), tol, prec);//Do in opposite order to get correct angle
    area=gsub(area, ang);
  }
  ang=circle_angle(gel(circles, 1), gel(circles, blen-1), gel(vertices, blen-1), tol, prec);//The last side wraps around.
  area=gsub(area, ang);
  return gerepileupto(top, area);
}


//FUNDAMENTAL DOMAIN COMPUTATION


/*Let Gamma be a Fuschian group; the following methods allow us to compute a fundamental domain for Gamma, following the paper of John Voight. We are assuming that we have an "exact" way to deal with the elements of Gamma. In otherwords, we need separate methods to:
i)  Multiply elements of Gamma: format as GEN eltmul(GEN *data, GEN x, GEN y), with data the extra data you need.
ii) Find the image of g in PSL(1, 1): format as GEN gamtopsl(GEN *data, GEN g, long prec), with data the extra data you need (e.g. the quaterniona algebra in the case of Shimura curves).
iii) Identify if g is trivial in Gamma: format as int istriv(GEN *data, GEN g). Need to take care as we care about being trivial in PSL, whereas most representations of Gamma are for SL (so need to check with +/-id).
iv) Pass in the identity element of Gamma and find the area of the fundamental domain. These methods are not passed in; just the values.
We pass references to these methods into the methods here.
We do computations mostly in PSU, and shift from PSL to PSU via phi(z):=(z-p)/(z-conj(p)) for some given upper half plane point p.
*/

//Some of the methods may not be that useful (i.e. assuming we have p only and not the matrices).

//Returns the edge pairing, as VECSMALL v where v[i]=j means i is paired with j. If not all edges can be paired, instead returns [v1, v2, ...] where vi is a VECSMALL that is either [gind, v1ind] or [gind, v1ind, v2ind]. gind is the index of the unpaired side, and the viind are the corresponding unpaired vertices (1 or 2). If rboth=1, returns [paired, unpaired], where this time paired is a vector of vecsmalls [i, j] meaning i is paired to j (differing to the output format when rboth=0 and everything is paired)
static GEN edgepairing(GEN U, GEN tol, int rboth, long prec){
  pari_sp top=avma;
  GEN vangles=gel(U, 4);//Vertex angles
  GEN baseangle=gel(vangles, 1);//Base angle
  GEN toldata=cgetg(3, t_VEC);//Stores the necessary info for searching with tolerance (tolcmp_sort)
  gel(toldata, 1)=tol;
  gel(toldata, 2)=cgetg(2, t_VECSMALL);
  gel(toldata, 2)[1]=prec;
  long lU=lg(gel(U, 1));
  GEN unpair=vectrunc_init(lU+1), vim, vimang, pair=vectrunc_init(lU);//Unpair stores the unpaired edges, pair stores the paired edges
  long ind1, ind2, i1, i2;
  for(long i=1;i<lU;i++){
    if(gequal0(gel(gel(U, 5), i))){vectrunc_append(pair, mkvecsmall2(i, i));continue;}//oo side, go next (we say it is paired with itself)
    ind1=i;
    vim=mat_eval(gel(gel(U, 5), i), gel(gel(U, 3), ind1));//The new vertex
    vimang=shiftangle(garg(vim, prec), baseangle, tol, prec);//The new angle
    i1=gen_search(vangles, vimang, 0, &toldata, &tolcmp_sort);
    if(i1!=0) if(!toleq(vim, gel(gel(U, 3), i1), tol, prec)) i1=0;//Just because the angles are equal, the points don't have to be (though this occurence is expected to be extremely rare).
    if(i==1) ind2=lU-1;
    else ind2=i-1;//The two vertices of the side, this is the second one
    vim=mat_eval(gel(gel(U, 5), i), gel(gel(U, 3), ind2));//The second new vertex
    if(i1!=0){//If ind2 is paired, it MUST be paired to vertex i1+1
      i2=i1+1;
      if(i2==lU) i2=1;//Wrap around
      if(!toleq(vim, gel(gel(U, 3), i2), tol, prec)) i2=0;
      if(i2!=0){
        if(i<=i2) vectrunc_append(pair, mkvecsmall2(i, i2));//i<=i2 is put so that we are not pairing things twice.
      }
      else vectrunc_append(unpair, mkvecsmall2(i, ind2));//Second vertex not paired
    }
    else{//ind1 not paired
      vimang=shiftangle(garg(vim, prec), baseangle, tol, prec);//The second new angle
      i2=gen_search(vangles, vimang, 0, &toldata, &tolcmp_sort);
      if(i2!=0) if(!toleq(vim, gel(gel(U, 3), i2), tol, prec)) i2=0;//Just because the angles are equal, the points don't have to be (though this occurence is expected to be extremely rare).
      if(i2!=0) vectrunc_append(unpair, mkvecsmall2(i, ind1));//First vtx not paired
       else vectrunc_append(unpair, mkvecsmall3(i, ind1, ind2));//Neither vtx paired
    }
  }
  if(rboth) return gerepilecopy(top, mkvec2(pair, unpair));
  if(lg(unpair)==1){
    GEN pairvs=cgetg(lU, t_VECSMALL);
    for(long i=1;i<lg(pair);i++){
      pairvs[gel(pair, i)[1]]=gel(pair, i)[2];
      pairvs[gel(pair, i)[2]]=gel(pair, i)[1];
    }
    return gerepilecopy(top, pairvs);//No unpaired vertices
  }
  return gerepilecopy(top, unpair);//Unpaired vertices!
}

//Computes the isometric circle for g, returning [g, image in PSU(1, 1), circle]. Must pass in mats (psltopsu_transmats(p)), and a method that translates g to an element of PSL(2, R).
GEN isometriccircle_mats(GEN g, GEN mats, GEN *data, GEN (*gamtopsl)(GEN *, GEN, long), GEN tol, long prec){
  pari_sp top=avma;
  GEN ginpsl=gamtopsl(data, g, prec);
  GEN ret=cgetg(4, t_VEC);
  gel(ret, 1)=gcopy(g);
  gel(ret, 2)=psltopsu_mats(ginpsl, mats);
  gel(ret, 3)=isometriccircle_psu(gel(ret, 2), tol, prec);
  if(!gequal(gel(ret, 3), gen_m1)) return gerepileupto(top, ret);//Everything A-OK
  //If we reach here, we did not have enough precision
  long newprec=prec;
  pari_CATCH(CATCH_ALL){
    avma=top;
    pari_CATCH_reset();
    pari_err(e_MISC,"Could not increase precision enough. Please increase precision/memory");
    return gen_0;
  }
  pari_TRY{
    do{
      if(newprec-prec==5) pari_err(e_MISC,"Throw");
      avma=top;
      newprec++;//Increase precision
      tol=deftol(newprec);
      if(precision(gel(mats, 3))>0){//p is inexact
        GEN p=gtofp(gel(mats, 3), newprec);
        mats=psltopsu_transmats(p);//Updating mats to more precision
      }
      ginpsl=gamtopsl(data, g, newprec);
      ret=cgetg(4, t_VEC);
      gel(ret, 1)=gcopy(g);
      gel(ret, 2)=psltopsu_mats(ginpsl, mats);
      gel(ret, 3)=isometriccircle_psu(gel(ret, 2), tol, newprec);
    }
    while(gequal(gel(ret, 3), gen_m1));
    pari_CATCH_reset();
    return gerepileupto(top, ret);
  }
  pari_ENDCATCH
}

//Given an element g of PSU(1, 1), this returns the isometric circle associated to it.
GEN isometriccircle_psu(GEN g, GEN tol, long prec){
  pari_sp top=avma;
  if(toleq(gen_0, gcoeff(g, 2, 1), tol, prec)) return gen_0;//Isometric circle is everything, don't want to call it here.
  GEN geod=zerovec(8);
  gel(geod, 2)=gdivsg(1, gcoeff(g, 2, 1));//Need to take absolute value
  gel(geod, 1)=gneg(gmul(gcoeff(g, 2, 2), gel(geod, 2)));//-g[2,2]/g[2,1], the centre of the circle
  gel(geod, 2)=gabs(gel(geod, 2), prec);//We do things in this order to save a division.
  gel(geod, 7)=gen_1;//Always oriented counterclockwise
  pari_CATCH(CATCH_ALL){
    avma=top;
    pari_CATCH_reset();
    return gen_m1;//We increase precision in g and retry.
  }
  pari_TRY{
    GEN ipts=circle_int(geod, mkvec3s(0, 1, 0), tol, prec);//Intersect with x^2+y^2=1
    GEN ang=anglediff(garg(gsub(gel(ipts, 2), gel(geod, 1)), prec), garg(gsub(gel(ipts, 1), gel(geod, 1)), prec), tol, prec);
    if(gcmp(ang, mppi(prec))==1){//Properly orienting the start and endpoints
      gel(geod, 3)=gel(ipts, 2);
    gel(geod, 4)=gel(ipts, 1);
    }
    else{
      gel(geod, 3)=gel(ipts, 1);
      gel(geod, 4)=gel(ipts, 2);
    }
    gel(geod, 5)=radialangle(geod, gel(geod, 3), tol, prec);//Start angle
    gel(geod, 6)=shiftangle(radialangle(geod, gel(geod, 4), tol, prec), gel(geod, 5), tol, prec);//End angle
    pari_CATCH_reset();
    return gerepilecopy(top, geod);
  }
  pari_ENDCATCH
}

//Returns the normalized basis of G. Follows Algorithm 4.7 of Voight. Can pass in Ubase as a normalized boundary to append to, or Ubase=0 means we just start with G.
GEN normalizedbasis(GEN G, GEN Ubase, GEN mats, GEN gamid, GEN *data, GEN (*gamtopsl)(GEN *, GEN, long), GEN (*eltmul)(GEN *, GEN, GEN), GEN (*eltinv)(GEN *, GEN), int (*istriv)(GEN *, GEN), GEN tol, long prec){
  pari_sp top=avma, mid;
  long w=lg(G);
  GEN Gwithinv=cgetg(2*w-1, t_VEC);
  for(long i=1;i<lg(G);i++){
    gel(Gwithinv, i)=gel(G, i);
    gel(Gwithinv, w)=eltinv(data, gel(G, i));
    w++;
  }
  GEN U;
  if(gequal0(Ubase)) U=normalizedboundary(Gwithinv, mats, gamid, data, gamtopsl, tol, prec);//Step 2 when Ubase=0
  else U=normalizedboundary_givenU(Ubase, Gwithinv, mats, gamid, data, gamtopsl, tol, prec);//Step 2
  if(lg(gel(U, 1))==1) return gerepileupto(top, U);//No iso circles, go back now.
  GEN Gadd=Gwithinv, gbardat, gbar, unpair, v, g, scale=gdivgs(stoi(9), 10), Gaddnew, Uold=gel(U, 1);//Uold will track if we need to go back to step 3 (from step 4) or increase the scale (in step 5)
  long lunp, gind, vind;
  for(;;){
    if(gc_needed(top, 2)){//Garbage collection
      mid=avma;
      Uold=gcopy(Uold);
      U=gcopy(U);
      scale=gcopy(scale);
      Gadd=gcopy(Gadd);
      gerepileallsp(top, mid, 4, &Uold, &U, &scale, &Gadd);
    }
    Gaddnew=vectrunc_init(2*lg(Gadd));//The reductions to add to G
    for(long i=1;i<lg(Gadd);i++){//Doing step 3
      gbardat=reduceelt_givennormbound(U, gel(Gadd, i), gen_0, gamid, data, gamtopsl, eltmul, tol, prec);//Finding gbar=red_U(g).
      gbar=gel(gbardat, 1);
      if(!istriv(data, gbar)){
        vectrunc_append(Gaddnew, eltinv(data, gbar));//If not trivial, add the inverse.
      }
    }
    if(lg(Gaddnew)!=1){//We add Gaddnew to G, compute the normalized boundary, and go back to step 3 if U is changed
      Gadd=Gaddnew;
      U=normalizedboundary_givenU(U, Gadd, mats, gamid, data, gamtopsl, tol, prec);//Adding Gaddnew to U, recomputing normalized boundary
      if(lg(Uold)==lg(gel(U, 1))){
        if(!gequal(Uold, gel(U, 1))){Uold=gel(U, 1);continue;}//Go back to step 3.
      }
    }
    //Step 5. Find unpaired vertices, add them, go back!
    unpair=edgepairing(U, tol, 0, prec);//We only want the unpaired edges if they exist, and the paired only if there are no unpaired
    lunp=lg(unpair);
    if(typ(unpair)==t_VECSMALL){gel(U, 7)=unpair;break;}//Done!
    Gaddnew=vectrunc_init(2*lunp-1);//We add up to 2 new elements per unpaired vertex (there are lunp-1 such)
    for(long i=1;i<lunp;i++){
      gind=gel(unpair, i)[1];
      for(long k=2;k<lg(gel(unpair, i));k++){//The 1 or 2 vertices
        vind=gel(unpair, i)[k];
        v=gel(gel(U, 3), vind);
        if(toleq(gabs(v, prec), gen_1, tol, prec)){//Infinite vertex
          if(gind==vind) v=normalizedbasis_shiftpoint(gel(gel(U, 2), gind), scale, 1, prec);
          else v=normalizedbasis_shiftpoint(gel(gel(U, 2), gind), scale, 0, prec);
        }
        g=gel(gel(U, 1), gind);
        gbardat=reduceelt_givennormbound(U, g, v, gamid, data, gamtopsl, eltmul, tol, prec);//Reduce with respect to v (and not 0 like in step 3).
        vectrunc_append(Gaddnew, gel(gbardat, 1));
      } 
    }
    Gadd=Gaddnew;
    U=normalizedboundary_givenU(U, Gadd, mats, gamid, data, gamtopsl, tol, prec);//Adding Gadd to U, recomputing normalized boundary.
    if(lg(gel(U, 1))==lg(Uold)){//They might be equal
      if(gequal(gel(U, 1), Uold)){//They are equal. There must be an error where we have an oo vertex which gets shifted to being in the INTERIOR, i.e. we scaled down too much. To fix, we increase scale
        scale=gdivgs(gaddgs(scale, 9), 10);
      }
    }
    Uold=gel(U, 1);
  }
  return gerepilecopy(top, U);
}

//Given a circle arc c and either the initial(=1) or terminal(=0) point, this returns a nearby point on the arc (the new angle is r of the way between them; want r to be close to 1).
static GEN normalizedbasis_shiftpoint(GEN c, GEN r, int initial, long prec){
  pari_sp top=avma;
  GEN newang;
  if(initial==0) newang=gadd(gel(c, 5), gmul(gsub(gel(c, 6), gel(c, 5)), r));//a+(b-a)*r with a the initial and b the terminal angle
  else newang=gsub(gel(c, 6), gmul(gsub(gel(c, 6), gel(c, 5)), r));//b-(b-a)*r with a the initial and b the terminal angle
  return gerepileupto(top, gadd(gel(c, 1), gmul(expIr(gtofp(newang, prec)), gel(c, 2))));//c[1]+c[2]*e^(I*newang), the shifted point.
}

//Given a normalized boundary U, this appends the elements of G to it (G is the set of [elt, image in PSU(1, 1), isometric circle].
static GEN normalizedboundary_append(GEN Ubase, GEN G, GEN mats, GEN id, GEN tol, long prec){
  pari_sp top=avma, mid;
  long nGp1=lg(G), nG=nGp1-1;//nG=number of elements in G
  long nUp1=lg(gel(Ubase, 1)), nU=nUp1-1;//nUp1=number of elts in Ubase+1
  if(nU==0) return normalizedboundary_givencircles(G, mats, id, tol, prec);//Ubase is trivial.
  long maxsides=2*nG+nUp1+4;//maxsides is 5+the maximal number of sides we can have (every side exists and there are two sides added for each element of G (the side and an infinite side). The 5 is for added security (we loop around)
  GEN U=cgetg(maxsides, t_VECSMALL);//Stores the indices in U in order. The initial ordering may not have v_1 correct (may need to do a cyclic shift at the end).
  GEN vertices=cgetg(maxsides, t_VEC);//Stores the vertices in order; an entry is [vertex, radial angle]
  GEN pi=mppi(prec);//Pi
  GEN moo=mkmoo();
  GEN inter, ang, ang1, ang2, sidecirc, sidecirctermang, Ltermang, ten=stoi(10);
  GEN L=gel(gel(Ubase, 2), 1);//The current segment we are looking for intersections with.
  GEN baseang=garg(gel(L, 4), prec);//Angle to the terminal point of the first element of Ubase
  GEN Utermangles=cgetg(nUp1, t_VEC);//The angles to the terminal points in Ubase
  for(long i=1;i<nUp1;i++){
    L=gel(gel(Ubase, 2), i);
    if(gequal0(L)) gel(Utermangles, i)=gen_0;//infinite side, ignoring
    else gel(Utermangles, i)=shiftangle(garg(gel(L, 4), prec), baseang, tol, prec);//Angle to origin in [baseind, baseind+2*Pi)
  }
  GEN Gtermangles=cgetg(nGp1, t_VEC);
  for(long i=1;i<nGp1;i++){
    if(gequal0(gel(gel(G, i), 3))) gel(Gtermangles, i)=moo;//Does not give rise to a circle, want to ignore. -oo will ALWAYS be less than the start angle.
    else gel(Gtermangles, i)=shiftangle(garg(gel(gel(gel(G, i), 3), 4), prec), baseang, tol, prec);//Angle to origin in [baseind, baseind+2*Pi)
  }
  GEN Gord=indexsort(Gtermangles);//The order in which we want to look at the elements of G.
  long Gordind=0;
  for(long i=1;i<nGp1;i++){
    if(!gequal(gel(Gtermangles, Gord[i]), moo)){Gordind=i;break;}
  }//Moving past the -2's, i.e. elements of G giving no circle. These occur first as the other angles are >-2.
  if(Gordind==0) return gerepilecopy(top, Ubase);//No new circles.
  U[1]=1;
  gel(vertices, 1)=mkvec2(gel(gel(gel(Ubase, 2), 1), 4), baseang);//The first vertex is initially set to be the start angle of the first side.
  long ulen=1;//The current length of U and vertices. If we have to delete some vertices, this can decrease.
  long side, sidem1;//We basically re-insert the first side back at the end
  GEN gang=gel(Gtermangles, Gord[Gordind]);//The angle to the terminal side of the current element of G.
  int lastsidenew=0, newsidefromG, finalstretch=0, leftoverGs=0;
  long startpt=1;//Stores what will be the start point.
  for(long sid=2;sid<=2*nU;sid++){//We move through Ubase, stopping to insert the elements of G if need be.
    if(sid>nU){
      if(!gequal(gang, ten)) leftoverGs=1;//There are still G's to do, and we have looped around
      else leftoverGs=0;
      side=sid-nU;
      if(side==1) sidem1=nU;
      else sidem1=side-1;
      finalstretch=1;//We stop once the (looped around) side gets re-inserted in.
    }
    else{
      side=sid;
      sidem1=sid-1;
    }
    mid=avma;
	if(finalstretch && !leftoverGs){//We go back and hit the start of U
	  if(U[side]==0){startpt++;continue;}
	  else if(U[side]>0){
		sidecirc=gel(gel(Ubase, 2), U[side]);
		sidecirctermang=gel(Utermangles, U[side]);
		newsidefromG=0;
	  }
	  else{
		sidecirc=gel(gel(G, -U[side]), 3);
		sidecirctermang=gel(Gtermangles, -U[side]);
		newsidefromG=1;
	  }
	  if(U[ulen]>0){
		L=gel(gel(Ubase, 2), U[ulen]);
		Ltermang=gel(Utermangles, U[ulen]);
	  }
	  else{
		L=gel(gel(G, -U[ulen]), 3);
		Ltermang=gel(Gtermangles, -U[ulen]);
	  }
	}
	else if(gequal0(gel(gel(Ubase, 2), side))) continue;//We skip past infinite sides of Ubase.
    else if(lastsidenew==0){//Working on a non-infinite side of Ubase, and the last side was also a side of Ubase
      if(tolcmp(gang, gel(Utermangles, side), tol, prec)==1 && !leftoverGs){//Consecutive old sides
	    if(gequal0(gel(gel(Ubase, 2), sidem1))){//The last side was infinite and skipped over; must re-insert it.
          ulen++;
          U[ulen]=0;
          gel(vertices, ulen)=mkvec2(gel(gel(Ubase, 3), sidem1-1), gel(gel(Ubase, 4), sidem1-1));//The infinite side
		  ulen++;
          U[ulen]=side;
          gel(vertices, ulen)=mkvec2(gel(gel(Ubase, 3), sidem1), gel(gel(Ubase, 4), sidem1));
          continue;//Go on
        }
		//Now we check if the new vertex is beyond the old one
	    ang=gel(gel(Ubase, 4), sidem1);//Angle to the intersection point
		ang1=anglediff(ang, gel(gel(vertices, ulen), 2), tol, prec);//ang-angle to the previous vertex.
        if(gequal0(ang1) || gcmp(ang1, pi)>=0){//Delete last side and go backwards. The previous side MUST be a new side.
		  avma=mid;
		  ulen--;
		  lastsidenew=1;
		  sid--;
		  continue;
		}
        ulen++;
        U[ulen]=side;
        gel(vertices, ulen)=mkvec2(gel(gel(Ubase, 3), sidem1), gel(gel(Ubase, 4), sidem1));
        continue;//Go on
      }
      //Here, this means that we need to try to insert G[Gord[Gordind]] BEFORE the current side.
      sidecirc=gel(gel(G, Gord[Gordind]), 3);//The new side we are trying to insert
      sidecirctermang=gang;//Terminal angle
      L=gel(gel(Ubase, 2), sidem1);//The previous side
      if(gequal0(L)){L=gel(gel(Ubase, 2), sidem1-1);Ltermang=gel(Utermangles, sidem1-1);}//If the previous side was oo, we go back one further, as this side is not oo.
      else Ltermang=gel(Utermangles, sidem1);
      newsidefromG=1;
    }
    else{//Now we check if we are inserting a new side or an old side
      if(tolcmp(gang, gel(Utermangles, side), tol, prec)==1 && !leftoverGs){//Inserting old side
        sidecirc=gel(gel(Ubase, 2), side);//The old side we are trying to insert
        sidecirctermang=gel(Utermangles, side);//Terminal angle
        L=gel(gel(G, -U[ulen]), 3);//The newly inserted side
        Ltermang=gel(Gtermangles, -U[ulen]);
        newsidefromG=0;
      }
      else{
        sidecirc=gel(gel(G, Gord[Gordind]), 3);//The new side we are trying to insert
        sidecirctermang=gang;//Terminal angle
        L=gel(gel(G, -U[ulen]), 3);//The newly inserted side
        Ltermang=gel(Gtermangles, -U[ulen]);
        newsidefromG=1;
      }
    }
    //At this point, we are either trying to insert a new on an old, or an old on a new side, or a new on a new
    inter=arc_int(L, sidecirc, tol, prec);//Intersection of L and the next side.
    if(lg(inter)==1){//It did NOT intersect.
      ang2=garg(gel(L, 3), prec);
      ang1=anglediff(ang2, Ltermang, tol, prec);//The angle between the terminal and initial points of L
      ang=anglediff(garg(gel(sidecirc, 3), prec), Ltermang, tol, prec);//Angle to the initial point of sidecirc from the terminal angle of L.
      if(tolcmp(ang, ang1, tol, prec)<=0){//sidecirc is contained within L
		if(finalstretch && !leftoverGs){
		  if(!gequal(L, sidecirc)){
		    avma=mid;
			startpt++;//We add 1 to startpt to signal that we must start at a different point.
			continue;
		  }
		  //L=sidecirc, and so we are actually done. This should only happen when Ubase has 1 non-trivial side, and we don't end up adding in any sides that can do better (i.e. the output is the same as the input). So we just continue on and let the rest do its thing.
		}
        else if(newsidefromG){//Failed to insert since it did not help.
		  avma=mid;
          Gordind++;
          if(Gordind==nGp1){gang=ten;}//We are done with G, so we just want to start appending old sides.
          else gang=gel(Gtermangles, Gord[Gordind]);
          sid--;//We need to try again with the current side since we "jumped the line" with the element of G.
		  continue;
        }//We also want to leave lastsidenew unchanged, as we did not insert
        else{avma=mid;continue;}
      }
      //We have two new sides: a side at infinity, and this side.
      ulen++;
      U[ulen]=0;
      gel(vertices, ulen)=mkvec2(gel(L, 3), ang2);
      ulen++;
	  gel(vertices, ulen)=mkvec2(gel(sidecirc, 4), sidecirctermang);
	  if(finalstretch && !leftoverGs){
	    U[ulen]=U[side];
		break;//Done!
	  }
      else if(newsidefromG){
        U[ulen]=-Gord[Gordind];
        Gordind++;
        if(Gordind==nGp1){gang=ten;}//We are done with G, so we just want to start appending old sides.
        else gang=gel(Gtermangles, Gord[Gordind]);
        sid--;//Redo this side
        lastsidenew=1;
      }
      else{
        U[ulen]=side;
        lastsidenew=0;
      }
    }
    else{//It DID intersect
      //It may be that this new side actually comes into U[1] from the bottom. Then we need a side at oo
      ang=anglediff(sidecirctermang, Ltermang, tol, prec);//Angle from the terminal angle of the last side to the terminal angle of the new side.
      if(gequal0(ang)){
        if(tolcmp(gel(sidecirc, 2), gel(L, 2), tol, prec)<=0){//This side lies inside the previous one, continue on (compared radii).
          avma=mid;
		  if(finalstretch && !leftoverGs) startpt++;
          else if(newsidefromG){//Failed to insert since it did not help.
            Gordind++;
            if(Gordind==nGp1){gang=ten;}//We are done with G, so we just want to start appending old sides.
            else gang=gel(Gtermangles, Gord[Gordind]);
            sid--;//We need to try again with the current side since we "jumped the line" with the element of G.
          }//We also want to leave lastsidenew unchanged, as we did not insert
          continue;
        }
		//Only need to update U[ulen]; the intersection point is the same.
		if(newsidefromG){U[ulen]=-Gord[Gordind];lastsidenew=1;}
		else{U[ulen]=side;lastsidenew=0;}
		continue;
      }
      else if(gcmp(ang, pi)==1){//We DID come in from below
        ulen++;
        U[ulen]=0;//Side at oo
        gel(vertices, ulen)=mkvec2(gel(L, 3), garg(gel(L, 3), prec));//Side at oo
        ulen++;
		gel(vertices, ulen)=mkvec2(gel(sidecirc, 4), sidecirctermang);
		if(finalstretch && !leftoverGs){//Done!
		  U[ulen]=U[side];
		  break;
		}
        if(newsidefromG){
          U[ulen]=-Gord[Gordind];
          Gordind++;
          if(Gordind==nGp1){gang=ten;}//We are done with G, so we just want to start appending old sides.
          else gang=gel(Gtermangles, Gord[Gordind]);
          sid--;//Redo this side
          lastsidenew=1;
        }
        else{
          U[ulen]=side;
          lastsidenew=0;
        }
        continue;
      }
      //Now we are sure did not come in from below.
      inter=gel(inter, 1);//The point
      if(toleq(inter, gel(L, 3), tol, prec)){//The side lies entirely in the previous side OR touches it at the end
	    if(toleq(inter, gel(sidecirc, 3), tol, prec)){//Lies inside
          avma=mid;
		  if(finalstretch && !leftoverGs) startpt++;
          else if(newsidefromG){//Failed to insert since it did not help.
            Gordind++;
            if(Gordind==nGp1){gang=ten;}//We are done with G, so we just want to start appending old sides.
            else gang=gel(Gtermangles, Gord[Gordind]);
            sid--;//We need to try again with the current side since we "jumped the line" with the element of G.
          }//We also want to leave lastsidenew unchanged, as we did not insert
          continue;
		}
      }
      //Now we have a proper "normal" intersection
      ang1=garg(inter, prec);
      ang=anglediff(ang1, gel(gel(vertices, ulen), 2), tol, prec);//Angle to the new vtx from the previous as a bases
      if(gcmp(ang, pi)!=-1 || toleq(ang, gen_0, tol, prec)){//We must go backwards!
        while(ulen>1){
          ulen--;
          if(U[ulen]>0) L=gel(gel(Ubase, 2), U[ulen]);
          else L=gel(gel(G, -U[ulen]), 3);//U[ulen]=0 is impossible; it is guarenteed that we do not backtrack past an infinite side in this algorithm.
          inter=gel(arc_int(L, sidecirc, tol, prec), 1);//They MUST intersect
          ang1=garg(inter, prec);
          ang=anglediff(ang1, gel(gel(vertices, ulen), 2), tol, prec);
          if(gcmp(ang, pi)!=-1 || toleq(ang, gen_0, tol, prec)) continue;//Keep going back
          break;//At this point we have reached where we need to insert the new side.
        }
      }
      //Now we are ready to insert it.
	  if(finalstretch && !leftoverGs){//Done if we intersected before the original intersection, and must continue on if not.
		if(gcmp(pi, anglediff(ang1, gel(gel(vertices, side+1), 2), tol, prec))==1){startpt++;continue;}//We have superseeded the previous side, increase start point.
		ulen++;
		U[ulen]=U[side];
	    gel(vertices, ulen)=mkvec2(inter, ang1);
		break;
	  }
	  ulen++;
	  gel(vertices, ulen)=mkvec2(inter, ang1);
      if(newsidefromG){
        U[ulen]=-Gord[Gordind];
        Gordind++;
        if(Gordind==nGp1){gang=ten;}//We are done with G, so we just want to start appending old sides.
        else gang=gel(Gtermangles, Gord[Gordind]);
        sid--;//Redo this side
        lastsidenew=1;
      }
      else{
        U[ulen]=side;
        lastsidenew=0;
      }
    }
  }
  //We potentially adjust the start position.
  long best=ulen;//How many sides we have (plus 1).
  mid=avma;
  long k=ulen-1;
  GEN L0=zerovec(8);gel(L0, 4)=gen_1;gel(L0, 7)=gen_1;gel(L0, 8)=gen_1;//L0=line segment [0, 1].
  GEN mininter=gen_2;
  ang1=anglediff(gel(Utermangles, 1), gen_0, tol, prec);
  while(U[k]<=0){//U[1] in the start position can only be superseeded by a new edge.
    if(U[k]==0){k--;continue;}//We may have added in a 0 edge
    L=gel(gel(G, -U[k]), 3);
    inter=arcseg_int(L, L0, tol, prec);
    if(lg(inter)!=1){//Intersection
      inter=real_i(gel(inter, 1));//Intersected [0, 1]
      if(tolcmp(inter, mininter, tol, prec)==-1){best=k;mininter=inter;}
	  else break;//We cannot do better than this.
    }
    else{//No intersection
      if(gequal(mininter, gen_2)) best=k;//Might have to take the side with smallest angle
      else break;//Done, as we had intersections, and now there are none.
      ang=anglediff(gel(Gtermangles, -U[k]), gen_0, tol, prec);
      if(tolcmp(ang, ang1, tol, prec)!=-1){best++;break;}//Done, and the best was actually one ago.
    }
    k--;
  }
  if(U[startpt]>0) L=gel(gel(Ubase, 2), U[startpt]);
  else L=gel(gel(G, -U[startpt]), 3);
  inter=arcseg_int(L, L0, tol, prec);
  if(lg(inter)>1){//Intersect. Now we need to go forward, as it is possible that the best intersection with [0, 1] was added right after the first side.
    inter=real_i(gel(inter, 1));
    if(tolcmp(mininter, inter, tol, prec)!=-1){best=ulen;mininter=inter;}//The first arc was better
	k=2;
	while(U[k]<0){
	  L=gel(gel(G, -U[k]), 3);
      inter=arcseg_int(L, L0, tol, prec);
	  if(lg(inter)==1) break;//No, intersection, done
      inter=real_i(gel(inter, 1));//Intersected [0, 1]
      if(tolcmp(inter, mininter, tol, prec)!=1){best=k;mininter=inter;}
	  else break;//We cannot do better than this.
	  k++;
	}
  }
  else{
	if(gequal(mininter, gen_2) && U[k]<0) best++;//There was NO intersection with [0, 1], so we are ending up on an infinite side! This is not right, so we must increment it by one. If best=ulen then we did not boop it backward, so don't need to increment if forward
  }
  
  avma=mid;
  long np1=ulen-startpt+1;
  GEN firstang=gel(gel(vertices, startpt+1), 2);//The angle to the first vertex
  //By wrapping back around, we have ulen-startpt sides: the last side is the same as the first.
  GEN ret=cgetg(NORMBOUND, t_VEC);
  for(long i=1;i<=5;i++) gel(ret, i)=cgetg(np1, t_VEC);//elements, icircs, vertices, matrices, term angles. Places 6, 7 will store the area and side pairing (0 for now)
  long j=1, h;
  for(long i=best;i<ulen;i++){
    if(U[i]==0){//Side at oo
      gel(gel(ret, 1), j)=gcopy(id);//Element
      gel(gel(ret, 2), j)=gen_0;//No circle
      gel(gel(ret, 3), j)=gcopy(gel(gel(vertices, i+1), 1));//Vertex
      gel(gel(ret, 4), j)=shiftangle(gel(gel(vertices, i+1), 2), firstang, tol, prec);//Vertex angle
      gel(gel(ret, 5), j)=gen_0;//No matrix
    }//Now we have a real side
    else if(U[i]<0){
      h=-U[i];
      gel(gel(ret, 1), j)=gcopy(gel(gel(G, h), 1));//Element
      gel(gel(ret, 2), j)=gcopy(gel(gel(G, h), 3));//Circle
      gel(gel(ret, 3), j)=gcopy(gel(gel(vertices, i+1), 1));//Vertex
      gel(gel(ret, 4), j)=shiftangle(gel(gel(vertices, i+1), 2), firstang, tol, prec);//Vertex angle
      gel(gel(ret, 5), j)=gcopy(gel(gel(G, h), 2));//Matrix
    }
    else{
      gel(gel(ret, 1), j)=gcopy(gel(gel(Ubase, 1), U[i]));//Element
      gel(gel(ret, 2), j)=gcopy(gel(gel(Ubase, 2), U[i]));//Circle
      gel(gel(ret, 3), j)=gcopy(gel(gel(vertices, i+1), 1));//Vertex
      gel(gel(ret, 4), j)=shiftangle(gel(gel(vertices, i+1), 2), firstang, tol, prec);//Vertex angle
      gel(gel(ret, 5), j)=gcopy(gel(gel(Ubase, 5), U[i]));//Matrix
    }
    j++;
  }
  for(long i=startpt;i<best;i++){
    if(U[i]==0){//Side at oo
      gel(gel(ret, 1), j)=gcopy(id);//Element
      gel(gel(ret, 2), j)=gen_0;//No circle
      gel(gel(ret, 3), j)=gcopy(gel(gel(vertices, i+1), 1));//Vertex
      gel(gel(ret, 4), j)=shiftangle(gel(gel(vertices, i+1), 2), firstang, tol, prec);//Vertex angle
      gel(gel(ret, 5), j)=gen_0;//No matrix
    }//Now we have a real side
    else if(U[i]<0){
      h=-U[i];
      gel(gel(ret, 1), j)=gcopy(gel(gel(G, h), 1));//Element
      gel(gel(ret, 2), j)=gcopy(gel(gel(G, h), 3));//Circle
      gel(gel(ret, 3), j)=gcopy(gel(gel(vertices, i+1), 1));//Vertex
      gel(gel(ret, 4), j)=shiftangle(gel(gel(vertices, i+1), 2), firstang, tol, prec);//Vertex angle
      gel(gel(ret, 5), j)=gcopy(gel(gel(G, h), 2));//Matrix
    }
    else{
      gel(gel(ret, 1), j)=gcopy(gel(gel(Ubase, 1), U[i]));//Element
      gel(gel(ret, 2), j)=gcopy(gel(gel(Ubase, 2), U[i]));//Circle
      gel(gel(ret, 3), j)=gcopy(gel(gel(vertices, i+1), 1));//Vertex
      gel(gel(ret, 4), j)=shiftangle(gel(gel(vertices, i+1), 2), firstang, tol, prec);//Vertex angle
      gel(gel(ret, 5), j)=gcopy(gel(gel(Ubase, 5), U[i]));//Matrix
    }
    j++;
  }
  gel(ret, 6)=hpolygon_area(gel(ret, 2), gel(ret, 3), tol, prec);//Area
  gel(ret, 7)=gen_0;
  gel(ret, 8)=gcopy(mats);
  return gerepileupto(top, ret);
}

//Initializes the inputs for normalizedboundary_append. Works BEST if p is given as an EXACT number.
static GEN normalizedboundary_givenU(GEN Ubase, GEN G, GEN mats, GEN id, GEN *data, GEN (*gamtopsl)(GEN *, GEN, long), GEN tol, long prec){
  pari_sp top=avma;
  long lx;
  GEN Gnew=cgetg_copy(G, &lx);
  int skipped=0;
  for(long i=1;i<lx;i++){
    pari_CATCH(CATCH_ALL){
      gel(Gnew, i)=mkvec3(gen_0, gen_0, gen_0);//The third entry being 0 means this gets ignored
      skipped++;
    }
    pari_TRY{
      gel(Gnew, i)=isometriccircle_mats(gel(G, i), mats, data, gamtopsl, tol, prec);
    }
    pari_ENDCATCH
  }
  if(skipped>0){
    char *warningmsg=pari_sprintf("%d isometric circles were skipped due to insufficient precision", skipped);
    pari_warn(warner, warningmsg);
    pari_free(warningmsg);
  }
  return gerepileupto(top, normalizedboundary_append(Ubase, Gnew, mats, id, tol, prec));
}

//G is the set of [elt, image in PSU(1, 1), isometric circle]. This returns the normalized boundary of the exterior domain. The output is [elements, icircs, vertices, vertex angles, matrices, area, 0, mats]. The circle corresponding to elements[i] is icircs[i], and the vertices are vertices[i-1] and vertices[i]. matrices[i] is the image in PSU(1,1) of elements[i]. The element 1 corresponds to a section on the unit circle, which also corresponds to a circle of 0. Vertex angles stores the radial angle to the ith vertex (with base angle being the first one). The area is the area, and the 0 stores the side pairing when we have a fundamental domain (so a priori stores nothing).
static GEN normalizedboundary_givencircles(GEN G, GEN mats, GEN id, GEN tol, long prec){
  pari_sp top=avma, mid;
  long np1=lg(G), n=np1-1, twonp3=2*n+3;
  GEN U=cgetg(twonp3, t_VECSMALL);//Stores the indices in U in order. The initial ordering may not have v_1 correct (may need to do a cyclic shift at the end). We double since there are at most n sides coming from G and n sides coming from the unit circle.
  GEN vertices=cgetg(twonp3, t_VEC);//Stores the vertices in order; an entry is [vertex, radial angle]
  GEN pi=mppi(prec), moo=mkmoo();//Pi, -oo
  GEN inter, ang, sidecirc, ang1;
  //Finding the first element.
  mid=avma;
  GEN L0=zerovec(8);gel(L0, 4)=gen_1;gel(L0, 7)=gen_1;gel(L0, 8)=gen_1;//L0=line segment [0, 1].
  long hminind=0;
  int isnew;
  GEN Hmin=gen_1;
  for(long i=1;i<np1;i++){//We start by finding intersections with L.
    sidecirc=gel(gel(G, i), 3);
    if(gequal0(sidecirc)) continue;//Ignore elts of G giving no circle.
    inter=arcseg_int(sidecirc, L0, tol, prec);
    if(lg(inter)==1) continue;//No intersection
    inter=real_i(gel(inter, 1));//Take the real part since it is actually real, and need to get the type correct.
    isnew=tolcmp(inter, Hmin, tol, prec);
    if(isnew==-1){hminind=i;Hmin=inter;continue;}//New min
    else if(isnew==1) continue;//Not new
    //Now we = the min, so need to see which initial point is longer
    ang=anglediff(garg(gel(gel(gel(G, i), 3), 3), prec), garg(gel(gel(gel(G, hminind), 3), 3), prec), tol, prec);
    if(gcmp(ang, pi)==-1) hminind=i;//Better min.
  }
  avma=mid;//Don't need these calcs other than hminind
  GEN baseang;
  if(hminind>0) baseang=garg(gel(gel(gel(G, hminind), 3), 4), prec);//The base angle to the terminal point.
  else baseang=gen_0;//We start from the 0 angle instead.
  GEN termangles=cgetg(np1, t_VEC);
  for(long i=1;i<np1;i++){
    if(gequal0(gel(gel(G, i), 3))) gel(termangles, i)=moo;//Does not give rise to a circle, want to ignore. -oo will ALWAYS be less than the start angle.
    else gel(termangles, i)=shiftangle(garg(gel(gel(gel(G, i), 3), 4), prec), baseang, tol, prec);//Angle to origin in [baseind, baseind+2*Pi)
  }
  //We now order G by the angle to the terminal points.
  GEN ordering=indexsort(termangles);//The order in which we want to look at the elements of G.
  long startind=0;
  for(long i=1;i<np1;i++){
    if(!gequal(gel(termangles, ordering[i]), moo)){startind=i;break;}
  }//Moving past the -2's, i.e. elements of G giving no circle. These occur first as the other angles are >=0. If hind!=0, then ordering[startind]=hind.
  //Now we start at G[ordering[ind]]. For the first element, we ONLY need to store the index.
  if(startind==0){//NO valid isometric circles inputted.
    avma=top;
    GEN retempty=cgetg(NORMBOUND, t_VEC);
    for(long i=1;i<=5;i++) gel(retempty, i)=cgetg(1, t_VEC);//elements, icircs, vertices, matrices, area, sidepairing
    gel(retempty, 6)=mkoo();
    gel(retempty, 7)=gen_0;
	gel(retempty, 8)=cgetg(1, t_VEC);
    return retempty;
  }
  if(hminind>0 && ordering[startind]!=hminind){//The first side has multiple circles coming out of its terminal point. I'm pretty sure this happens if and only if (well, in the Shimura curve case) Q is unramified everywhere.
    for(long i=startind+1;i<np1;i++){
      if(ordering[i]==hminind){ordering[i]=ordering[startind];ordering[startind]=hminind;break;}//Fix it
    }
  }
  U[1]=ordering[startind];
  gel(vertices, 1)=mkvec2(gel(gel(gel(G, ordering[startind]), 3), 4), gel(termangles, ordering[startind]));//The first vertex is initially set to be the start angle of the first side.
  long ulen=1;//The current length of U and vertices. If we have to delete some vertices, this can decrease.
  GEN L=gel(gel(G, ordering[startind]), 3), ang2;//The current segment we are looking for intersections with.
  long side;//We basically re-insert the first side back at the end
  for(long sid=startind+1;sid<=np1;sid++){//Doing things for the current side.
    if(sid==np1) side=startind;
    else side=sid;
    mid=avma;
    sidecirc=gel(gel(G, ordering[side]), 3);//The circle arc of the next side
    inter=arc_int(L, sidecirc, tol, prec);//Intersection of L and the next side.
    if(lg(inter)==1){//It did NOT intersect.
      ang1=garg(gel(L, 3), prec);
      ang=anglediff(ang1, gel(termangles, ordering[side]), tol, prec);//Angle to the initial point of L from the terminal point of the new side.
      if(gcmp(ang, pi)==-1 && ordering[side]!=U[ulen]){//the new side is contained entirely in L, OR comes in from below (the last check guarentees that it isn't ths same side; applicable when the final normalized boundary has 1 iso circle.
        ang2=garg(gel(L, 4), prec);
        ang=anglediff(gel(termangles, ordering[side]), ang2, tol, prec);//Angle to the terminal point of the new side from the terminal point of L
        if(gcmp(ang, pi)==-1){
          avma=mid;//the new side is contained entirely in L, discard and continue on (may as well reset avma).
          continue;
        }
        //Now we are actually coming in from below, so can continue on as normal.
      }
      //We have two new sides: a side at infinity, and this side.
      ulen++;
      U[ulen]=-1;
      gel(vertices, ulen)=mkvec2(gel(L, 3), ang1);
      ulen++;
      U[ulen]=ordering[side];
      gel(vertices, ulen)=mkvec2(gel(gel(gel(G, ordering[side]), 3), 4), gel(termangles, ordering[side]));//The terminal point of sidecirc is a vertex.
      L=gel(gel(G, ordering[side]), 3);//Setting L
    }
    else{//It DID intersect
      //It may be that this new side actually comes into U[1] from the bottom. Then we need a side at oo
      ang=anglediff(gel(termangles, U[ulen]), gel(termangles, ordering[side]), tol, prec);//Angle from the terminal angle of the last side to the terminal angle of the new side.
      if(gequal0(ang)){
        if(tolcmp(gel(gel(gel(G, ordering[side]), 3), 2), gel(gel(gel(G, U[ulen]), 3), 2), tol, prec)<=0){avma=mid;continue;}//This side lies inside the previous one, continue on (compared radii).
      }
      else if(gcmp(ang, pi)==-1){//We DID come in from below
        ulen++;
        U[ulen]=-1;//Side at oo
        gel(vertices, ulen)=mkvec2(gel(L, 3), garg(gel(L, 3), prec));//Side at oo
        ulen++;
        U[ulen]=ordering[side];
        gel(vertices, ulen)=mkvec2(gel(gel(gel(G, ordering[side]), 3), 4), gel(termangles, ordering[side]));
        L=gel(gel(G, ordering[side]), 3);
        continue;
      }
      //Now we are sure did not come in from below.
      inter=gel(inter, 1);//The point
      if(toleq(inter, gel(gel(gel(G, ordering[side]), 3), 3), tol, prec)){avma=mid;continue;}//The side lies entirely in the previous side
      ang1=garg(inter, prec);
      ang=anglediff(ang1, gel(gel(vertices, ulen), 2), tol, prec);//Angle to the new vtx from the previous as a bases
      if(gcmp(ang, pi)!=-1 || toleq(ang, gen_0, tol, prec)){//We must go backwards!
        while(ulen>1){
          ulen--;
          L=gel(gel(G, U[ulen]), 3);//The old side we look at
          inter=gel(arc_int(L, sidecirc, tol, prec), 1);//They MUST intersect
          ang1=garg(inter, prec);
          ang=anglediff(ang1, gel(gel(vertices, ulen), 2), tol, prec);
          if(gcmp(ang, pi)!=-1 || toleq(ang, gen_0, tol, prec)) continue;//Keep going back
          break;//At this point we have reached where we need to insert the new side.
        }
      }
      //Now we are ready to insert it.
      ulen++;
      L=sidecirc;
      U[ulen]=ordering[side];
      gel(vertices, ulen)=mkvec2(inter, ang1);
    }
  }
  GEN firstang=gel(gel(vertices, 2), 2);//The angle to the first vertex
  //By wrapping back around, we have ulen-1 sides: the last side is the same as the first.
  GEN ret=cgetg(NORMBOUND, t_VEC);
  for(long i=1;i<=5;i++) gel(ret, i)=cgetg(ulen, t_VEC);//elements, icircs, vertices, matrices, term angles. Places 6, 7 will store the area and side pairing (0 for now)
  for(long i=1;i<ulen;i++){
    if(U[i]==-1){//Side at oo
      gel(gel(ret, 1), i)=gcopy(id);//Element
      gel(gel(ret, 2), i)=gen_0;//No circle
      gel(gel(ret, 3), i)=gcopy(gel(gel(vertices, i+1), 1));//Vertex
      gel(gel(ret, 4), i)=shiftangle(gel(gel(vertices, i+1), 2), firstang, tol, prec);//Vertex angle
      gel(gel(ret, 5), i)=gen_0;//No matrix
      continue;
    }//Now we have a real side
    gel(gel(ret, 1), i)=gcopy(gel(gel(G, U[i]), 1));//Element
    gel(gel(ret, 2), i)=gcopy(gel(gel(G, U[i]), 3));//Circle
    gel(gel(ret, 3), i)=gcopy(gel(gel(vertices, i+1), 1));//Vertex
    gel(gel(ret, 4), i)=shiftangle(gel(gel(vertices, i+1), 2), firstang, tol, prec);//Vertex angle
    gel(gel(ret, 5), i)=gcopy(gel(gel(G, U[i]), 2));//Matrix
  }
  gel(ret, 6)=hpolygon_area(gel(ret, 2), gel(ret, 3), tol, prec);//Area
  gel(ret, 7)=gen_0;
  gel(ret, 8)=gcopy(mats);
  return gerepileupto(top, ret);
}

//Initializes the inputs for normalizedboundary_givencircles. Works BEST if p is given as an EXACT number.
GEN normalizedboundary(GEN G, GEN mats, GEN id, GEN *data, GEN (*gamtopsl)(GEN *, GEN, long), GEN tol, long prec){
  pari_sp top=avma;
  long lx;
  GEN Gnew=cgetg_copy(G, &lx);
  int skipped=0;
  for(long i=1;i<lx;i++){
    pari_CATCH(CATCH_ALL){
      gel(Gnew, i)=mkvec3(gen_0, gen_0, gen_0);//The third entry being 0 means this gets ignored
      skipped++;
    }
    pari_TRY{
      gel(Gnew, i)=isometriccircle_mats(gel(G, i), mats, data, gamtopsl, tol, prec);
    }
    pari_ENDCATCH
  }
  if(skipped>0){
    char *warningmsg=pari_sprintf("%d isometric circles were skipped due to insufficient precision", skipped);
    pari_warn(warner, warningmsg);
    pari_free(warningmsg);
  }
  return gerepileupto(top, normalizedboundary_givencircles(Gnew, mats, id, tol, prec));
}

//This returns [v, vecsmall(ind)], where the first side of U that c intersects has index ind, and the point of intersection is v. If start=1, we are searching near the start point of c, else we seach near the end point of c (c is either an arc or a segment with start/end points on the unit circle. If we intersect an oo side, returns [v, vecsmall(-1)] instead.
static GEN normalizedboundary_sideint(GEN U, GEN c, int start, GEN tol, long prec){
  pari_sp top=avma;
  GEN v, ret, inter, d1, d2;
  long ind;
  if(gequal(gel(c, 8), gen_1)){//Line segment; the index found from normalizedboundary_outside is correct guarenteed.
    if(start==1) v=gel(c, 3);
    else v=gel(c, 4);
    ind=normalizedboundary_outside(U, v, tol, prec);//This is the index
    if(ind==-1){//Infinite side
      ret=cgetg(3, t_VEC);
      gel(ret, 1)=gcopy(v);
      gel(ret, 2)=mkvecsmall(-1);
      return gerepileupto(top, ret); 
    }
    inter=arcseg_int(gel(gel(U, 2), ind), c, tol, prec);
    if(lg(inter)==2) v=gel(inter, 1);//1 intersection point, the most common
    else if(lg(inter)==3){
      d1=gabs(gsub(gel(inter, 1), v), prec);
      d2=gabs(gsub(gel(inter, 2), v), prec);
      if(gcmp(d1, d2)==-1) v=gel(inter, 1);//Finding which intersection is closer to v. Tolerance not required.
      else v=gel(inter, 2);
    }
    else pari_err_TYPE("Should have been intersections, but there were none", inter);//This should never happen, but here in case.
    ret=cgetg(3, t_VEC);
    gel(ret, 1)=gcopy(v);
    gel(ret, 2)=mkvecsmall(ind);
    return gerepileupto(top, ret);
  }
  GEN v2;//The other point
  if(start==1){
    if(gequal(gel(c, 7), gen_m1)){v=gel(c, 4);v2=gel(c, 3);}//If directed backwards, need the terminal point as the start
    else {v=gel(c, 3);v2=gel(c, 4);}//If undirected or directed forwards, the initial point is the start
  }
  else{
    if(gequal(gel(c, 7), gen_m1)){v=gel(c, 3);v2=gel(c, 4);}//If directed backwards, need the initial point.
    else{v=gel(c, 4);v2=gel(c, 3);}//Undirected/forwards, the terminal point is the start.
  }
  ind=normalizedboundary_outside(U, v, tol, prec);
  long lU=lg(gel(U, 1));//This is the index
  if(ind==-1){//Infinite side
    ret=cgetg(3, t_VEC);
    gel(ret, 1)=gcopy(v);
    gel(ret, 2)=mkvecsmall(-1);
    return gerepileupto(top, ret); 
  }
  //Circle arc (which should be the much more common case).
  GEN vang, ang, pi=mppi(prec);//Base angle+Pi
  int where, ind2;
  for(;;){//We intersect, see if it is on the sement or left or right of it. If left go to next index, if right go to previous
    inter=arc_int(gel(gel(U, 2), ind), c, tol, prec);
    if(lg(inter)==2) v=gel(inter, 1);//Only 1 intersection point, the most common case
    else if(lg(inter)==3){//Two intersection points
      d1=gabs(gsub(gel(inter, 1), v), prec);
      d2=gabs(gsub(gel(inter, 1), v2), prec);
      if(gcmp(d1, d2)==-1) v=gel(inter, 1);//Intersection[1] is closer to v than v2, so this is right.
      else v=gel(inter, 2);
    }
    else pari_err_TYPE("Should have been intersections, but there were none", inter);//This should never happen, but here in case.
    vang=garg(v, prec);//Angle to v
    ang=anglediff(gel(gel(U, 4), ind), vang, tol, prec);//Angle to the vertex from v.
    where=tolcmp(ang, pi, tol, prec);
    if(where==1){//Left of the vertex, continue
      ind++;
      if(ind==lU) ind=1;
      continue;
    }
    if(ind==1) ind2=lU-1;
    else ind2=ind-1;
    ang=anglediff(vang, gel(gel(U, 4), ind2), tol, prec);//Angle to v with reference to vertex 2
    where=tolcmp(pi, ang, tol, prec);//No need for cases this time
    if(where>=0) break;//Correct side!
    ind=ind2;//We must go backwards to the side ind2
  }
  ret=cgetg(3, t_VEC);
  gel(ret, 1)=gcopy(v);
  gel(ret, 2)=mkvecsmall(ind);
  return gerepileupto(top, ret);
}

//Returns -1 if z is in the interior of the normalized boundary or on the edge, and ind if in the exterior (i.e. between the boundary and the unit circle), where ind is the index of the side it is on (projecting from the origin to z). Does not behave well if z is a vertex of U that is on the unit disc.
static long normalizedboundary_outside(GEN U, GEN z, GEN tol, long prec){
  pari_sp top=avma;
  int outside;
  pari_CATCH(CATCH_ALL){//Catching if U is trivial OR z=0
    avma=top;
    pari_CATCH_reset();
    return -1;
  }
  pari_TRY{
    GEN ang=shiftangle(garg(z, prec), gel(gel(U, 4), 1), tol, prec);//Shifting to base angle
    long ind=gen_search(gel(U, 4), ang, 1, NULL, &gcmp_strict);//Index to put z. We ONLY need to search for this cicle.
    if(ind==lg(gel(U, 1))) ind=1;//Insert at the end means the first circle.
    GEN circle=gel(gel(U, 2), ind);
    if(gequal0(circle)){pari_CATCH_reset();avma=top;return -1;}//Intersects with the edge of the unit disc.
    outside=tolcmp(gel(circle, 2), gabs(gsub(z, gel(circle, 1)), prec), tol, prec);//Are we outside?
    if(outside==0) outside=-1;
    else if(outside==1) outside=ind;
  }
  pari_ENDCATCH
  avma=top;
  return outside;//There is no tolerance issues with our search for ind; they are taken care of by the tolerance check with inside (the only possible issues occur if z and v[ind] or v[ind-1] are equal up to tolerance, but of course that is solved by the one tolcmp).
}

//Given g in PSL(2, R) and p in the upper half plane, this returns the image of g in PSU(1, 1) via phi(z)=(z-p)/(z-conj(p)).
GEN psltopsu(GEN g, GEN p){
  pari_sp top=avma;
  GEN M=psltopsu_transmats(p);
  return gerepileupto(top, psltopsu_mats(g, M));
}

//Given g in PSL(2, R), M=[m1, m2, p] with m1=1/(p-conj(p))[1,-p;1,-conj(p)], m2=[conj(p), -p;1, -1], this returns m1*g*m2.
GEN psltopsu_mats(GEN g, GEN M){
  pari_sp top=avma;
  return gerepileupto(top, gmul(gel(M, 1), gmul(g, gel(M, 2))));
}

//Returns [m1, m2, p]: matrices so that if g is in PSL(2, R), then the image of g in PSU(1, 1) corresponding to phi(z)=(z-p)/(z-conj(p)) is m1*g*m2.
GEN psltopsu_transmats(GEN p){
  pari_sp top=avma;
  GEN pneg=gneg(p), pconj=conj_i(p);
  GEN m1=gdiv(mkmat22(gen_1, pneg, gen_1, gneg(pconj)), gsub(p, pconj));//m1=1/(p-conj(p))[1,-p;1,-conj(p)]
  GEN m2=mkmat22(pconj, pneg, gen_1, gen_m1);//m2=[conj(p), -p;1, -1]
  GEN ret=cgetg(4, t_VEC);
  gel(ret, 1)=gcopy(m1);
  gel(ret, 2)=gcopy(m2);
  gel(ret, 3)=gcopy(p);
  return gerepileupto(top, ret);
}

//Returns the roots of the hyperbolic matrix M in PSL(2, R) in order.
static GEN psl_roots(GEN M, GEN tol, long prec){
  pari_sp top = avma;
  GEN trace=gadd(gcoeff(M, 1, 1), gcoeff(M, 2, 2));
  int sgn=tolcmp(trace, gen_0, tol, prec);
  if(sgn==0) pari_err_TYPE("Please enter a hyperbolic matrix", M);
  GEN a, b, c;
  if(sgn==1){//Positive trace, correct order
    a=gcoeff(M, 2, 1);
    b=gsub(gcoeff(M, 2, 2), gcoeff(M, 1, 1));
    c=gneg(gcoeff(M, 1, 2));
  }
  else{
    a=gneg(gcoeff(M, 2, 1));
    b=gsub(gcoeff(M, 1, 1), gcoeff(M, 2, 2));
    c=gcoeff(M, 1, 2);
  }//[a',b';c',d'] -> c'x^2+(d'-a')x-b'=0, but for the roots to be in proper order, we need the trace to be positive.
  if(toleq(a, gen_0, tol, prec)){//a=0, roots are oo and -c/b (b!=0 else M=[+/-1, x;0;+/-1], not hyperbolic.
    GEN rnum=gneg(c);
    int bsgn=tolcmp(b, gen_0, tol, prec);
    GEN ret=cgetg(3, t_VEC);
    if(bsgn==1){//b>0, first root is finite
      gel(ret, 1)=gdiv(rnum, b);
      gel(ret, 2)=mkoo();
    }
    else if(bsgn==-1){//b<0, first root is oo
      gel(ret, 1)=mkoo();
      gel(ret, 2)=gdiv(rnum, b);
    }
    else pari_err_TYPE("Please enter a hyperbolic matrix", M);
    return gerepileupto(top, ret);
  }
  //Now both roots are finite.
  GEN twoa=gmulsg(2, a);//2a
  GEN rtD=gsqrt(gsub(gsqr(b), gmulsg(2, gmul(twoa, c))), prec);//b^2-4ac
  GEN mb=gneg(b);
  GEN r1num=gadd(mb, rtD);//first root is (-b+sqrt(D))/2a
  GEN r2num=gsub(mb, rtD);//Second root is (-b-sqrt(D))/2a
  GEN ret=cgetg(3, t_VEC);
  gel(ret, 1)=gdiv(r1num, twoa);
  gel(ret, 2)=gdiv(r2num, twoa);
  return gerepileupto(top, ret);
}

//Returns a random point z in the unit disc, uniform inside the ball of radius R. See page 19 of Page (before section 2.5).
GEN randompoint_ud(GEN R, long prec){
  pari_sp top=avma;
  GEN arg=gmul(randomr(prec), Pi2n(1, prec));//Random angle
  GEN zbound=expIr(arg);//The boundary point. Now we need to scale by a random hyperbolic distance in [0, R]
  //a(r)=Area of hyperbolic disc radius r=4*Pi*sinh^2(r/2).
  GEN dist=gmul(gsqr(gsinh(gdivgs(R, 2), prec)), randomr(prec));//A random element in [0, a(R)/4Pi].
  GEN r=gmulsg(2, gasinh(gsqrt(dist, prec), prec));//The radius
  GEN e2r=gexp(r, prec);
  return gerepileupto(top, gmul(zbound, gdiv(gsubgs(e2r, 1), gaddgs(e2r, 1))));
}

//Returns a random point z in the unit disc, uniform inside the ball of radius R. See page 19 of Page (before section 2.5).
GEN randompoint_udarc(GEN R, GEN ang1, GEN ang2, long prec){
  pari_sp top=avma;
  GEN arg=gadd(ang1, gmul(randomr(prec), gsub(ang2, ang1)));//Random angle in [ang1, ang2]
  GEN zbound=expIr(arg);//The boundary point. Now we need to scale by a random hyperbolic distance in [0, R]
  //a(r)=Area of hyperbolic disc radius r=4*Pi*sinh^2(r/2).
  GEN dist=gmul(gsqr(gsinh(gdivgs(R, 2), prec)), randomr(prec));//A random element in [0, a(R)/4Pi].
  GEN r=gmulsg(2, gasinh(gsqrt(dist, prec), prec));//The radius
  GEN e2r=gexp(r, prec);
  return gerepileupto(top, gmul(zbound, gdiv(gsubgs(e2r, 1), gaddgs(e2r, 1))));
}

//Reduces g with respect to z as in reduceelt_givenpsu, but does so much more efficiently using the normalized boundary provided.
GEN reduceelt_givennormbound(GEN U, GEN g, GEN z, GEN gamid, GEN *data, GEN (*gamtopsl)(GEN *, GEN, long), GEN (*eltmul)(GEN *, GEN, GEN), GEN tol, long prec){
  pari_sp top=avma;
  GEN delta=gamid;
  llist *decomp=NULL;
  GEN gmat=psltopsu_mats(gamtopsl(data, g, prec), gel(U, 8));
  z=mat_eval(gmat, z);//gz, the real start point
  long count=0, outside;
  for(;;){
    outside=normalizedboundary_outside(U, z, tol, prec);
    if(outside==-1) break;//Done:Either reached inside or the boundary.
    z=mat_eval(gel(gel(U, 5), outside), z);//Update z
    delta=eltmul(data, gel(gel(U, 1), outside), delta);//update delta
    llist_putstart(&decomp, outside);//add outside to the list
    count++;
  }
  GEN ret=cgetg(4, t_VEC);
  gel(ret, 1)=eltmul(data, delta, g);//gbar=delta*g
  gel(ret, 2)=gcopy(delta);
  gel(ret, 3)=llist_tovecsmall(decomp, count, 1);
  return gerepileupto(top, ret);
}

//Algorithm 4.3 of Voight. Inputs G, a finite subset of Gamma, corresponding to Gmats in PSU(1, 1), g (->gmat) an element of Gamma, z in the unit disc. This G-reduces g, i.e. translating gz to the exterior domain of G. Returns [gbar, delta, decomp], where gbar=delta*g. gbar is (G, z)-reduced, delta is in <G>, and delta=G[i1]*G[i2]*...*G[in] with decomp=[i1, i2, ..., in] (vecsmall).
GEN reduceelt_givenpsu(GEN G, GEN Gmats, GEN g, GEN gmat, GEN z, GEN gamid, GEN *data, GEN (*eltmul)(GEN *, GEN, GEN), GEN tol, long prec){
  pari_sp top=avma;
  GEN delta=gamid;
  llist *decomp=NULL;
  z=mat_eval(gmat, z);//gz, the real start point
  GEN mindist, curdist=hdist_ud(z, gen_0, prec), znew, zmin, dist;
  long ind, n=lg(G), count=0;
  for(;;){
    zmin=mat_eval(gel(Gmats, 1), z);
    mindist=hdist_ud(zmin, gen_0, prec);
    if(typ(mindist)==t_COMPLEX) mindist=mkoo();//Rounding error forced point outside the unit disc.
    ind=1;
    for(long i=2;i<n;i++){
      if(gequal0(gel(Gmats, i))) continue;//oo side, ignore
      znew=mat_eval(gel(Gmats, i), z);
      dist=hdist_ud(znew, gen_0, prec);//Distance to g_i*z
      if(typ(dist)==t_COMPLEX) continue;//Rounding error forced us onto point is too close to being outside the unit disc.
      if(tolcmp(dist, mindist, tol, prec)==-1){//Strictly smaller distance
        zmin=znew;
        mindist=dist;
        ind=i;
      }
    }
    if(tolcmp(mindist, curdist, tol, prec)!=-1) break;//Done
    count++;
    z=zmin;
    curdist=mindist;
    delta=eltmul(data, gel(G, ind), delta);
    llist_putstart(&decomp, ind);
  }
  GEN ret=cgetg(4, t_VEC);
  gel(ret, 1)=eltmul(data, delta, g);//gbar=delta*g
  gel(ret, 2)=gcopy(delta);
  gel(ret, 3)=llist_tovecsmall(decomp, count, 1);
  return gerepileupto(top, ret);
}

//Reduces z to the interior of U (Almost identical to reduceelt_givennormbound). Returns [g, z'], where g is the transition element and z' is the new point.
GEN reducepoint(GEN U, GEN z, GEN gamid, GEN *data, GEN (*eltmul)(GEN *, GEN, GEN), GEN tol, long prec){
  pari_sp top=avma;
  GEN g=gamid;
  long outside;
  for(;;){
    outside=normalizedboundary_outside(U, z, tol, prec);
    if(outside==-1) break;//Done:Either reached inside or the boundary.
    z=mat_eval(gel(gel(U, 5), outside), z);//Update z
    g=eltmul(data, gel(gel(U, 1), outside), g);//update g
  }
  GEN ret=cgetg(3, t_VEC);
  gel(ret, 1)=gcopy(g);
  gel(ret, 2)=gcopy(z);
  return gerepileupto(top, ret);
}

//Returns the root geodesic in the unit disc corresponding to M in PSL(2, R) and p the reference point to mapping the upper half plane to the unit disc (mats=[m1, m2, p], with m1 being the mapping).
GEN rootgeodesic_ud(GEN M, GEN mats, GEN tol, long prec){
  pari_sp top=avma;
  GEN geod=rootgeodesic_uhp(M, tol, prec);
  return gerepileupto(top, mobius(gel(mats, 1), geod, tol, prec));
}

//Returns the upper half plane root geodesic of the hyperbolic element M in PSL(2, R)
GEN rootgeodesic_uhp(GEN M, GEN tol, long prec){
  pari_sp top=avma;
  GEN rts=psl_roots(M, tol, prec), arc;
  if(typ(gel(rts, 1))==t_INFINITY){//First root infinite; we have a vertical segment from rts[2] to oo
    arc=cgetg(ARCLEN, t_VEC);
    gel(arc, 1)=mkoo();
    gel(arc, 2)=gcopy(gel(rts, 2));
    gel(arc, 3)=gcopy(gel(rts, 2));//Starts at finite (second) root
    gel(arc, 4)=mkoo();//Ends at oo
    gel(arc, 5)=gen_0;
    gel(arc, 6)=gen_1;//Vertically up
    gel(arc, 7)=gen_0;//No dir since endpt is oo
    gel(arc, 8)=gen_1;//Segment
    return gerepileupto(top, arc);
  }
  else if(typ(gel(rts, 2))==t_INFINITY){//Second root infinite; we have a vertical segment fr
    arc=cgetg(ARCLEN, t_VEC);
    gel(arc, 1)=mkoo();
    gel(arc, 2)=gcopy(gel(rts, 1));
    gel(arc, 3)=mkoo();//Starts at oo
    gel(arc, 4)=gcopy(gel(rts, 1));//Ends at finite (second) root
    gel(arc, 5)=gen_0;
    gel(arc, 6)=gen_m1;//Vertically down
    gel(arc, 7)=gen_0;//No dir since endpt is oo
    gel(arc, 8)=gen_1;//Segment
    return gerepileupto(top, arc);
  }
  //Now both roots are finite, so we have a nice circle.
  GEN centre=gdivgs(gadd(gel(rts, 1), gel(rts, 2)), 2);
  int firstbigger=gcmp(gel(rts, 1), gel(rts, 2));//=1 if the first root is bigger, and -1 if not (no tolerance needed).
  if(firstbigger==1){
    GEN radius=gsub(gel(rts, 1), centre);
    GEN c=mkvec3(centre, radius, gen_0);//The full circle
    arc=arc_init(c, gel(rts, 1), gel(rts, 2), -1, prec);//Arc points from second to first root.
  }
  else{
    GEN radius=gabs(gsub(gel(rts, 2), centre), prec);
    GEN c=mkvec3(centre, radius, gen_0);//The full circle
    arc=arc_init(c, gel(rts, 2), gel(rts, 1), 1, prec);//Arc points from second to first root.
  }
  return gerepileupto(top, arc);
}


//FUNDAMENTAL DOMAIN OTHER COMPUTATIONS


//Returns the set of minimal cycles of the side pairing pair. A cycle is a vecsmall [i1,i2,...,in] so that the cycle is v_i1, v_i2, ..., v_in. A cycle [-i] means that the "vertex" on side i is is a one element cycle (happens when a side is fixed).
static GEN minimalcycles(GEN pair){
  pari_sp top=avma;
  long np1=lg(pair), n=np1-1, vleft=np1;//Number of sides/vertices (not counting vertices that occur on the middle of a side).
  GEN vind=cgetg(np1, t_VECSMALL);
  for(long i=1;i<np1;i++) vind[i]=1;//Tracking if the vertices have run out or not
  GEN cycles=vectrunc_init(2*np1), cyc;//Max number of cycles, since each side could have a middle vertex. In reality the number is probably much smaller, but this is safe.
  long startind=1, ind;
  for(long i=1;i<np1;i++){//We sort the fixed sides first, as later on we would miss the ones that get removed before checking.
    if(pair[i]==i){//Side fixed!
	  vectrunc_append(cycles, mkvecsmall(-i));//Middle of the side is fixed.
	} 
  }
  do{
	cyc=vecsmalltrunc_init(vleft);
	vecsmalltrunc_append(cyc, startind);//Starting the cycle.
	vind[startind]=0;
	vleft--;
	ind=smodss(pair[startind]-2, n)+1;//Hit it with the side pairing and subtract 1 to reach the paired vertex.
	while(ind!=startind){//Move along the cycle.
	  vind[ind]=0;
	  vleft--;//One less vertex
	  vecsmalltrunc_append(cyc, ind);//Append it
	  ind=smodss(pair[ind]-2, n)+1;//Update
	}
	vectrunc_append(cycles, cyc);//New cycle.
    while(startind<np1){//Finding the next vertex we haven't eliminated.
	  startind++;
	  if(vind[startind]==1) break;
	}
  }
  while(startind<np1);
  return gerepilecopy(top, cycles);
}

//Returns [cycles, types], where cycles[i] has type types[i]. Type 0=parabolic, 1=accidental, m>=2=elliptic of order m.
GEN minimalcycles_bytype(GEN U, GEN gamid, GEN *data, GEN (*eltmul)(GEN *, GEN, GEN), GEN (*elttrace)(GEN *, GEN), int (*istriv)(GEN *, GEN)){
  pari_sp top=avma;
  GEN G=gel(U, 1);
  GEN cycles=minimalcycles(gel(U, 7)), cyc, g, trd;//The minimal cycles.
  long ncyc=lg(cycles);
  GEN types=cgetg(ncyc, t_VECSMALL);//The types
  for(long i=1;i<ncyc;i++){
	cyc=gel(cycles, i);
	if(lg(cyc)==2 && cyc[1]<0) g=gel(G, -cyc[1]);//Minimal cycle that was the middle of a side.
	else{
	  g=gamid;
	  for(long j=1;j<lg(cyc);j++) g=eltmul(data, gel(G, cyc[j]), g);//Multiply on the left by G[cyc[j]]
	}
	if(istriv(data, g)){types[i]=1;continue;}//Accidental cycle, continue on.
	trd=elttrace(data, g);//The trace
	if(gequal(trd, gen_2) || gequal(trd, gen_m2)){types[i]=0;continue;}//Parabolic cycle.
    long ord=1;
	GEN gpow=g;
	do{//Finding the order of g
	  ord++;
	  gpow=eltmul(data, g, gpow);
	}
	while(!istriv(data, gpow));
	types[i]=ord;
  }
  GEN ordering=vecsmall_indexsort(types);
  return gerepilecopy(top, mkvec2(vecpermute(cycles, ordering), vecsmallpermute(types, ordering)));//The return, [cycles, types]
}

//Returns the vecsmall of indices of the infinite sides of U.
GEN normalizedboundary_oosides(GEN U){
  long n=lg(gel(U, 1));
  GEN sides=vecsmalltrunc_init(n);
  for(long i=1;i<n;i++) if(gequal0(gel(gel(U, 2), i))) vecsmalltrunc_append(sides, i);//Append the sides
  return sides;
}

//Returns the group presentation of the fundamental domain U. The return is a vector, where the 1st element is the list of indices of the generators, 2nd element is the vector of relations, whose ith element is a relation of the form [indices, powers], where indices and powers are vecsmall's. If indices=[i1,i2,...,ik] and powers=[p1,p2,...,pk], then this corresponds to g_{i1}^p1*...*g_{ik}^{pk}=1.
GEN presentation(GEN U, GEN gamid, GEN *data, GEN (*eltmul)(GEN *, GEN, GEN), GEN (*elttrace)(GEN *, GEN), int (*istriv)(GEN *, GEN)){
  pari_sp top=avma;
  GEN mcyc=minimalcycles_bytype(U, gamid, data, eltmul, elttrace, istriv);//Minimal cycles by type.
  GEN cyc=gel(mcyc, 1), cyctype=gel(mcyc, 2);
  long lgelts=lg(gel(U, 1)), lgcyc=lg(cyc), ngens=0;
  GEN H=cgetg(lgelts, t_VECSMALL);
  for(long i=1;i<lgelts;i++){//H[i]=1 if g=g^(-1), and for exactly one of [g,g^(-1)] for all other g.
	if(gel(U, 7)[i]>=i){H[i]=1;ngens++;}
	else H[i]=0;
  }
  long ind=1, k;
  while(ind<lgcyc && cyctype[ind]==0) ind++;//Find the first accidental cycle.
  long ellind=ind;
  while(ellind<lgcyc && cyctype[ellind]==1) ellind++;//Find the first elliptic cycle.
  long naccident=ellind-ind;//How many accidental cycles
  GEN r=gen_0;//The accidental relation, if it exists.
  if(naccident>0){//We have some accidental cycles!
	r=gcopy(gel(cyc, ind));
	r=vecsmall_reverse(r);//The relation is r backwards.
	if(naccident>1){//More than one relation.
	  ngens=ngens-naccident+1;//Updating the number of generators.
	  long lastrel=ellind-1;//The last relation we consider in a cycle. We swap this with the relation we find every step and decrease it, so we only need to consider relations from ind+1 to lastrel at each step.
	  long indrep=1, torep;//Stores the index replaced in r. On the next pass, we may as well start from there, as we have already checked the previous indices for replacement.
	  GEN repind=gen_0, cycle, newr;//Stores [porm, j, l, m], where term j in relation r is replaced by using term m of relation l. If porm=-1, we need to replace the inverse, otherwise we do not.
	  for(long i=1;i<naccident;i++){//Each step we solve the relation.
	    for(long j=indrep;j<lg(r);j++){//Trying to replace index j
		  torep=r[j];
		  if(torep<0) torep=-torep;//In case power is -1.
		  for(long l=ind+1;l<=lastrel;l++){//Looking at cycle l
			for(long m=1;m<lg(gel(cyc, l));m++){//element m of cycle l
			  if(gel(cyc, l)[m]==torep){//Replace it!
			    H[torep]=0;
				H[gel(U, 7)[torep]]=0;//Make sure it has been deleted from H.
				repind=mkvecsmall4(1, j, l, m);
				if(r[j]<0) repind[1]=-1;//Actually the inverse.
				l=lastrel+1;//Break loop
				j=lg(r);//Break loop
				break;//Break loop
			  }
			}
		  }
		  if(gel(U, 7)[torep]==torep) continue;//g=g^(-1), so no need to search for the inverse.
		  torep=gel(U, 7)[torep];//Now we try to replace the inverse
		  for(long l=ind+1;l<=lastrel;l++){//Looking at cycle l
			for(long m=1;m<lg(gel(cyc, l));m++){//element m of cycle l
			  if(gel(cyc, l)[m]==torep){//Replace it!
			    H[torep]=0;
				H[gel(U, 7)[torep]]=0;//Make sure it has been deleted from H.
				repind=mkvecsmall4(1, j, l, m);
				if(r[j]>0) repind[1]=-1;//Actually the inverse.
				l=lastrel+1;//Break loop
				j=lg(r);//Break loop
				break;//Break loop
			  }
			}
		  }
		}
		//Now, repind gives us all the information we need to do the replacement.
		cycle=gel(cyc, repind[3]);
		newr=cgetg(lg(r)+lg(cycle)-3, t_VECSMALL);
		for(long j=1;j<repind[2];j++) newr[j]=r[j];//First part is the same.
		k=repind[2];//Keeps track of the index of newr that we are working on.
		if(repind[1]==1){//Replace it normally. Cycle [a1,...,an] -> an*...*a1=1 -> a_j=a_{j+1}^(-1)*...*a_n^(-1)*a_1^(-1)*...*a_{j-1}^(-1).
		  for(long j=repind[4]+1;j<lg(cycle);j++){newr[k]=-cycle[j];k++;}
		  for(long j=1;j<repind[4];j++){newr[k]=-cycle[j];k++;}
		}
		else{//Replace its inverse. Cycle [a1,...,an] -> an*...*a1=1 -> a_j^(-1)=a_{j-1}*...*a_1*a_n*...*a_{j+1}.
		  for(long j=repind[4]-1;j>0;j--){newr[k]=cycle[j];k++;}
		  for(long j=lg(cycle)-1;j>repind[4];j--){newr[k]=cycle[j];k++;}
		}
		for(long j=repind[2]+1;j<lg(r);j++){newr[k]=r[j];k++;}//The final part is the same.
		r=newr;
		indrep=repind[2];//The index we replaced.
		gel(cyc, repind[3])=gel(cyc, lastrel);//Replace the relation we replaced with the last one.
		lastrel--;//One less relation.
	  }
	}
  }
  //Now we have to finalize things. First, we prepare for replacement of g and g^(-1)
  long oppind;
  GEN indused=cgetg(ngens+1, t_VECSMALL);//Which indices were used.
  k=0;
  for(long i=1;i<lgelts;i++){
	if(H[i]<=0) continue;//Nothing to do
	k++;
	indused[k]=i;
    H[i]=i;
	oppind=gel(U, 7)[i];
	if(oppind==i) continue;//g=g^(-1), continue on
	H[oppind]=-i;//The inverse of g.
  }//Now when we see an element with index i, we should replace it with H[i] if H[i]>0, and (-H[i])^(-1) if H[i]<0. Elements with H[i]=0 are guarenteed to not appear now.
  GEN relations;
  if(naccident==0) relations=cgetg(lgcyc-ellind+1, t_VEC);//We get 1 relation for every elliptic cycle, and an additional relation if we have >=1 accidental cycles.
  else relations=cgetg(lgcyc-ellind+2, t_VEC);
  k=1;
  long w;
  for(long i=ellind;i<lgcyc;i++){//Sorting the elliptic cycles first. We assume that they are all of length 1.
    w=gel(cyc, i)[1];
	if(w<0) w=-w;//Swap the sign back.
	w=H[w];//The real element.
	gel(relations, k)=cgetg(3, t_VEC);
	if(w>0) gmael(relations, k, 1)=mkvecsmall(w);
	else gmael(relations, k, 1)=mkvecsmall(-w);
	gmael(relations, k, 2)=mkvecsmall(cyctype[i]);
	k++;
  }
  if(naccident>0){
    gel(relations, k)=cgetg(3, t_VEC);//Now we have the last relation, coming from r.
    long lenr;
    gmael(relations, k, 1)=cgetg_copy(r, &lenr);
    gmael(relations, k, 2)=cgetg_copy(r, &lenr);
    for(long i=1;i<lenr;i++){
	  w=H[r[i]];
	  if(w>0){gmael(relations, k, 1)[i]=w;gmael(relations, k, 2)[i]=1;}
	  else{gmael(relations, k, 1)[i]=-w;gmael(relations, k, 2)[i]=-1;}
	}
  }
  return gerepilecopy(top, mkvec2(indused, relations));
}

//Finds the image of the root geodesic of g in the fundamental domain specified by U.
GEN rootgeodesic_fd(GEN U, GEN g, GEN gamid, GEN *data, GEN (*gamtopsl)(GEN *, GEN, long), GEN (*eltmul)(GEN *, GEN, GEN), GEN (*eltinv)(GEN *, GEN), GEN tol, long prec){
  pari_sp top=avma;
  GEN gpsl=gamtopsl(data, g, prec), mats=gel(U, 8);//The image in PSL
  GEN geod=rootgeodesic_ud(gpsl, mats, tol, prec);
  GEN z;
  if(gequal0(gel(geod, 8))) z=arc_midpoint(geod, gel(geod, 3), gel(geod, 4), tol, prec);//First move the midpoint to the fundamental domain;
  else z=gdivgs(gadd(gel(geod, 3), gel(geod, 4)), 2);//When geod is a line segment, do this instead (Q ram at 11, 13 and g=[11/2, 3/2, 0, 0] e.g.)
  GEN red=reducepoint(U, z, gamid, data, eltmul, tol, prec);
  g=eltmul(data, eltmul(data, gel(red, 1), g), eltinv(data, gel(red, 1)));//Conjugating g by gel(red, 1);
  gpsl=gamtopsl(data, g, prec);//The image in PSL of the new g
  geod=rootgeodesic_ud(gpsl, mats, tol, prec);//The new geodesic, which necessarily goes through z, and hence the interior.
  //Now we need to find the start vertex, and run along.
  GEN vbaseinfo=normalizedboundary_sideint(U, geod, 1, tol, prec);//Find the vertex nearest the start of the geodesic.
  GEN vbase=gel(vbaseinfo, 1);
  GEN vstart=vbase, vend, startcentre=gel(geod, 1), starttype=gel(geod, 8);
  glist *Gs=NULL;//Tracking the g's corresponding to the arcs
  glist *circs=NULL;//Tracking the circle arcs.
  llist *sides=NULL;//Tracking the sides hit.
  llist *othersides=NULL;//Tracking the sides we are leaving from.
  llist_putstart(&othersides, gel(vbaseinfo, 2)[1]);
  long count=0;
  for(;;){
    vend=normalizedboundary_sideint(U, geod, 0, tol, prec);
    glist_putstart(&Gs, g);
    if(gequal0(gel(geod, 8))){//Arc
      if(gequal(gel(geod, 7), gen_m1)) glist_putstart(&circs, arc_init(geod, gel(vend, 1), vstart, -1, prec));//geod travelling backwards
      else glist_putstart(&circs, arc_init(geod, vstart, gel(vend, 1), 1, prec));//geod travelling normally
    }
    else{//Segment
      gel(geod, 3)=vstart;//Start
      gel(geod, 4)=gel(vend, 1);//End
      glist_putstart(&circs, geod);
    }
    llist_putstart(&sides, gel(vend, 2)[1]);//Adding the side hit.
    count++;
    vstart=mat_eval(gel(gel(U, 5), gel(vend, 2)[1]), gel(vend, 1));//The new start vertex
    g=eltmul(data, eltmul(data, gel(gel(U, 1), gel(vend, 2)[1]), g), eltinv(data, gel(gel(U, 1), gel(vend, 2)[1])));//Conjugating g by the side hit
    gpsl=gamtopsl(data, g, prec);//The image in PSL of the new g
    geod=rootgeodesic_ud(gpsl, mats, tol, prec);//The new geodesic
    if(toleq(vbase, vstart, tol, prec)){
      if(gequal(starttype, gel(geod, 8)) && toleq(startcentre, gel(geod, 1), tol, prec)) break;//Done! We need this second check since an embedding CAN have a self-intersection on the boundary (e.g. [Q, order]=qa_init_2primes(2, 3), p=I/2, g=[3674890, -1623699, 463914, -1391742]). This is sufficient because the two arcs/segments share a point AND the centre/slope => unique.
    }
    llist_putstart(&othersides, gel(U, 7)[gel(vend, 2)[1]]);//Hit it with the side pairing.
  }
  GEN ret=cgetg(5, t_VEC);
  gel(ret, 1)=glist_togvec(Gs, count, -1);
  gel(ret, 2)=glist_togvec(circs, count, -1);
  gel(ret, 3)=llist_tovecsmall(sides, count, -1);
  gel(ret, 4)=llist_tovecsmall(othersides, count, -1);
  return gerepileupto(top, ret);
}

//Computes the signature of the fundamental domain U. The return in [g, V, s], where g is the genus, V=[m1,m2,...,mt] (vecsmall) are the orders of the elliptic cycles (all >=2), and s is the number of parabolic cycles. The signature is normally written as (g;m1,m2,...,mt;s).
GEN signature(GEN U, GEN gamid, GEN *data, GEN (*eltmul)(GEN *, GEN, GEN), GEN (*elttrace)(GEN *, GEN), int (*istriv)(GEN *, GEN)){
  pari_sp top=avma;
  GEN mcyc=minimalcycles_bytype(U, gamid, data, eltmul, elttrace, istriv);//The minimal cycles and their types.
  long nfixed=0, lgcyc=lg(gel(mcyc, 1));//The number of fixed sides, and number of cycles+1
  for(long i=1;i<lgcyc;i++){if(gmael(mcyc, 1, i)[1]<0) nfixed++;}
  long genus=((lg(gel(U, 1))-1+nfixed)/2-lgcyc+2)/2;//2*g+t+s+1_{t+s>0}=min number of generators. Initially, we have (n+k)/2, where k is the number of sides fixed by the element (nfixed) and n is the number of sides of the fdom (b/c we take one of g and g^(-1) every time). Then each accidental cycle AFTER THE FIRST removes exactly one generator. This gives the formula (if there are no accidental cycles then we are off by 1/2, but okay as / rounds down in C. We are actually overcounting by 1 if t+s=0 and there are no accidental cycles, but this is impossible as there is >=1 cycle.
  long s, firstell=lgcyc;//s counts the number of parabolic cycles, and firstell is the first index of an elliptic cycle.
  int foundlastpar=0;
  for(long i=1;i<lgcyc;i++){
	if(!foundlastpar){
	  if(gel(mcyc, 2)[i]==0) continue;
	  s=i-1;//Found the last parabolic.
	  foundlastpar=1;
	}
	if(gel(mcyc, 2)[i]==1){continue;}
	firstell=i;//Found the last accidental.
	break;
  }
  long lgell=lgcyc-firstell+1;
  GEN rvec=cgetg(4, t_VEC);//[g, V, s]
  gel(rvec, 1)=stoi(genus);
  gel(rvec, 2)=cgetg(lgell, t_VECSMALL);
  for(long i=1;i<lgell;i++){gel(rvec, 2)[i]=gel(mcyc, 2)[firstell];firstell++;}
  gel(rvec, 3)=stoi(s);
  return gerepileupto(top, rvec);
}


//GEOMETRIC HELPER METHODS


//Returns the angle ang-bot in the range [0, 2*Pi)
static GEN anglediff(GEN ang, GEN bot, GEN tol, long prec){
  pari_sp top=avma;
  GEN twopi=Pi2n(1, prec);
  GEN angdiff=gmod(gsub(ang, bot), twopi);
  if(toleq(angdiff, twopi, tol, prec) || toleq(angdiff, gen_0, tol, prec)){avma=top;return gen_0;}
  return gerepileupto(top, angdiff);
}

//Arctan of x, where oo is allowed as input (returning Pi/2)
static GEN atanoo(GEN x, long prec){
  pari_sp top=avma;
  if(typ(x)==t_INFINITY) return gerepileupto(top, gdivgs(mppi(prec), 2));
  return gatan(x, prec);
}

//Returns gcmp(x, y), except if x==y, returns -1. Useful for gen_search when you ALWAYS want to return the index where it should be inserted/is
static int gcmp_strict(void *data, GEN x, GEN y){
  int c=gcmp(x, y);
  if(c==0) c=-1;
  return c;
}

//Returns 0 if c is a circle, 1 if c is a line, 2 if c is a circle arc, 3 if c is line segment, and -1 if none of the above
static int geom_check(GEN c){
  if(typ(c)!=t_VEC) return -1;
  int len=lg(c);
  if(len==CIRCLEN){//Circle or line
    if(gequal0(gel(c, 3))) return 0;//Circle
    if(gequal(gel(c, 3), gen_1)) return 1;//Line
    return -1;//Neither
  }
  else if(len==ARCLEN){
    if(gequal0(gel(c, 8))) return 2;//Circle arc
    if(gequal(gel(c, 8), gen_1)) return 3;//Line segment
    return -1;//Neither
  }
  return -1;//Not a circle/line
}

//Returns the default tolerance given the precision.
GEN deftol(long prec){
  pari_sp top=avma;
  return gerepileupto(top, gtofp(powis(gen_2, 32*(2-prec)), prec));
}

//Shifts the given angle ang by multiples of 2*Pi into the range [bot, bot+2*Pi).
static GEN shiftangle(GEN ang, GEN bot, GEN tol, long prec){
  pari_sp top=avma;
  GEN diff=anglediff(ang, bot, tol, prec);
  if(gequal0(diff)){avma=top;return gcopy(bot);}
  return gerepileupto(top, gadd(bot, diff));
}

//Returns -1 if x<y, 0 if x==y, 1 if x>y (x, y are t_REAL). Accounts for the tolerance, so will deem x==y if they are equal up to tol AND at least one is inexact
static long tolcmp(GEN x, GEN y, GEN tol, long prec){
  if(typ(x)==t_INFINITY || typ(y)==t_INFINITY) return gcmp(x, y);//No precision concerns
  pari_sp top=avma;
  GEN d=gsub(x, y);
  if(precision(d)==0){
    long ret=gsigne(d);
    avma=top;
    return ret;
  }//Exact objects
  if(gcmp(gabs(d, prec), tol)==-1){avma=top;return 0;}//Within tolerance
  long ret=gsigne(d);
  avma=top;return ret;
}

//Data points to [tol, vecsmall(prec)]. Used to sort/search a list with tolerance.
static int tolcmp_sort(void *data, GEN x, GEN y){
  return tolcmp(x, y, gel(*(GEN*)data, 1), gel(*(GEN*)data, 2)[1]);
}

//Returns 1 if x==y, and 0 if x!=y. If x or y is not a precise objects (e.g. t_REAL), will return 1 if they are equal up to the tolerance tol.
static int toleq(GEN x, GEN y, GEN tol, long prec){
  if(typ(x)==t_INFINITY || typ(y)==t_INFINITY || gequal0(tol)) return gequal(x, y);//No precision concerns
  pari_sp top=avma;
  GEN d=gsub(x, y);
  if(gequal0(d)){avma=top;return 1;}//Deemed equal already
  if(precision(d)==0){avma=top;return 0;}//Exact objects
  if(gcmp(gabs(d, prec), tol)==-1){avma=top;return 1;}//Within tolerance
  avma=top;return 0;
}



//SECTION 3: QUATERNIONIC METHODS



/*
To compute the fundamental domain, we store the quaternion algebra as [A, ramdat, varnos, roots], where A=the algebra, ramdat=finite ramified places of A, varnos are the variable numbers for K, L (L the splitting field of K), and roots are approximations to the defining variables for K, L that correspond to the split real embedding. Such an entry is called a "qalg" and uses the variable name Q, so any methods with the name "qalg" expect this as input (whereas "alg" expects an algebra that was output by alginit, and this uses the variable name A).
*/


//QUATERNION ALGEBRA METHODS


//Returns the absolute reduced norm with respect to z1 and z2, i.e. the quadratic from Q_{z_1, z_2}(g) for g in the algebra. If the output is q, then g (written in basis form) has value g~*q*g.
GEN algabsrednorm(GEN A, GEN p, GEN z1, GEN z2, long prec){
  pari_sp top=avma;
  GEN Q=qalg_fdominitialize(A, prec);
  GEN mats=psltopsu_transmats(p);
  return gerepileupto(top, qalg_absrednormqf(Q, mats, z1, z2, gen_0, prec));
}

//Initializes and checks the inputs, and computes the fundamental domain. Can supply constants as 0 or [C, R, passes, type]. Any entry that is 0 is auto-set.
GEN algfdom(GEN A, GEN p, int dispprogress, int dumppartial, GEN partialset, GEN constants, long prec){
  pari_sp top=avma;
  GEN tol=deftol(prec);
  GEN Q=qalg_fdominitialize(A, prec), newA=A, U;
  if(typ(constants)!=t_VEC || lg(constants)<5) constants=zerovec(4);
  if(gequal0(p)){p=cgetg(3, t_COMPLEX);gel(p, 1)=gen_0;gel(p, 2)=ghalf;}
  long newprec=prec;
  unsigned int precinc=0;
  pari_CATCH(e_PREC){
	pari_warn(warner, "Increasing precision");
	newA=algmoreprec(newA, 1, newprec);
	newprec++;
	precinc++;
	tol=deftol(newprec);
	Q=qalg_fdominitialize(newA, newprec);
  }
  pari_RETRY{
    U=qalg_fdom(Q, p, dispprogress, dumppartial, partialset, gel(constants, 1), gel(constants, 2), gel(constants, 3), itos(gel(constants, 4)), tol, newprec);
  }pari_ENDCATCH
  if(precinc) pari_warn(warner, "Precision increased %d times to %d (the number of times it increased may be wrong, but the final precision should be correct). Please recompile your number field and algebra with precision \\p%Pd", precinc, newprec, precision00(U, NULL));
  return gerepileupto(top, U);
}

//Generate the optimal C value
GEN algfdom_bestC(GEN A, long prec){
  pari_sp top=avma;
  GEN nf=alg_get_center(A);//Field
  long n=nf_get_degree(nf);//Degree of field
  GEN algdisc=algnormdisc(A);//Norm of disc
  GEN discpart=gmul(nf_get_disc(nf), gsqrt(algdisc, prec));//disc(F)*sqrt(algdisc)
  GEN discpartroot=gpow(discpart, gdivgs(gen_1, n), prec);//discpart^(1/n)=disc(F)^(1/n)*algdisc^(1/2n)
  GEN npart;
  if(n==1) npart=dbltor(2.8304840896);
  else if(n==2) npart=dbltor(0.9331764427);
  else if(n==3) npart=dbltor(0.9097513831);
  else if(n==4) npart=dbltor(0.9734563346);
  else if(n==5) npart=dbltor(1.0195386113);
  else if(n==6) npart=dbltor(1.0184814342);
  else if(n==7) npart=dbltor(0.9942555240);
  else if(n==8) npart=dbltor(0.9644002039);
  else npart=gen_1;
  GEN best=gerepileupto(top, gmul(npart, discpartroot));//npart*disc(F)^(1/n)*N_F/Q(algebra disc)^(1/2n)
  if(gcmpgs(best, n)<=0) best=gerepileupto(top, gaddsg(n, gen_2));//Make sure best>n. If it is not, then we just add 2 (I doubt this will ever occur, but maybe in a super edge case).
  return best;
}

//Returns the area of the fundamental domain of the order stored in A.
GEN algfdomarea(GEN A, int lessprec, long prec){
  pari_sp top=avma;
  GEN Q=qalg_fdominitialize(A, prec);
  long lp = lessprec? 3:prec;
  return gerepileupto(top, qalg_fdomarea(Q, lp, prec));
}

//Returns the minimal cycles in the fundamental domain U of the algebra A.
GEN algfdomminimalcycles(GEN A, GEN U, long prec){
  pari_sp top=avma;
  GEN id=gel(alg_get_basis(A), 1);//The identity
  GEN Q=qalg_fdominitialize(A, prec);
  return gerepileupto(top, minimalcycles_bytype(U, id, &Q, &qalg_fdommul, &qalg_fdomtrace, &qalg_istriv));
}

//Returns the presentation of the algebra A, obtained from the fundamental domain U.
GEN algfdompresentation(GEN A, GEN U, long prec){
  pari_sp top=avma;
  GEN id=gel(alg_get_basis(A), 1);//The identity
  GEN Q=qalg_fdominitialize(A, prec);
  return gerepileupto(top, presentation(U, id, &Q, &qalg_fdommul, &qalg_fdomtrace, &qalg_istriv));
}

//Reduces the norm 1 element x with respect to the fundamental domain fdom and the point z (default z=0)
GEN algfdomreduce(GEN A, GEN U, GEN g, GEN z, long prec){
  pari_sp top=avma;
  GEN tol=deftol(prec);
  GEN id=gel(alg_get_basis(A), 1);//The identity
  GEN Q=qalg_fdominitialize(A, prec);
  return gerepileupto(top, reduceelt_givennormbound(U, g, z, id, &Q, &qalg_fdomm2rembed, &qalg_fdommul, tol, prec));
}

//Returns the root geodesic of g translated to the fundamental domain U. Output is [elements, arcs, vecsmall of sides hit, vecsmall of sides left from].
GEN algfdomrootgeodesic(GEN A, GEN U, GEN g, long prec){
  pari_sp top=avma;
  GEN tol=deftol(prec);
  GEN id=gel(alg_get_basis(A), 1);//The identity
  GEN Q=qalg_fdominitialize(A, prec);
  return gerepileupto(top, rootgeodesic_fd(U, g, id, &Q, &qalg_fdomm2rembed, &qalg_fdommul, &qalg_fdominv, tol, prec));
	
}

//Returns the signature of the quaternion algebra A with fundamental domain U.
GEN algfdomsignature(GEN A, GEN U, long prec){
  pari_sp top=avma;
  GEN id=gel(alg_get_basis(A), 1);//The identity
  GEN Q=qalg_fdominitialize(A, prec);
  return gerepileupto(top, signature(U, id, &Q, &qalg_fdommul, &qalg_fdomtrace, &qalg_istriv));
}

//Returns the same quaternion algebra, just with more precision. Assumes it currently has prec precision, then adds increment to the precision (or 1 if increment=0).
GEN algmoreprec(GEN A, long increment, long prec){
  pari_sp top=avma;
  if(increment<=0) increment=1;
  GEN nf=alg_get_center(A);
  GEN L=alg_get_splittingfield(A), pol=rnf_get_pol(L);//L=nf(sqrt(a))
  long varn=rnf_get_varn(L);
  GEN a=gneg(gsubst(pol, varn, gen_0));//Polynomial is of the form x^2-a, so plug in 0 and negate to get a
  GEN b=lift(alg_get_b(A));
  GEN aden=Q_denom(a), bden=Q_denom(b);//When initializing the algebra, PARI insists that a, b have no denominators
  if(!isint1(aden)){
	a=gmul(a, gsqr(aden));//We need to get rid of the denominator of a
	pari_warn(warner,"Changed the value of a, since it had a denominator.");
  }
  if(!isint1(bden)){
	b=gmul(b, gsqr(bden));//We need to get rid of the denominator of b
    pari_warn(warner,"Changed the value of a, since it had a denominator.");
  }
  long newprec=prec+increment;
  GEN newnf=nfinit(nf_get_pol(nf), newprec);
  return gerepileupto(top, alginit(newnf, mkvec2(a, b), -1, 1));
}

//Returns the normalized basis of the set of elements G
GEN algnormalizedbasis(GEN A, GEN G, GEN p, long prec){
  pari_sp top=avma;
  GEN mats=psltopsu_transmats(p);//Transition matrices
  GEN id=gel(alg_get_basis(A), 1);//The identity
  GEN Q=qalg_fdominitialize(A, prec);
  GEN tol=deftol(prec);
  return gerepileupto(top, normalizedbasis(G, gen_0, mats, id, &Q, &qalg_fdomm2rembed, &qalg_fdommul, &qalg_fdominv, &qalg_istriv, tol, prec));
}

//Returns the normalized boundary of the set of elements G
GEN algnormalizedboundary(GEN A, GEN G, GEN p, long prec){
  pari_sp top=avma;
  GEN mats=psltopsu_transmats(p);//Transition matrices
  GEN id=gel(alg_get_basis(A), 1);//The identity
  GEN Q=qalg_fdominitialize(A, prec);
  GEN tol=deftol(prec);
  return gerepileupto(top, normalizedboundary(G, mats, id, &Q, &qalg_fdomm2rembed, tol, prec));
}

//Returns the norm to Q of the discriminant of A
GEN algnormdisc(GEN A){
  pari_sp top=avma;
  GEN nf=alg_get_center(A);//Field
  GEN rams=algramifiedplacesf(A);
  GEN algdisc=gen_1;
  for(long i=1;i<lg(rams);i++) algdisc=mulii(algdisc, idealnorm(nf, gel(rams, i)));//Norm to Q of the ramification
  return gerepileupto(top, algdisc);
}

//Returns the vector of finite ramified places of the algebra A.
GEN algramifiedplacesf(GEN A){
  pari_sp top=avma;
  GEN hass=alg_get_hasse_f(A);//Shallow
  long nhass=lg(gel(hass, 2));
  GEN rp=vectrunc_init(nhass);
  for(long i=1;i<nhass;i++){
    if(gel(hass, 2)[i]==0) continue;//Unramified
	vectrunc_append(rp, gel(gel(hass, 1), i));//Ramified
  }
  return gerepilecopy(top, rp);
}

//Returns small norm 1 elements (Q_{z1,z2}(x)<=C) of the order in A
GEN algsmallnorm1elts(GEN A, GEN p, GEN C, GEN z1, GEN z2, int type, long prec){
  pari_sp top=avma;
  GEN Q=qalg_fdominitialize(A, prec);
  GEN nf=alg_get_center(A);
  long nfdeg=nf_get_degree(nf), fourn=4*nfdeg;
  GEN nformpart=qalg_normform(Q);
  GEN nfdecomp, nform;//nfdecomp used if type=1, nform if type=2
  if(type==1 || (type==0 && nfdeg>=2)){//qfminim
    type=1;//Setting in case type=0
	nfdecomp=mat_nfcholesky(nf, nformpart);
  }
  else{//type=2 or type=0 and nfdeg=1, i.e. nf=Q. condition
    type=2;//Setting in case type=0
	nform=gcopy(nformpart);
  }
  for(long i=1;i<=fourn;i++){
	for(long j=1;j<=fourn;j++){
	  gcoeff(nformpart, i, j)=nftrace(nf, gcoeff(nformpart, i, j));//Taking the trace to Q
	}
  }//Tr_{nf/Q}(nrd(elt));
  if(type==1) return gerepileupto(top, qalg_smallnorm1elts_qfminim(Q, p, C, z1, z2, 0, nfdecomp, nformpart, prec));
  return gerepileupto(top, qalg_smallnorm1elts_condition(Q, p, C, z1, z2, 0, nform, nformpart, prec)); 
}



//FUNDAMENTAL DOMAIN COMPUTATION


//Generate the fundamental domain for a quaternion algebra initialized with alginit. We can pass C, R, type, and they will be auto-set if passed as 0.
static GEN qalg_fdom(GEN Q, GEN p, int dispprogress, int dumppartial, GEN partialset, GEN C, GEN R, GEN passes, int type, GEN tol, long prec){
 pari_sp top=avma, mid;
  GEN mats=psltopsu_transmats(p);
  GEN A=qalg_get_alg(Q);
  GEN nf=alg_get_center(A);
  long nfdeg=nf_get_degree(nf);
  GEN area=qalg_fdomarea(Q, 3, prec);//Smallest precision possible.
  GEN areabound=gdivgs(gmulgs(area, 3), 2);//Times 1.5.

  if(gequal0(passes)){//Setting passes
    if(nfdeg==1) passes=gen_2;
	else passes=stoi(12);
  }
  if(gequal0(C)) C=algfdom_bestC(A, prec);//Setting C
  GEN N=gceil(gdiv(gsqr(area), gmul(gmul(Pi2n(3, prec), gsubgs(C, nfdeg)), passes)));//Area^2/(8*Pi*(C-n)*#Passes)
  if(gcmpgs(N, 1)<=0) N=gen_2;//Make sure N>=2
  GEN gamma=dbltor(2.1);//2.1
  if(gequal0(R)) R=hdiscradius(gpow(area, gamma, prec), prec);//Setting R
  GEN epsilon=gdivgs(gen_1, 6);
  
  if(dispprogress) pari_printf("Initial constants:\n   C=%P.8f\n   N=%Ps\n   R=%P.8f\nTarget Area: %P.8f\n\n", C, N, R, area);
  
  GEN id=gel(alg_get_basis(A), 1);//The identity  
  long iN;
  GEN points, w;
  
  long fourn=4*nfdeg;//The lg of a normal entry
  GEN nformpart=qalg_normform(Q);
  GEN nfdecomp, nform;//nfdecomp used if type=1, nform if type=2
  if(lg(qalg_get_rams(Q))==1) type=1;//We must use qfminim, since our implementation of condition won't work in a non-division algebra.
  long maxelts=1;
  if(type==1 || (type==0 && nfdeg>=2)){//qfminim
    type=1;//Setting in case type=0
	nfdecomp=mat_nfcholesky(nf, nformpart);
  }
  else{//type=2 or type=0 and nfdeg=1, i.e. nf=Q. condition
    type=2;//Setting in case type=0
	nform=gcopy(nformpart);
  }
  for(long i=1;i<=fourn;i++){
	for(long j=1;j<=fourn;j++){
	  gcoeff(nformpart, i, j)=nftrace(nf, gcoeff(nformpart, i, j));//Taking the trace to Q
	}
  }//Tr_{nf/Q}(nrd(elt));
  
  GEN U=cgetg(2, t_VEC);
  gel(U, 1)=cgetg(1, t_VEC);//Setting U=[[]], so that the first time normalizedbasis is called, it works
  if(gequal0(partialset)){
	if(type==1) partialset=qalg_smallnorm1elts_qfminim(Q, p, C, gen_0, gen_0, 0, nfdecomp, nformpart, prec);
	else partialset=qalg_smallnorm1elts_condition(Q, p, C, gen_0, gen_0, 0, nform, nformpart, prec);
  }
  U=normalizedbasis(partialset, U, mats, id, &Q, &qalg_fdomm2rembed, &qalg_fdommul, &qalg_fdominv, &qalg_istriv, tol, prec);

  FILE *f;
  if(dumppartial) f=fopen("algfdom_partialdata_log.txt", "w");

  mid=avma;
  long pass=0, ooend=0, nskip, nsidesp1=lg(gel(U, 1));//How many sides+1; pass tracks which pass we are on
  GEN oosides=cgetg(1, t_VEC), ang1, ang2;//Initialize, so that lg(oosides)=1 to start.
  int anyskipped=0;

  for(;;){
	pass++;
	if(dispprogress){pari_printf("Pass %d with %Ps random points in the ball of radius %P.8f\n", pass, N, R);}
	if(nsidesp1>1){//We have a partial domain.
	  oosides=normalizedboundary_oosides(U);
	  ooend=lg(oosides)-1;
	  if(dispprogress) pari_printf("%d infinite sides\n", ooend);
	}
	iN=itos(gfloor(N))+1;
	nskip=0;//How many points are skipped due to poor precision
	points=cgetg(iN, t_VEC);
	for(long i=1;i<iN;i++){//Random points in ball of radius R
	  if(0<ooend){//Going near infinite side
	    long iside=((i-1)%ooend)+1;
	    ang2=gel(gel(U, 4), oosides[iside]);
	    if(oosides[iside]==1) ang1=gel(gel(U, 4), lg(gel(U, 1))-1);//Last side, which is the previous side
	    else ang1=gel(gel(U, 4), oosides[iside]-1);
	    w=randompoint_udarc(R, ang1, ang2, prec);
	  }
	  else w=randompoint_ud(R, prec);//Random point
	  pari_CATCH(CATCH_ALL){
	    GEN err=pari_err_last();//Get the error
		long errcode=err_get_num(err);
		pari_printf("");//Seems to be required, else I get a segmentation fault? weird...
		if(errcode==e_TYPE){//If R is large
	      if(!anyskipped) pari_printf("Point skipped, consider stopping the job increasing the precision (or decreasing R) to avoid this\n");
		  anyskipped=1;
		  nskip++;
		  gel(points, i)=cgetg(1, t_VEC);
		}
		else pari_err(errcode);
	  }
	  pari_TRY{
		GEN smallelts;
		if(type==1) smallelts=qalg_smallnorm1elts_qfminim(Q, p, C, gen_0, w, maxelts, nfdecomp, nformpart, prec);
	    else smallelts=qalg_smallnorm1elts_condition(Q, p, C, gen_0, w, maxelts, nform, nformpart, prec);
		if(smallelts) gel(points, i)=smallelts;
		else gel(points, i)=cgetg(1, t_VEC);//There was an issue (possibly precision induced)
	  }
	  pari_ENDCATCH
	}
	points=shallowconcat1(points);
	if(dispprogress){
	  if(nskip) pari_printf("%d points skipped due to lack of precision\n", nskip);
	  pari_printf("%d elements found\n", lg(points)-1);
	}
	U=normalizedbasis(points, U, mats, id, &Q, &qalg_fdomm2rembed, &qalg_fdommul, &qalg_fdominv, &qalg_istriv, tol, prec);
	if(dispprogress) pari_printf("Current normalized basis has %d sides\n\n", lg(gel(U, 1))-1);
	if(gcmp(gel(U, 6), areabound)==-1){
	  if(dumppartial) fclose(f);
	  return gerepileupto(top, U);
	}
	if(pass>1 && (ooend==0 || nsidesp1==lg(gel(U, 1)))) R=gadd(R, epsilon);//Updating R_n
	nsidesp1=lg(gel(U, 1));//How many sides+1
	if(gc_needed(top, 2)) gerepileall(mid, 3, &U, &N, &R);
	if(dumppartial) pari_fprintf(f, "%Ps\n", gel(U, 1));
  }
}


//HELPER METHODS


//If A is an algebra over nf, let decomp be the Cholesky decomposition of the norm form. This returns the norm of x given the decomposition (this is ~10x faster than algnorm).
GEN algnorm_chol(GEN nf, GEN decomp, GEN x){
  pari_sp top=avma;
  GEN part=gen_0, U;
  long n=lg(x);
  for(long i=n-1;i>0;i--){
	if(gequal0(gcoeff(decomp, i, i))) continue;//This will happen for all but 4 indices
	U=gel(x, i);
	for(long j=i+1;j<n;j++) U=nfadd(nf, U, nfmul(nf, gcoeff(decomp, i, j), gel(x, j)));
	U=nfmul(nf, gcoeff(decomp, i, i), nfsqr(nf, U));
	part=nfadd(nf, part, U);
  }
  return gerepileupto(top, part);
}

//If the algebra A has a unique split infinite place, this returns the index of that place. Otherwise, returns 0.
static long algsplitoo(GEN A){
  GEN infram=alg_get_hasse_i(A);//shallow
  long split=0;
  for(long i=1;i<lg(infram);i++){//Finding the split place
	if(infram[i]==0){//Split place
	  if(split>0) return 0;
	  split=i;
	}
  }
  return split;//No garbage!!
}

//Returns the qf absrednorm as a matrix. If order has basis v1, ..., vn, then this is absrednorm_{z1, z2}(e1*v1+...+en*vn), defined in the Enumeration section of my paper. This is the combination of the Definition on page 478 of Voight, as well as Page's definition for Kleinian groups. It satisfies Q_{z1, z2}(g)=cosh(d(gz1, z2))+n-1 if g has norm 1 in the order. Thus, this allows us to see if z_1 and z_2 are "close". on the quotient.
GEN qalg_absrednormqf(GEN Q, GEN mats, GEN z1, GEN z2, GEN normformpart, long prec){
  pari_sp top=avma;

  //First we find the part coming from |f_g(p)|
  GEN A=qalg_get_alg(Q);
  long n=lg(gel(alg_get_basis(A), 1));//The lg of a normal entry
  GEN basisimage=cgetg(n, t_VEC);//The image of the basis elements in M_2(R)
  GEN belt=zerocol(n-1);
  gel(belt, 1)=gen_1;
  gel(basisimage, 1)=qalg_fdomm2rembed(&Q, belt, prec);

  for(long i=2;i<n;i++){
	gel(belt, i-1)=gen_0;
	gel(belt, i)=gen_1;//Updating belt to have a 1 only in the ith place
	gel(basisimage, i)=qalg_fdomm2rembed(&Q, belt, prec);
  }

  GEN tvars=cgetg(n, t_VECSMALL);
  GEN xvars=cgetg(n, t_VEC);
  for(long i=1;i<n;i++){
	tvars[i]=fetch_var();//The temporary variable numbers
	gel(xvars, i)=pol_x(tvars[i]);//The variables
  }

  GEN elt=zeromatcopy(2, 2);//Will represent the general element
  for(long i=1;i<n;i++) elt=gadd(elt, gmul(gel(basisimage, i), gel(xvars, i)));//The general element
  if(!gequal0(z2)){//Shift by h2^(-1) on the left, where h2=phi^(-1)*W*phi, and W=[a, b; conj(b), a] with a=1/sqrt(1-|z|^2) and b=za
	GEN a=gdivsg(1, gsqrt(gsubsg(1, gsqr(gabs(z2, prec))), prec));//1/sqrt(1-|z2|^2)
	GEN mb=gneg(gmul(z2, a));//-z2*a=-b
	GEN Winv=mkmat22(a, mb, conj_i(mb), a);//W^(-1)
	GEN h2inv=gmul(gel(mats, 2), gmul(Winv, gel(mats, 1)));//h2^(-1); We have to translate back to PSL, hence why we do mats[2]*Winv*mats[1] (mats[1],[2] are inverses to each other.
	elt=gmul(h2inv, elt);//Shifting it by the appropriate amount.
  }
  if(!gequal0(z1)){//Shift by h1 on the right, where h1=phi^(-1)*W*phi, and W=[a,b;conj(b), a] with a=1/sqrt(1-|z|^2) and b=za
	GEN a=gdivsg(1, gsqrt(gsubsg(1, gsqr(gabs(z1, prec))), prec));//1/sqrt(1-|z1|^2)
	GEN b=gmul(z1, a);//z1*a=b
	GEN W=mkmat22(a, b, conj_i(b), a);//W
	GEN h1=gmul(gel(mats, 2), gmul(W, gel(mats, 1)));//h1; We have to translate back to PSL, hence why we do mats[2]*Winv*mats[1] (mats[1],[2] are inverses to each other.
	elt=gmul(elt, h1);//Shifting it by the appropriate amount.
  }
  GEN p=gel(mats, 3);//p
  GEN fg=gmul(gcoeff(elt, 2, 1), p);//elt[2, 1]*p
  fg=gadd(fg, gsub(gcoeff(elt, 2, 2), gcoeff(elt, 1, 1)));//elt[2, 1]*p+elt[2, 2]-elt[2, 1]
  fg=gsub(gmul(fg, p), gcoeff(elt, 1, 2));//elt[2, 1]*p^2+(elt[2, 2]-elt[2, 1])*p-elt[1, 2], i.e. f_g(p) (page 477 of Voight)
  GEN invrad1=gadd(gsqr(real_i(fg)), gsqr(imag_i(fg)));//real(fg)^2+imag(fg)^2

  GEN K=alg_get_center(A);
  GEN invrad2;
  if(gequal0(normformpart)){//If not zero we pass it into the method.
    invrad2=qalg_normform(Q);
    for(long i=1;i<n;i++){
	  for(long j=1;j<n;j++){
	    gcoeff(invrad2, i, j)=nftrace(K, gcoeff(invrad2, i, j));//Taking the trace to Q
	  }
    }//invrad2=Tr_{K/Q}(nrd(elt));
  }
  else invrad2=normformpart;
  //At this point, invrad1 is a polynomial in the temporary variables, and invrad2 is a matrix in the correct form.*/

  GEN pscale=gdivsg(1, gmulsg(2, gsqr(imag_i(p))));//1/2y^2 for p=x+iy
  GEN qf=cgetg(n, t_MAT);//Creating the return matrix
  for(long i=1;i<n;i++) gel(qf, i)=cgetg(n, t_COL);
  GEN part1=invrad1;//This stores the part we are working on.
  long var=qalg_get_varnos(Q)[1];//The variable number for K
  GEN Kx=gel(qalg_get_roots(Q), 1);//The approximation of K
  for(long i=1;i<n;i++){//Working with the ith variable
	gcoeff(qf, i, i)=gadd(gmul(gsubst(lift(polcoef_i(part1, 2, tvars[i])), var, Kx), pscale), gcoeff(invrad2, i, i));//iith coefficient
	GEN linpart=polcoef_i(part1, 1, tvars[i]);//The part that is linear in the ith coefficient.
	for(long j=i+1;j<n;j++){
	  gcoeff(qf, i, j)=gadd(gmul(gdivgs(gsubst(lift(polcoef_i(linpart, 1, tvars[j])), var, Kx), 2), pscale), gcoeff(invrad2, i, j));//the ijth coefficient
	  gcoeff(qf, j, i)=gcoeff(qf, i, j);//Okay as we will copy it later
	}
	part1=polcoef_i(part1, 0, tvars[i]);//The part that has no v_i
  }
  long delv;
  do{delv=delete_var();} while(delv && delv<=tvars[1]);//Delete the temporary variables
  return gerepilecopy(top, qf);
}

//Given the (assumed) basis representation of an element, this returns the conjugate in the basis representation.
static GEN qalg_basis_conj(GEN Q, GEN x){
  pari_sp top=avma;
  GEN A=qalg_get_alg(Q);
  GEN varnos=qalg_get_varnos(Q);
  GEN Lx=pol_x(varnos[2]);//The variable for L
  GEN algx=liftall(algbasistoalg(A, x));
  gel(algx, 1)=gsubst(gel(algx, 1), varnos[2], gneg(Lx));//Conjugating
  gel(algx, 2)=gneg(gel(algx, 2));//The conjugate of l1+jl2 is sigma(l1)-jl2
  return gerepileupto(top, algalgtobasis(A, algx));//Converting back
}

//Returns the area of the fundamental domain of Q, computed to computeprec. Assumes the order is MAXIMAL FOR NOW (see Voight Theorem 39.1.18 to update for general; very easy). We also submit the old precision since the loss of precision causes errors later.
GEN qalg_fdomarea(GEN Q, long computeprec, long prec){
  pari_sp top=avma;
  long bits=bit_accuracy(computeprec);
  GEN A=qalg_get_alg(Q);
  GEN F=alg_get_center(A), pol=nf_get_pol(F);
  GEN zetaval=lfun(pol, gen_2, bits);//Zeta_F(2)
  GEN rams=qalg_get_rams(Q);
  GEN norm=gen_1;
  for(long i=1;i<lg(rams);i++){
	if(typ(gel(rams, i))==t_INT) continue;//We don't want to count the infinite places
	norm=mulii(norm, subis(idealnorm(F, gel(rams, i)), 1));//Product of N(p)-1 over finite p ramifying in A
  }
  GEN ar=gmul(gpow(nfdisc(pol), gdivsg(3, gen_2), computeprec), norm);//d_F^(3/2)*phi(D)
  long n=nf_get_degree(F), twon=2*n;
  ar=gmul(ar, zetaval);//zeta_F(2)*d_F^(3/2)*phi(D)
  ar=gmul(ar, gpowgs(gen_2, 3-twon));
  ar=gmul(ar, gpowgs(mppi(computeprec), 1-twon));//2^(3-2n)*Pi^(1-2n)*zeta_F(2)*d_F^(3/2)*phi(D)
  return gerepileupto(top, gtofp(ar, prec));
}

//Initializes the quaternion algebra Q split at one real place using the algebras framework. Assume that A is input as a quaternion algebra with pre-computed maximal order. This is not suitable for gerepile.
GEN qalg_fdominitialize(GEN A, long prec){
  GEN K=alg_get_center(A);//The centre, i.e K where A=(a,b/K)
  GEN L=alg_get_splittingfield(A);//L=K(sqrt(a)).
  GEN ramdat=algramifiedplacesf(A);//Finite ramified places
  GEN varnos=mkvecsmall2(nf_get_varn(K), rnf_get_varn(L));//Stores the variable numbers for K and L, in that order
  long split=algsplitoo(A);//The split real place
  if(split==0) pari_err_TYPE("Quaternion algebra has 0 or >=2 split real infinite places, not valid for fundamental domains.", A);  
  GEN Kroot=gel(nf_get_roots(K), split);//The split root
  GEN aval=poleval(gneg(polcoef_i(rnf_get_pol(L), varnos[2], 0)), Kroot);//Find the defining eqn of L (x^2+u for u in K), find u, then take -u=a
  if(gsigne(aval)!=1) pari_err_TYPE("We require a>0 at the split real place. Please swap a, b.", A);
  GEN roots=mkvec2(Kroot, gsqrt(aval, prec));//The approximate values for the variables defining K and L.
  return mkvec4(A, ramdat, varnos, roots);
}

//Returns the norm form as a matrix, i.e. M such that nrd(e_1*v_1+...+e_nv_n)=(e_1,...,e_n)*M*(e_1,...,e_n)~, where v_1, v_2, ..., v_n is the given basis of A. The iith coefficient is nrd(v_i), and the ijth coefficient (i!=j) is .5*trd(v_iv_j)
GEN qalg_normform(GEN Q){
  pari_sp top=avma;
  GEN A=qalg_get_alg(Q);
  long n=lg(gel(alg_get_basis(A), 1));//The lg of a normal entry
  GEN basis=matid(n-1);
  return gerepileupto(top, qalg_normform_givenbasis(Q, basis));
}

//Basis is a matrix whose columns span a lattice, say v_1, ..., v_k. This returns M such that nrd(e_1*v_1+...+e_k*v_k)=(e1,...,e_k)*M*(e_1,...,e_k)~. The iith coefficient is nrd(v_i) and the ijth coefficient if i!=j is .5*trd(v_iv_j).
static GEN qalg_normform_givenbasis(GEN Q, GEN basis){
  pari_sp top=avma;
  GEN A=qalg_get_alg(Q);
  long n;
  GEN basisconj=cgetg_copy(basis, &n);
  for(long i=1;i<n;i++) gel(basisconj, i)=qalg_basis_conj(Q, gel(basis, i));//The conjugate of the basis element.
  GEN M=cgetg_copy(basis, &n);//Initialize the matrix
  for(long i=1;i<n;i++) gel(M, i)=cgetg(n, t_COL);
  for(long i=1;i<n;i++) gcoeff(M, i, i)=lift0(algnorm(A, gel(basis, i), 0), -1);
  for(long i=1;i<n;i++){
	for(long j=i+1;j<n;j++){
	  GEN prod=algmul(A, gel(basis, i), gel(basisconj, j));
	  GEN tr=gdivgs(algtrace(A, prod, 0), 2);
	  gcoeff(M, i, j)=tr;
	  gcoeff(M, j, i)=tr;//OK since we copy at the end.
	}
  }
  return gerepilecopy(top, M);
}

//Computes all norm 1 elements for which Q_{z_1,z_2}(x)<=C. If z1 and z2 are close on the Shimura curve, then this should return a point. maxelts is the maximum number of return elements (or 0 for all norm 1 elts)
GEN qalg_smallnorm1elts_qfminim(GEN Q, GEN p, GEN C, GEN z1, GEN z2, long maxelts, GEN nfdecomp, GEN nformpart, long prec){
  pari_sp top=avma, mid;
  GEN A=qalg_get_alg(Q);
  GEN nf=alg_get_center(A);
  GEN mats=psltopsu_transmats(p);
  GEN absrednorm=qalg_absrednormqf(Q, mats, z1, z2, nformpart, prec);
  GEN vposs=gel(qfminim0(absrednorm, C, NULL, 2, prec), 3);
  long nvposs=lg(vposs), mret;
  if(maxelts) mret=maxelts+1;
  else mret=nvposs;
  GEN ret=vectrunc_init(mret), norm;
  for(long i=1;i<nvposs;i++){
	mid=avma;
	norm=algnorm_chol(nf, nfdecomp, gel(vposs, i));
	if(gequal(norm, gen_1)){
	  avma=mid;
	  vectrunc_append(ret, gel(vposs, i));//Don't append a copy, will copy at the end.
	  if(lg(ret)>=mret) break;//Done
	  continue;
	}
	avma=mid;
  }
  return gerepilecopy(top, ret);
}

//Computes G(C) ala Voight, i.e. elements of O_{N=1} with a large radius near v.
GEN qalg_smallnorm1elts_condition(GEN Q, GEN p, GEN C, GEN z1, GEN z2, long maxelts, GEN nform, GEN nformpart, long prec){
  pari_sp top=avma;
  GEN mats=psltopsu_transmats(p);
  GEN absrednorm=qalg_absrednormqf(Q, mats, z1, z2, nformpart, prec);
  GEN A=qalg_get_alg(Q);
  return gerepileupto(top, smallvectors_nfcondition(absrednorm, C, maxelts, mkvec3(alg_get_center(A), nform, gen_1), prec));
}


//BASIC OPERATIONS FOR NORMALIZED BASIS ET AL


//Must pass *data as a quaternion algebra. This just formats things correctly for the fundamental domain.
GEN qalg_fdominv(GEN *data, GEN x){
  return alginv(qalg_get_alg(*data), x);
}

//Must pass *data as a quaternion algebra. This embeds the element x into M_2(R), via l1+jl2->[l1, sigma(l2)b;l2, sigma(l1)].
GEN qalg_fdomm2rembed(GEN *data, GEN x, long prec){
  pari_sp top=avma;
  GEN A=qalg_get_alg(*data);
  GEN rts=qalg_get_roots(*data);
  GEN varnos=qalg_get_varnos(*data); 
  GEN b=gsubst(lift(alg_get_b(A)), varnos[1], gel(rts, 1));//Inputting the real value of K into b
  if(lg(x)!=3) x=algbasistoalg(A, x);//Converting it to algebraic representation
  x=liftall(x);//Lifiting x
  GEN l1K=gsubst(gel(x, 1), varnos[1], gel(rts, 1));//Inputting the real value of K into l1
  GEN l2K=gsubst(gel(x, 2), varnos[1], gel(rts, 1));//Inputting the real value of K into l2
  GEN mLrt=gneg(gel(rts, 2));//The negative root of L
  GEN M=cgetg(3, t_MAT);
  gel(M, 1)=cgetg(3, t_COL), gel(M, 2)=cgetg(3, t_COL);//Initializing the matrix
  gcoeff(M, 1, 1)=gsubst(l1K, varnos[2], gel(rts, 2));//l1
  gcoeff(M, 1, 2)=gmul(gsubst(l2K, varnos[2], mLrt), b);//sigma(l2)*b, i.e. substitute in the negative root of L into l2
  gcoeff(M, 2, 1)=gsubst(l2K, varnos[2], gel(rts, 2));//l2
  gcoeff(M, 2, 2)=gsubst(l1K, varnos[2], mLrt);//sigma(l1)
  return gerepilecopy(top, M);
}

//Must pass *data as a quaternion algebra. This just formats things correctly for the fundamental domain.
GEN qalg_fdommul(GEN *data, GEN x, GEN y){
  return algmul(qalg_get_alg(*data), x, y);
}

//Must pass *data as a quaternion algebra. Returns the trace of x.
GEN qalg_fdomtrace(GEN *data, GEN x){
  return algtrace(qalg_get_alg(*data), x, 0);
}

//Returns 1 if x==+/-1. x must be in the basis representation (note that the first element of the basis is always 1).
int qalg_istriv(GEN *data, GEN x){
  if(!gequal(gel(x, 1), gen_1) && !gequal(gel(x, 1), gen_m1)) return 0;
  for(long i=2;i<lg(x);i++) if(!gequal0(gel(x, i))) return 0;
  return 1;
}


//SHALLOW RETRIEVAL METHODS


//Shallow method to return the algebra
GEN qalg_get_alg(GEN Q){return gel(Q, 1);}

//Shallow method to get the ramification
GEN qalg_get_rams(GEN Q){return gel(Q, 2);}

//Shallow method to return the variable numbers of K and L
GEN qalg_get_varnos(GEN Q){return gel(Q, 3);}

//Shallow method to return the roots of K and L.
GEN qalg_get_roots(GEN Q){return gel(Q, 4);}

