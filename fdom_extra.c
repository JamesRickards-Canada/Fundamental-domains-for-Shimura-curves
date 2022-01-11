//These methods are useful with the fdom.c package, but not essential.

//INCLUSIONS

#ifndef PARILIB
#define PARILIB
#include <pari/pari.h>
#endif

#ifndef METHDECL
#define METHDECL
#include "fdomdecl.h"
#endif

//STATIC DECLARATIONS

//2: QUATERNION ALGEBRA METHODS
static GEN algd(GEN A, GEN a);
static GEN algfromnormdisc(GEN F, GEN D, GEN infram);
static int subset_lexincrem(GEN S, long n);

//3: OPTIMIZING THE VALUE OF C FOR ENUMERATION
static GEN enum_bestC_givenAB(long n, GEN A, GEN B, long prec);
static void enum_bestC_plot(GEN reg, GEN Cmin, GEN Cmax, long n, char *fdata, int isArange);
static long enum_nontrivial(GEN L);
static GEN enum_successrate_givendata(GEN Q, GEN p, GEN C, long Ntests, GEN R, GEN area, GEN normdecomp, GEN normformpart, long prec);
static void enum_successrate_plot(GEN A, GEN reg, GEN Cmin, GEN Cmax, char *fdata, int WSL);
static void enum_time_plot(GEN A, GEN reg, GEN Cmin, GEN Cmax, char *fdata);
static long enum_timeforNelts_givendata(GEN Q, GEN p, GEN C, long nelts, GEN R, int type, GEN nform, GEN nfdecomp, GEN nformpart, long prec);
static void enum_timeforNelts_plot(GEN A, GEN bestC, double bestC_val, char *fname);

//3: NUMBER OF ELEMENTS REQUIRED TO GENERATE ALGEBRA
static GEN qalg_fdom_nelts(GEN Q, GEN p, GEN CNRdata, int type, GEN tol, long prec);


//SECTION 1: GEOMETRY METHODS


//Writes the set of arcs in the given colour to filename, so that it is ready to be executed by python.
void
python_printarcs(GEN arcs, char *filename, int view, char *extrainput, long prec)
{
  pari_sp top=avma;
  if(!pari_is_dir("fdoms")){//Checking the directory
    int s=system("mkdir -p fdoms");
    if(s==-1) pari_err(e_MISC, "ERROR CREATING DIRECTORY fdoms");
  }
  char *fullfile=pari_sprintf("fdoms/%s.dat", filename);
  FILE *f=fopen(fullfile, "w");
  pari_free(fullfile);//Now we have created the output file f.
  GEN arc, fact=gdiv(stoi(180), mppi(prec));//fact=180/Pi
  for(long i=1;i<lg(arcs);i++){
    arc=gel(arcs, i);
    if(gequal0(arc)) continue;//Not a circle
    if(gequal0(gel(arc, 8))){//Arc
      pari_fprintf(f, "0 %lf %lf %lf %lf %lf %d\n", rtodbl(gtofp(real_i(gel(arc, 1)), prec)), rtodbl(gtofp(imag_i(gel(arc, 1)), prec)), rtodbl(gtofp(gel(arc, 2), prec)), rtodbl(gtofp(gmul(gel(arc, 5), fact), prec)), rtodbl(gtofp(gmul(gel(arc, 6), fact), prec)), itos(gel(arc, 7)));
    }
    else{//Segment
      pari_fprintf(f, "1 %lf %lf %lf %lf\n", rtodbl(gtofp(real_i(gel(arc, 3)), prec)), rtodbl(gtofp(imag_i(gel(arc, 3)), prec)), rtodbl(gtofp(real_i(gel(arc, 4)), prec)), rtodbl(gtofp(imag_i(gel(arc, 4)), prec)));
    }
  }
  fclose(f);
  if(view==1){
    char *line;
    if(extrainput==NULL) line=pari_sprintf("%s", filename);
    else line=pari_sprintf("%s %s", extrainput, filename);
    python_plotviewer(line);
    pari_free(line);
  }
  set_avma(top);
}

//Launches the plotviewer with the given inputs.
void
python_plotviewer(char *input)
{
  char *command;
  command=pari_sprintf("cmd.exe /C start py fdviewer.py %s", input);
  int s=system(command);
  if(s==-1) pari_err(e_MISC, "ERROR EXECUTING COMMAND");
  pari_free(command);
}

//Writes the fundamental domain corresponding to U
void
python_printfdom(GEN U, char *filename, long prec)
{
  pari_sp top=avma;
  if(!pari_is_dir("fdoms")){//Checking the directory
    int s=system("mkdir -p fdoms");
    if(s==-1) pari_err(e_MISC, "ERROR CREATING DIRECTORY fdoms");
  }
  char *fullfile=pari_sprintf("fdoms/%s.dat", filename);
  FILE *f=fopen(fullfile, "w");
  pari_free(fullfile);//Now we have created the output file f.
  GEN pair=gel(U, 7);
  pari_fprintf(f, "%d", pair[1]);
  for(long i=2;i<lg(pair);i++) pari_fprintf(f, " %d", pair[i]);//Print side pairing.
  //GEN vangs=gel(U, 4);
  //pari_fprintf(f, "\n%lf", rtodbl(gtofp(gel(vangs, 1), prec)));//First angle
  //for(long i=2;i<lg(pair);i++) pari_fprintf(f, " %lf", rtodbl(gtofp(gel(vangs, i), prec)));//Print angles.
  pari_fprintf(f, "\n");
  GEN arcs=gel(U, 2);
  GEN arc, fact=gdiv(stoi(180), mppi(prec)), v1, v2;//fact=180/Pi
  for(long i=1;i<lg(arcs);i++){
    arc=gel(arcs, i);
    if(gequal0(arc)) continue;//Not a circle
    if(i==1) v1=gmael(U, 3, lg(arcs)-1);
    else v1=gmael(U, 3, i-1);
    v2=gmael(U, 3, i);//The two vertices
    pari_fprintf(f, "%lf %lf %lf %lf %lf\n", rtodbl(gtofp(real_i(gel(arc, 1)), prec)), rtodbl(gtofp(imag_i(gel(arc, 1)), prec)), rtodbl(gtofp(gel(arc, 2), prec)), rtodbl(gmul(garg(gsub(v1, gmael3(U, 2, i, 1)), prec), fact)), rtodbl(gmul(garg(gsub(v2, gmael3(U, 2, i, 1)), prec), fact)));
  }
  fclose(f);
  set_avma(top);
}



//SECTION 2: QUATERNIONIC METHODS


//Given x=[e, f, g, h], returns x in the algebra representation.
GEN
alg1ijktoalg(GEN A, GEN x)
{
  pari_sp top=avma;
  GEN L=alg_get_splittingfield(A);//L=F(i),
  long Lvar=rnf_get_varn(L);
  GEN Lx=pol_x(Lvar);//x
  GEN e1=gadd(gel(x, 1), gmul(gel(x, 2), Lx));//e+fi
  GEN e2=gsub(gel(x, 3), gmul(gel(x, 4), Lx));//g-hi
  return gerepilecopy(top, mkcol2(e1, e2));
}

//Given x=[e, f, g, h], returns x in the basis representation.
GEN
alg1ijktobasis(GEN A, GEN x)
{
  pari_sp top=avma;
  GEN xalg=alg1ijktoalg(A, x);
  return gerepileupto(top, algalgtobasis(A, xalg));
}

//Given an element in the algebra representation of A, returns [e, f, g, h], where x=e+fi+gj+hk. e, f, g, h live in the centre of A.
GEN
algalgto1ijk(GEN A, GEN x)
{
  pari_sp top=avma;
  x=liftall(x);//Make sure there are no mods.
  GEN L=alg_get_splittingfield(A);//L=F(i),
  long Lvar=rnf_get_varn(L);
  GEN f=polcoef_i(gel(x, 1), 1, Lvar);//Extract f
  GEN e=polcoef_i(gel(x, 1), 0, Lvar);//Extract e
  GEN g=polcoef_i(gel(x, 2), 0, Lvar);//Extract g
  GEN mh=polcoef_i(gel(x, 2), 1, Lvar);//Extract -h
  return gerepilecopy(top, mkvec4(e, f, g, gneg(mh)));
}

//Given an element in the basis representation of A, returns [e, f, g, h], where x=e+fi+gj+hk. e, f, g, h live in the centre of A.
GEN
algbasisto1ijk(GEN A, GEN x)
{
  pari_sp top=avma;
  GEN xalg=liftall(algbasistoalg(A, x));
  return gerepileupto(top, algalgto1ijk(A, xalg));
}

//Given a1, ..., a4 in A, this returns d(a1,...,a4), i.e. det(trd(a_i*a_j))
static GEN
algd(GEN A, GEN a)
{
  pari_sp top=avma;
  GEN M=cgetg(5, t_MAT);
  for(long i=1;i<5;i++) gel(M, i)=cgetg(5, t_COL);
  for(long i=1;i<5;i++){//Setting the coefficients of M
    gcoeff(M, i, i)=algtrace(A, algsqr(A, gel(a, i)), 0);
    for(long j=i+1;j<5;j++) gcoeff(M, i, j)=gcoeff(M, j, i)=algtrace(A, algmul(A, gel(a, i), gel(a, j)), 0);
  }
  return gerepileupto(top, nfM_det(alg_get_center(A), M));
}

//Returns the Eichler order given by O intersect x O x^(-1); the output is [Z-basis of order, level].
GEN
algeichler_conj(GEN A, GEN x)
{
  pari_sp top=avma;
  long n=lg(x)-1;
  GEN O=matid(n);//The original order.
  GEN Oconj=algorderconj(A, x, O);
  GEN combine=shallowconcat(O, Oconj);//Concatenate the orders
  GEN ker=kerint(Q_primpart(combine));//The kernel
  GEN H=hnf(ker);//hnf
  GEN newO=hnf(matslice(H, 1, n, 1, n));//new order
  return gerepilecopy(top, mkvec2(newO, algorderlevel(A, newO, 0)));
}

//Returns a quaternion algebra over F (of degree n) with |N_{F/Q}(discriminant)|=D and infinite ramification prescribed by infram (a length n vector of 0's/1's), if it exists. If it does not, this returns 0. This is not gerepile suitable, it leaves a dirty stack. Actually, this MISSES the q-algebras with multiple ramifying primes that have the SAME norm. But that's okay.
static GEN
algfromnormdisc(GEN F, GEN D, GEN infram)
{
  pari_sp top=avma;
  if(typ(D)!=t_INT || signe(D)!=1) pari_err_TYPE("D should be a positive integer", D);//The absolute norm to Q should be a positive integer
  GEN pfac=Z_factor(D);
  long nfacsp1=lg(gel(pfac, 1));//# prime factors+1
  long nramplaces=nfacsp1-1;
  for(long i=1;i<lg(infram);i++) if(!gequal0(gel(infram, i))) nramplaces++;//Adding the number of oo ramified places
  if(nramplaces%2==1){set_avma(top);return gen_0;}//Odd number of ramification places, BAD
  GEN pfacideals=zerovec(nfacsp1-1), hass=cgetg(nfacsp1, t_VEC), possideals;
  long expon;
  for(long i=1;i<nfacsp1;i++){
    expon=itos(gcoeff(pfac, i, 2));//Coefficient of p we desire.
    possideals=idealprimedec(F, gcoeff(pfac, i, 1));
    for(long j=1;j<lg(possideals);j++){
      if(pr_get_f(gel(possideals, j))==expon){
        gel(pfacideals, i)=gel(possideals, j);//We win!
        break;
      }
    }
    if(gequal0(gel(pfacideals, i))){set_avma(top);return gen_0;}//Nope, return 0
    gel(hass, i)=gen_1;//The vector of 1's
  }
  return alginit(F, mkvec3(gen_2, mkvec2(pfacideals, hass), infram), -1, 1);
}

//Returns the order formed by conjugating O by x. The columns are written in terms of the stored basis of A. If O=NULL, assumes we are conjugating the stored maximal order
GEN
algorderconj(GEN A, GEN x, GEN O)
{
  pari_sp top=avma;
  long n=lg(x)-1;
  if(!O) O=matid(n);//Setting O if NULL
  GEN xinv=alginv(A, x);
  GEN Oconj=cgetg_copy(O, &n);
  for(long i=1;i<n;i++) gel(Oconj, i)=algmul(A, x, algmul(A, gel(O, i), xinv));
  return gerepilecopy(top, Oconj);
}

//Returns the level of the order O.
GEN
algorderlevel(GEN A, GEN O, int factored)
{
  pari_sp top=avma;
  GEN F=alg_get_center(A);
  GEN odisc=algorderdisc(A, O, 0, 0);
  GEN adisc=idealpow(F, algreduceddisc(A), gen_2);
  GEN levelsqr=idealdivexact(F, odisc, adisc);//Divides exactly
  GEN fact=idealfactor(F, levelsqr);
  for(long i=1;i<lg(gel(fact, 1));i++) gcoeff(fact, i, 2)=shifti(gcoeff(fact, i, 2), -1);//Divide by 2
  if(factored) return gerepilecopy(top, fact);
  return gerepileupto(top, idealfactorback(F, fact, NULL, 0));
}

//Given an order O in A, returns the discriminant of the order, which is disc(algebra)*level
GEN
algorderdisc(GEN A, GEN O, int reduced, int factored)
{
  pari_sp top=avma;
  GEN F=alg_get_center(A);//The centre
  GEN zdisc=sqri(mulii(ZM_det(O), algnormdisc(A)));//The norm to Q of the (non-reduced) discriminant
  GEN I=gen_0;//Stores the ideal generated by d(a1, a2, a3, a4) as we run across the basis.
  long lbas=lg(O);
  for(long i1=1;i1<lbas;i1++){
    for(long i2=i1+1;i2<lbas;i2++){
      for(long i3=i2+1;i3<lbas;i3++){
        for(long i4=1;i4<lbas;i4++){
          GEN d=algd(A, mkvec4(gel(O, i1), gel(O, i2), gel(O, i3), gel(O, i4)));
          I=idealadd(F, I, d);
          if(equalii(idealnorm(F, I), zdisc)){//We are done!
            if(!reduced && !factored) return gerepilecopy(top, I);
            GEN fact=idealfactor(F, I);
            if(!reduced && factored) return gerepileupto(top, fact);//From now on, we want it reduced
            for(long i=1;i<lg(gel(fact, 1));i++) gcoeff(fact, i, 2)=shifti(gcoeff(fact, i, 2), -1);//Divide by 2
            if(factored) return gerepilecopy(top, fact);
            return gerepileupto(top, idealfactorback(F, fact, NULL, 0));
          }
        }
      }
    }
  }
  pari_warn(warner, "We did not succeed in finding the discriminant.");
  set_avma(top);
  return gen_0;
}

//Returns the reduced discriminant of A as an ideal in the centre.
GEN
algreduceddisc(GEN A)
{
  pari_sp top=avma;
  GEN F=alg_get_center(A);
  GEN rplaces=algramifiedplacesf(A);
  long l=lg(rplaces);
  GEN e=const_vecsmall(l-1, 1);//Vecsmall of 1's
  return gerepileupto(top, idealfactorback(F, rplaces, e, 0));
}

//Returns a quaternion algebra over F (of degree n) with |N_{F/Q}(discriminant)|=D and split at the infinite place place only, if this exists. We also guarantee that a>0. F must be a totally real field.
GEN
algshimura(GEN F, GEN D, long place, long maxcomptime, int allowswap)
{
  pari_sp top=avma;
  if(nf_get_r2(F)>0) return gen_0;//Not totally real!
  long n=nf_get_degree(F);
  if(place<=0 || place>n) return gen_0;//Invalid place!
  GEN infram=cgetg(n+1, t_VEC);
  for(long i=1;i<=n;i++) gel(infram, i)=gen_1;
  gel(infram, place)=gen_0;
  GEN A;
  
  if(maxcomptime) pari_alarm(maxcomptime);
  pari_CATCH(CATCH_ALL){
    pari_warn(warner, "Time limit or memory exceeded, skipping.");
    set_avma(top);
    pari_CATCH_reset();
    return gen_0;
  }
  pari_TRY{
    A=algfromnormdisc(F, D, infram);
    if(gequal0(A)){//Nope
      set_avma(top);
      pari_CATCH_reset();
      return gen_0;
    }
    GEN L=alg_get_splittingfield(A), pol=rnf_get_pol(L);//L=K(sqrt(a))
    long varn=rnf_get_varn(L);
    GEN a=gneg(gsubst(pol, varn, gen_0));//Polynomial is of the form x^2-a, so plug in 0 and negate to get a
    if(gsigne(gsubst(a, nf_get_varn(F), gel(nf_get_roots(F), place)))!=1){//Must swap a, b
      if(!allowswap){set_avma(top);A=gen_0;}//We do not allow the swap.
      else{//Allowing the swap
        GEN b=lift(alg_get_b(A));
        GEN aden=Q_denom(a), bden=Q_denom(b);
        if(!isint1(aden)) a=gmul(a, gsqr(aden));//We need to get rid of the denominator of a
        if(!isint1(bden)) b=gmul(b, gsqr(bden));//We need to get rid of the denominator of b
        A=gerepileupto(top, alginit(F, mkvec2(b, a), -1, 1));//Swapping a, b
      }
    }
    else A=gerepilecopy(top, A);//No swap required.
  }pari_ENDCATCH
  pari_alarm(0);//Stop the alarm!
  return A;
}

//Returns a quaternion algebra over F (of degree n) with |N_{F/Q}(discriminant)|=D and split at the infinite place place only, if this exists. We also guarantee that a>0. F must be a totally real field. This returns the pair [a, b] that make the algebra. If allowswap=0, we do NOT allow the swapping of the a, b output by alginit, and either return the algebra (if it is valid) or 0 otherwise. Naively swapping a, b when deg(F) is large can be far too slow.
GEN
algshimura_ab(GEN F, GEN D, long place, long maxcomptime, int allowswap)
{
  pari_sp top=avma;
  if(nf_get_r2(F)>0) return gen_0;//Not totally real!
  long n=nf_get_degree(F);
  if(place<=0 || place>n) return gen_0;//Invalid place!
  GEN infram=cgetg(n+1, t_VEC);
  for(long i=1;i<=n;i++) gel(infram, i)=gen_1;
  gel(infram, place)=gen_0;
  GEN A;
  
  if(maxcomptime) pari_alarm(maxcomptime);
  pari_CATCH(CATCH_ALL){
    pari_warn(warner, "Time limit or memory exceeded, skipping.");
    set_avma(top);
    pari_CATCH_reset();
    return gen_0;
  }
  pari_TRY{
    A=algfromnormdisc(F, D, infram);
    if(gequal0(A)){//Nope
      set_avma(top);
      pari_CATCH_reset();
      return gen_0;
    }
    GEN L=alg_get_splittingfield(A), pol=rnf_get_pol(L);//L=K(sqrt(a))
    long varn=rnf_get_varn(L);
    GEN a=gneg(gsubst(pol, varn, gen_0));//Polynomial is of the form x^2-a, so plug in 0 and negate to get a
    GEN b=lift(alg_get_b(A));
    GEN aden=Q_denom(a), bden=Q_denom(b);//When initializing the algebra, PARI insists that a, b have no denominators
    if(!isint1(aden)) a=gmul(a, gsqr(aden));//We need to get rid of the denominator of a
    if(!isint1(bden)) b=gmul(b, gsqr(bden));//We need to get rid of the denominator of b
    if(gsigne(gsubst(a, nf_get_varn(F), gel(nf_get_roots(F), place)))!=1){
      if(allowswap) A=gerepilecopy(top, mkvec2(b, a));//Must swap a, b
      else{
        set_avma(top);//If we do not allow the swapping of a, b (recommended if deg(F)>=6), we do not return anything.
        A=gen_0;
      }
    }
    else A=gerepilecopy(top, mkvec2(a, b));//No swap required.
  }pari_ENDCATCH
  pari_alarm(0);//Stop the alarm!
  return A;
}

//Returns the algebra where a, b are swapped
GEN
algswapab(GEN A)
{
  pari_sp top=avma;
  GEN L=alg_get_splittingfield(A), pol=rnf_get_pol(L);//L=K(sqrt(a))
  long varn=rnf_get_varn(L);
  GEN a=gneg(gsubst(pol, varn, gen_0));//Polynomial is of the form x^2-a, so plug in 0 and negate to get a
  GEN b=lift(alg_get_b(A));
  GEN aden=Q_denom(a), bden=Q_denom(b);
  if(!isint1(aden)){
    a=gmul(a, gsqr(aden));
    pari_warn(warner, "b has to have no denominator, so we scaled it by its denominator squared.");
  }
  if(!isint1(bden)){
    b=gmul(b, gsqr(bden));
    pari_warn(warner, "a has to have no denominator, so we scaled it by its denominator squared.");
  }
  return gerepileupto(top, alginit(alg_get_center(A), mkvec2(b, a), -1, 1));
}

//Returns at most nwant pairs [a, b] for quaternion algebras over the field F that are suitable for fundamental domains. We start at N_F/Q(disc)=Dmin if supplied, and 2 otherwise. We stop at Dmax if supplied as non-zero, and otherwise continue until we fill up all nwant. Returns [discs, pairs]
GEN
smallalgebras(GEN F, long nwant, GEN Dmin, GEN Dmax, long maxcomptime, int allowswap)
{
  pari_sp top=avma;
  long nwp1=nwant+1;
  GEN Ds=vectrunc_init(nwp1);
  GEN v=vectrunc_init(nwp1);
  if(gequal0(Dmax)) Dmax=mkoo();
  long ind=1;
  GEN D=gceil(gmax(Dmin, gen_2));//Starting D
  while(ind<=nwant && gcmp(D, Dmax)<=0){
    GEN A=algshimura_ab(F, D, 1, maxcomptime, allowswap);
    if(!gequal0(A)){
      vectrunc_append(Ds, D);
      vectrunc_append(v, A);
      ind++;
    }
    D=gaddgs(D, 1);
  }
  return gerepilecopy(top, mkvec2(Ds, v));
}

//smallalgebras, but we provide a min and max area we want to find examples with.
GEN
smallalgebras_area(GEN nf, GEN Amin, GEN Amax, int retD, int maxcomptime, int allowswap)
{
  pari_sp top=avma, mid;
  long bits=bit_accuracy(3);
  GEN pol=nf_get_pol(nf);
  GEN zetaval=lfun(pol, gen_2, bits);//zeta_F(2)
  long n=nf_get_degree(nf), twon=2*n;
  GEN ar=gmul(gpow(nfdisc(pol), gdivsg(3, gen_2), 3), zetaval);//zeta_F(2)*d_F^(3/2)
  ar=gmul(ar, gpowgs(gen_2, 3-twon));
  ar=gmul(ar, gpowgs(mppi(3), 1-twon));//2^(3-2n)*Pi^(1-2n)*zeta_F(2)*d_F^(3/2). The area is ar*phi(D), where phi(p^e)=p^e-1 and phi is multiplicative
  GEN pDmin=gdiv(Amin, ar), pDmax=gdiv(Amax, ar);
  mid=avma;
  pDmin=gceil(pDmin);
  pDmax=gceil(pDmax);//The lower and upper bounds for phi(D). pDmax is actually the 1+the max
  gerepileallsp(top, mid, 2, &pDmin, &pDmax);
  GEN plist=primes0(mkvec2(gen_2, pDmax));//The possible prime divisors of D.
  long lp=lg(plist);
  GEN pposs=vectrunc_init(n*lp);
  for(long i=1;i<lp;i++){
    GEN pdec=idealprimedec(nf, gel(plist, i));
    for(long j=1;j<lg(pdec);j++){
      GEN nm=pr_norm(gel(pdec, j));
      if(gcmp(nm, pDmax)<=0) vectrunc_append(pposs, nm);//Adding the norms.
    }
  }
  pposs=gen_sort_uniq(pposs, (void*)cmpii, cmp_nodata);//Remove duplicates of the possible norms.
  GEN pposs_nm=cgetg_copy(pposs, &lp);
  for(long i=1;i<lp;i++) gel(pposs_nm, i)=subis(gel(pposs, i), 1);//The norms
  long par=(n-1)%2, ind=0;//The parity of the number of prime divisors we need.
  GEN pro=gen_1;
  while(cmpii(pro, pDmax)==-1 && ind<lp){//ind keeps track of the maximal number of prime divisors we can take
    ind++;
    pro=mulii(pro, gel(pposs_nm, ind));
  }
  if(ind%2!=par) ind--;//Make it line up with parity.
  glist *Ds=NULL;
  long count=0;
  if(par==0){
    if(cmpis(pDmin, 1)<=0) glist_putstart(&Ds, gen_1);//Disc 1
    par=2;
  }
  long lpm1=lp-1;//Number of prime powers we are looping over
  for(long ndiv=par;ndiv<=ind;ndiv=ndiv+2){//ndiv divisors
    GEN sub=cgetg(ndiv+1 ,t_VECSMALL);
    for(long i=1;i<=ndiv;i++) sub[i]=i;
    int cont=1;
    while(cont){
      mid=avma;
      GEN prod=gen_1;
      for(long i=1;i<=ndiv;i++) prod=mulii(prod, gel(pposs_nm, sub[i]));//Computing the norm.
      if(cmpii(prod, pDmin)==-1) set_avma(mid);//Too small, just go on.
      else{
        if(cmpii(prod, pDmax)==1){//Too large
          set_avma(mid);
          long nind=0;
          for(long i=ndiv;i>=2;i--){if(sub[i]>sub[i-1]+1){nind=i-1;break;}}
          if(nind==0) break;//Done, we will always be too large for now on out.
          long tochange=sub[nind]+1;//What we want to start with
          for(long i=nind;i<=ndiv;i++){sub[i]=tochange;tochange++;}
          continue;
        }
        else{
          prod=gen_1;
          for(long i=1;i<lg(sub);i++) prod=mulii(prod, gel(pposs, sub[i]));//Computing the D
          glist_putstart(&Ds, prod);
          count++;
        }
      }
      cont=subset_lexincrem(sub, lpm1);//Increment the subset
    }
  }
  GEN Discs=gerepileupto(top, sort(glist_togvec(Ds, count, -1)));//The discriminants.
  if(retD) return Discs;
  long lD=lg(Discs), sind=1;
  GEN Alist=vectrunc_init(lD), Dlist=vectrunc_init(lD);
  if(equali1(gel(Discs, 1))) sind=2;//We don't allow D=1, since this does not work with algshimura_ab.
  for(long i=sind;i<lD;i++){
    GEN A=algshimura_ab(nf, gel(Discs, i), 1, maxcomptime, allowswap);
    if(!gequal0(A)){
      vectrunc_append(Dlist, gel(Discs, i));
      vectrunc_append(Alist, A);
    }
  }
  return gerepilecopy(top, mkvec2(Dlist, Alist));
}

//S is a k-element subset of [1,...,n], given as a vecsmall. This updates S to the next lexicographical element. Returns 1 if successful, 0 if at the end.
static int
subset_lexincrem(GEN S, long n)
{
  long k=lg(S)-1, bound=n, ind=0;
  for(long i=k;i>0;i--){
    if(S[i]<bound){ind=i;S[i]=S[i]+1;break;}//Success
    bound--;
  }
  if(ind){//Success, fill in rest
    bound=S[ind]+1;
    for(long i=ind+1;i<=k;i++){S[i]=bound;bound++;}
    return 1;
  }
  else return 0;//No success
}


//Returns small norm nm elements (Q_{z1,z2}(x)<=C) of the order O
GEN
algsmallelts(GEN A, GEN O, GEN nm, GEN p, GEN C, GEN z1, GEN z2, long prec)
{
  pari_sp top=avma;
  GEN Q;
  if(!O) Q=qalg_fdominitialize(A, NULL, NULL, prec);//Maximal order in A
  else Q=qalg_fdominitialize(A, gel(O, 1), gel(O, 2), prec);//Supplied Eichler order
  GEN nf=alg_get_center(A);
  long nfdeg=nf_get_degree(nf), fourn=4*nfdeg;
  GEN nformpart=qalg_normform(Q);
  GEN nfdecomp=mat_nfcholesky(nf, nformpart);
  for(long i=1;i<=fourn;i++){
    for(long j=1;j<=fourn;j++){
      gcoeff(nformpart, i, j)=nftrace(nf, gcoeff(nformpart, i, j));//Taking the trace to Q
    }
  }//Tr_{nf/Q}(nrd(elt));
  return gerepileupto(top, qalg_smallelts_qfminim(Q, nm, p, C, z1, z2, 0, nfdecomp, nformpart, prec));
}

//Computes all norm nm elements for which Q_{z_1,z_2}(x)<=C. maxelts is the maximum number of return elements (or 0 for all norm 1 elts)
GEN
qalg_smallelts_qfminim(GEN Q, GEN nm, GEN p, GEN C, GEN z1, GEN z2, long maxelts, GEN nfdecomp, GEN nformpart, long prec)
{
  pari_sp top=avma, mid;
  GEN A=qalg_get_alg(Q);
  GEN nf=alg_get_center(A);
  GEN mats=psltopsu_transmats(p);
  GEN absrednorm=qalg_absrednormqf(Q, mats, z1, z2, nformpart, prec);
  GEN vposs=gel(qfminim0(absrednorm, C, NULL, 2, prec), 3);
  absrednorm=RgM_gtofp(absrednorm, prec);//absrednorm can have integer/rational entries, we want to make sure all entries are t_REAL, as otherwise it's a little slower, and a segmentation fault can occur in rare instances.
  GEN O=qalg_get_order(Q);
  int nonmax=0;
  if(!gequal1(qalg_get_level(Q))) nonmax=1;
  long nvposs=lg(vposs), mret;
  if(maxelts) mret=maxelts+1;
  else mret=nvposs;
  GEN ret=vectrunc_init(mret), norm;
  for(long i=1;i<nvposs;i++){
    mid=avma;
    norm=algnorm_chol(nf, nfdecomp, gel(vposs, i));
    if(gequal(norm, nm)){
      set_avma(mid);
      if(nonmax) vectrunc_append(ret, QM_QC_mul(O, gel(vposs, i)));//Change of basis backwards
      else vectrunc_append(ret, gel(vposs, i));//Don't append a copy, will copy at the end.
      if(lg(ret)>=mret) break;//Done
      continue;
    }
    set_avma(mid);
  }
  return gerepilecopy(top, ret);
}



//SECTION 3: PRODUCING DATA FOR MY PAPER


//3: OPTIMIZING THE VALUE OF C FOR ENUMERATION


//Finds and returns the optimal C for a given algebra. Actually, returns [A, B, C, R^2]. We do ntrials trials in a range [Cmin, scale^(1/2n)*Cmin], where this is centred at the conjectural best C value.
GEN
enum_bestC(GEN A, GEN p, GEN scale, long ntrials, long mintesttime, long prec)
{
  pari_sp top=avma;
  GEN Cbest=algfdom_bestC(A, NULL, prec);
  GEN nf=alg_get_center(A);
  long n=nf_get_degree(nf), twon=2*n, fourn=2*twon;
  GEN scalefact=gpow(scale, gdivgs(gen_1, fourn), prec);
  GEN Cmin=gdiv(Cbest, scalefact);
  GEN Cmax=gmul(Cbest, scalefact);
  if(ntrials<=1) ntrials=2;//Make sure it's at least 2.
  GEN blen=gdivgs(gsub(Cmax, Cmin), ntrials-1);//length.
  GEN Clist=cgetg(ntrials+1, t_VEC);
  GEN C=Cmin;
  gel(Clist, 1)=C;
  for(long i=2;i<=ntrials;i++){
    C=gadd(C, blen);
    gel(Clist, i)=C;
  }
  GEN times=enum_time(A, p, Clist, mintesttime, prec);//Executing the computation
  GEN Xdata=cgetg(ntrials+1, t_MAT);
  for(long i=1;i<=ntrials;i++) gel(Xdata, i)=mkcol2(gen_1, gpowgs(gel(Clist, i), twon));//The X data
  GEN reg=OLS(Xdata, times, 1);
  GEN Aco=gmael(reg, 1, 1), Bco=gmael(reg, 1, 2);
  if(gsigne(Aco)==-1 || gsigne(Bco)==-1) pari_err(e_MISC, "One of A, B are negative! This is not correct.");
  C=enum_bestC_givenAB(n, Aco, Bco, prec);
  return gerepilecopy(top, mkvec4(Aco, Bco, C, gel(reg, 2)));//[A, B, C, R^2].
}

//Returns the unique real solution >n to (2n-1)C^2n-2n^2c^(2n-1)=A/B.
static GEN
enum_bestC_givenAB(long n, GEN A, GEN B, long prec)
{
  pari_sp top=avma;
  long twon=2*n;
  long var=fetch_var();//The temporary variable number
  GEN P=cgetg(twon+3, t_POL);
  P[1]=evalvarn(var);
  gel(P, 2)=gdiv(gneg(A), B);//Constant coefficient.
  for(long i=3;i<=twon;i++) gel(P, i)=gen_0;//0's in the middle. P[i] corresp to x^(i-2)
  gel(P, twon+1)=stoi(-twon*n);
  gel(P, twon+2)=stoi(twon-1);
  P=normalizepol(P);//Normalize it in place.
  GEN rts=realroots(P, mkvec2(stoi(n), mkoo()), prec);//Real roots in [n, oo)
  delete_var();//Delete the variable.
  if(lg(rts)!=2) pari_err(e_MISC, "Could not find exactly one real root in [n, oo)!");
  return gerepilecopy(top, gel(rts, 1));
}

//Does enum_bestC for all algebras in Aset. If isArange, we assume that the base fields are all the same, and are regressing along Nm(disc)^(1/n). Otherwise, we assume n is fixed, and are regressing along disc(F).
GEN
enum_bestC_range(GEN Aset, GEN p, GEN scale, long ntrials, long mintesttime, char *fname, int isArange, int compile, int WSL, long prec)
{
  pari_sp top=avma;
  if(!pari_is_dir("plots/build")){//Checking the directory
    int s=system("mkdir -p plots/build");
    if(s==-1) pari_err(e_MISC, "ERROR CREATING DIRECTORY");
  }
  char *fname_full=pari_sprintf("plots/build/%s.dat", fname);
  FILE *f=fopen(fname_full, "w");
  pari_free(fname_full);
  pari_fprintf(f, "x y rsqr\n");

  GEN nf=alg_get_center(gel(Aset, 1));
  long n=nf_get_degree(nf);
  GEN oneovertwon=gdivgs(gen_1, 2*n), oneovern=gmulgs(oneovertwon, 2);
  GEN nfdisc;
  if(isArange) nfdisc=nf_get_disc(nf);//This is constant across all algebras when isArange=1.
  long lgAset=lg(Aset);
  GEN Xdat=vectrunc_init(lgAset), Cdat=coltrunc_init(lgAset), lastX=gen_0, firstX=NULL;
  long triespertrial=3;
  for(long i=1;i<lgAset;i++){
    GEN A=gel(Aset, i);
    GEN C=NULL;
    for(long trial=1;trial<=triespertrial;trial++){
      pari_CATCH(CATCH_ALL){
        C=NULL;
        pari_printf("Error in algebra %d trial %d, retrying\n", i, trial);
      }
      pari_TRY{C=enum_bestC(A, p, scale, ntrials, mintesttime, prec);}pari_ENDCATCH
      if(C) break;
    }
    if(!C){
      pari_printf("Algebra %d skipped due to error.\n", i);
      continue;
    }
    if(isArange){//Varying Adisc
      GEN Adisc=algnormdisc(A);
      pari_fprintf(f, "%Pd %Pf %Pf\n", Adisc, gel(C, 3), gel(C, 4));
      if(!firstX) firstX=Adisc;
      lastX=Adisc;
      vectrunc_append(Xdat, gpow(Adisc, oneovertwon, prec));
    }
    else{//Varying F
      nf=alg_get_center(A);
      nfdisc=nf_get_disc(nf);
      GEN Adisc=algnormdisc(A);
      gel(C, 3)=gdiv(gel(C, 3), gpow(Adisc, oneovertwon, prec));
      pari_fprintf(f, "%Pd %Pf %Pf\n", nfdisc, gel(C, 3), gel(C, 4));
      if(!firstX) firstX=nfdisc;
      lastX=nfdisc;
      vectrunc_append(Xdat, gpow(nfdisc, oneovern, prec));
    }
    vectrunc_append(Cdat, gel(C, 3));
    pari_printf("Algebra %d done.\n", i);
  }
  fclose(f);
  GEN reg=OLS_nointercept(Xdat, Cdat, 1);
  enum_bestC_plot(reg, firstX, lastX, n, fname, isArange);
  if(compile) plot_compile(fname, WSL);
  return gerepilecopy(top, reg);
}

//Prepares a basic latex plot of the data.
static void
enum_bestC_plot(GEN reg, GEN Cmin, GEN Cmax, long n, char *fdata, int isArange)
{
  pari_sp top=avma;
  if(!pari_is_dir("plots/build")){//Checking the directory
      int s=system("mkdir -p plots/build");
      if(s==-1) pari_err(e_MISC, "ERROR CREATING DIRECTORY");
  }
  char *plotmake=pari_sprintf("plots/build/%s_plotter.tex", fdata);
  FILE *f=fopen(plotmake, "w");
  pari_free(plotmake);
  long invpower;
  if(isArange) invpower=2*n;
  else invpower=n;
  pari_fprintf(f, "\\documentclass{article}\n\\usepackage{amsmath, amssymb, pgfplots}\n  \\usepgfplotslibrary{external}\n  \\tikzexternalize\n");
  pari_fprintf(f, "  \\pgfplotsset{compat=1.16}\n\\begin{document}\n\\tikzsetnextfilename{%s}\n\\begin{tikzpicture}\n  \\begin{axis}", fdata);
  if(isArange) pari_fprintf(f, "[xlabel=$\\text{Nm}_{F/\\mathbb{Q}}(\\mathfrak{D})$, ylabel=C,\n");
  else pari_fprintf(f, "[xlabel=$\\text{disc}(F)$, ylabel=$\\text{C}/\\text{Nm}_{F/\\mathbb{Q}}(\\mathfrak{D})^{1/%d}$,\n", 2*invpower);
  pari_fprintf(f, "    xmin=0, ymin=0,\n");
  pari_fprintf(f, "    scatter/classes={a={mark=o}}, clip mode=individual,]\n");
  pari_fprintf(f, "    \\addplot[scatter, blue, only marks, mark size=0.9]\n      table[x=x,y=y,col sep=space]{%s.dat};\n", fdata);
  pari_fprintf(f, "    \\addplot[red, ultra thick, samples=1000, domain=%Pf:%Pf]{%Pf*(x)^(1/%d)}", Cmin, Cmax, gel(reg, 1), invpower);
  pari_fprintf(f, ";%%R^2=%Pf\n  \\end{axis}\n\\end{tikzpicture}\n\\end{document}", gel(reg, 2));
  fclose(f);
  set_avma(top);
}

//Returns the number of non-trivial elements. Must be in the basis representation
static long
enum_nontrivial(GEN L)
{
  long nt=0;
  for(long i=1;i<lg(L);i++){
    if(!gequal1(gmael(L, i, 1)) && !gequalm1(gmael(L, i, 1))){nt++;continue;}
    for(long j=2;j<lg(gel(L, i));j++) if(!gequal0(gmael(L, i, j))){nt++;break;}
  }
  return nt;
}

//Does smallnorm1elts Ntests times, and sees how many successes we get. Returns [obs, exp], with obs the number of successes, and exp the expected number (2*Pi(C-n)/area)
GEN
enum_successrate(GEN A, GEN p, GEN C, long Ntests, GEN R, long prec)
{
  pari_sp top=avma;
  GEN area=algfdomarea(A, NULL, 1, prec);
  if(gequal0(R)){
    GEN gamma=dbltor(2.1);
    R=hdiscradius(gpow(area, gamma, prec), prec);
  }
  GEN Q=qalg_fdominitialize(A, NULL, NULL, prec);
  GEN nf=alg_get_center(A);
  GEN nformpart=qalg_normform(Q);
  GEN normdecomp=mat_nfcholesky(nf, nformpart);
  long n=4*nf_get_degree(nf);//The lg of a normal entry
  for(long i=1;i<=n;i++){
    for(long j=1;j<=n;j++){
      gcoeff(nformpart, i, j)=nftrace(nf, gcoeff(nformpart, i, j));//Taking the trace to Q
    }
  }//Tr_{K/Q}(nrd(elt));
  GEN retdat=enum_successrate_givendata(Q, p, C, Ntests, R, area, normdecomp, nformpart, prec);
  if(!retdat) pari_err_PREC("Precision is too low, please increase");
  return gerepileupto(top, retdat);
}

//enum_successrate, but over a range of C. Prints the observed data to file plots/build/fname.dat (if not NULL), and returns [A, B, R^2]: the predicted heuristic of A+BC and the R^2 value of this heuristic (B=2Pi/area, A=-2*pi*n/area)
GEN
enum_successrate_range(GEN A, GEN p, GEN Cmin, GEN Cmax, long ntrials, long Ntests, GEN R, char *fname, int compile, int WSL, long prec)
{
  pari_sp top=avma;
  GEN area=algfdomarea(A, NULL, 1, prec);//Initialize things
  if(gequal0(R)){
    GEN gamma=dbltor(2.1);
    R=hdiscradius(gpow(area, gamma, prec), prec);
  }
  GEN Q=qalg_fdominitialize(A, NULL, NULL, prec);
  GEN nf=alg_get_center(A);
  GEN nformpart=qalg_normform(Q);
  GEN normdecomp=mat_nfcholesky(nf, nformpart);
  long n=4*nf_get_degree(nf);//The lg of a normal entry
  for(long i=1;i<=n;i++){
    for(long j=1;j<=n;j++){
      gcoeff(nformpart, i, j)=nftrace(nf, gcoeff(nformpart, i, j));//Taking the trace to Q
    }
  }//Tr_{K/Q}(nrd(elt));
  //On to data collection
  if(ntrials<=1) ntrials=2;
  long ntp1=ntrials+1;
  GEN found=cgetg(ntp1, t_COL);
  GEN Cmat=cgetg(ntp1, t_MAT);
  GEN C=Cmin;
  GEN blen=gdivgs(gsub(Cmax, Cmin), ntrials-1);
  for(long i=1;i<=ntrials;i++){//Each trial
    GEN getdat=enum_successrate_givendata(Q, p, C, Ntests, R, area, normdecomp, nformpart, prec);
    if(!getdat) pari_err_PREC("Precision is too low, please increase");
    gel(found, i)=gel(getdat, 1);
    gel(Cmat, i)=mkcol2(gen_1, C);
    C=gadd(C, blen);
    if(i%10==0) pari_printf("Trial %d done.\n", i);
  }
  if(fname){
    if(!pari_is_dir("plots/build")){//Checking the directory
      int s=system("mkdir -p plots/build");
      if(s==-1) pari_err(e_MISC, "ERROR CREATING DIRECTORY");
    }
    char *fname_full=pari_sprintf("plots/build/%s.dat", fname);
    FILE *f=fopen(fname_full, "w");
    pari_free(fname_full);
    pari_fprintf(f, "x y\n");
    for(long i=1;i<=ntrials;i++) pari_fprintf(f, "%Pf %Pd\n", gcoeff(Cmat, 2, i), gel(found, i));
    fclose(f);
  }
  //R^2 time!
  GEN slope=gmulgs(gdiv(Pi2n(1, prec), area), Ntests);
  GEN constant=gmulgs(slope, -nf_get_degree(nf));
  GEN fitdat=mkcol2(constant, slope);
  GEN rsqr=rsquared(Cmat, found, fitdat);
  GEN dat=mkvec3(constant, slope, rsqr);
  if(compile && fname){
    enum_successrate_plot(A, dat, Cmin, Cmax, fname, WSL);
    plot_compile(fname, WSL);
  }
  return gerepilecopy(top, dat);
}

//enum_successrate with the data initialized.
static GEN
enum_successrate_givendata(GEN Q, GEN p, GEN C, long Ntests, GEN R, GEN area, GEN normdecomp, GEN nformpart, long prec)
{
  pari_sp top=avma, mid;
  GEN nf=alg_get_center(qalg_get_alg(Q));
  long found=0;
  for(long test=1;test<=Ntests;test++){
    mid=avma;
    GEN z=randompoint_ud(R, prec);
    GEN elts=qalg_smallnorm1elts_qfminim(Q, p, C, gen_0, z, 0, normdecomp, nformpart, prec);
    if(!elts) return gc_NULL(top);//Precision too low
    if(lg(elts)>1) found=found+enum_nontrivial(elts);
    set_avma(mid);
  }
  GEN expect=gdiv(gmul(gmulgs(Pi2n(1, prec), Ntests), gsubgs(C, nf_get_degree(nf))), area);
  return gerepilecopy(top, mkvec2(stoi(found), expect));
}

//Prepares a basic latex plot of the data.
static void
enum_successrate_plot(GEN A, GEN reg, GEN Cmin, GEN Cmax, char *fdata, int WSL)
{
  pari_sp top=avma;
  if(!pari_is_dir("plots/build")){//Checking the directory
      int s=system("mkdir -p plots/build");
      if(s==-1) pari_err(e_MISC, "ERROR CREATING DIRECTORY");
  }
  char *plotmake=pari_sprintf("plots/build/%s_plotter.tex", fdata);
  FILE *f=fopen(plotmake, "w");
  pari_free(plotmake);
  pari_fprintf(f, "\\documentclass{article}\n\\usepackage{pgfplots}\n  \\usepgfplotslibrary{external}\n  \\tikzexternalize\n");
  pari_fprintf(f, "  \\pgfplotsset{compat=1.16}\n\\begin{document}\n\\tikzsetnextfilename{%s}\n\\begin{tikzpicture}\n  \\begin{axis}", fdata);
  pari_fprintf(f, "[xlabel=C, ylabel=Elements found,\n");
  pari_fprintf(f, "    xmin=%Pf, xmax=%Pf, ymin=0,\n", gsubgs(Cmin, 1), gaddgs(Cmax, 1));
  pari_fprintf(f, "    scatter/classes={a={mark=o}}, clip mode=individual,]\n");
  pari_fprintf(f, "    \\addplot[scatter, blue, only marks, mark size=0.9]\n      table[x=x,y=y,col sep=space]{%s.dat};\n", fdata);
  pari_fprintf(f, "    \\addplot[red, ultra thick, samples=1000, domain=%Pf:%Pf] (x, %Pf+%Pf*x", Cmin, Cmax, gel(reg, 1), gel(reg, 2));
  pari_fprintf(f, ");%%R^2=%Pf\n  \\end{axis}\n\\end{tikzpicture}\n\\end{document}", gel(reg, 3));
  fclose(f);
  set_avma(top);
}

//Computes the average time to find algsmallnormelts(A, C, 0, z) for all C in Cset, and returns it as a column vector.
GEN
enum_time(GEN A, GEN p, GEN Cset, long mintesttime, long prec)
{
  pari_sp top=avma, mid;
  GEN Q=qalg_fdominitialize(A, NULL, NULL, prec);
  GEN nf=alg_get_center(A);
  GEN normformpart=qalg_normform(Q);
  GEN nfdecomp=mat_nfcholesky(nf, normformpart);
  long n=lg(gel(alg_get_basis(A), 1));
  for(long i=1;i<n;i++){
    for(long j=1;j<n;j++){
      gcoeff(normformpart, i, j)=nftrace(nf, gcoeff(normformpart, i, j));//Taking the trace to Q
    }
  }//Tr_{K/Q}(nrd(elt));
  GEN area=algfdomarea(A, NULL, 1, prec);
  GEN R=hdiscradius(gpow(area, dbltor(2.1), prec), prec);
  
  long lgCset=lg(Cset);
  GEN avgtimes=cgetg(lgCset, t_COL);
  pari_timer T;
  timer_start(&T);
  for(long i=1;i<lgCset;i++){
    long t=0, tries=0;
    timer_delay(&T);
    while(t<mintesttime){//Make sure we do at least mintesttime per test
      mid=avma;
      tries++;
      GEN z=randompoint_ud(R, prec);
      GEN theelts=qalg_smallnorm1elts_qfminim(Q, p, gel(Cset, i), gen_0, z, 0, nfdecomp, normformpart, prec);
      if(!theelts) return gc_NULL(top);
      t=timer_get(&T);
      set_avma(mid);
    }
    gel(avgtimes, i)=rdivss(t, 1000*tries, prec);
    if(i%10==0) pari_printf("%d test cases done.\n", i);
  }
  return gerepilecopy(top, avgtimes);
}

//enum_time, but we use ntrials>=2 from Cmin to Cmax. This also performs the regression on the data, and returns the result. We store the data in plots/build/fdata.dat, and compile the latex document to view it if this is not NULL.
GEN
enum_time_range(GEN A, GEN p, GEN Cmin, GEN Cmax, long ntrials, long mintesttime, char *fdata, int compile, int WSL, long prec)
{
  pari_sp top=avma;
  if(ntrials<=1) ntrials=2;
  GEN Clist=cgetg(ntrials+1, t_VEC);
  GEN C=Cmin;
  GEN blen=gdivgs(gsub(Cmax, Cmin), ntrials-1);
  for(long i=1;i<=ntrials;i++){
    gel(Clist, i)=C;
    C=gadd(C, blen);
  }
  GEN times=enum_time(A, p, Clist, mintesttime, prec);//Executing the computation
  if(!times) pari_err_PREC("Precision is too low, please increase");
  if(fdata){
    if(!pari_is_dir("plots/build")){//Checking the directory
      int s=system("mkdir -p plots/build");
      if(s==-1) pari_err(e_MISC, "ERROR CREATING DIRECTORY");
    }
    char *fdata_full=pari_sprintf("plots/build/%s.dat", fdata);
    FILE *f=fopen(fdata_full, "w");
    pari_free(fdata_full);
    pari_fprintf(f, "x y\n");
    for(long i=1;i<=ntrials;i++) pari_fprintf(f, "%Pf %Pf\n", gel(Clist, i), gel(times, i));
    fclose(f);
  }
  GEN nf=alg_get_center(A);
  long n=nf_get_degree(nf), twon=2*n;
  GEN Xdata=cgetg(ntrials+1, t_MAT);//Prep the regression
  for(long i=1;i<=ntrials;i++) gel(Xdata, i)=mkcol2(gen_1, gpowgs(gel(Clist, i), twon));//The X data
  GEN reg=gerepileupto(top, OLS(Xdata, times, 1));
  if(compile && fdata){
    enum_time_plot(A, reg, Cmin, Cmax, fdata);
    plot_compile(fdata, WSL);
  }
  return reg;
}

//Prepares a basic latex plot of the data.
static void
enum_time_plot(GEN A, GEN reg, GEN Cmin, GEN Cmax, char *fdata)
{
  pari_sp top=avma;
  if(!pari_is_dir("plots/build")){//Checking the directory
      int s=system("mkdir -p plots/build");
      if(s==-1) pari_err(e_MISC, "ERROR CREATING DIRECTORY");
  }
  GEN nf=alg_get_center(A);
  long n=nf_get_degree(nf);
  char *plotmake=pari_sprintf("plots/build/%s_plotter.tex", fdata);
  FILE *f=fopen(plotmake, "w");
  pari_free(plotmake);
  pari_fprintf(f, "\\documentclass{article}\n\\usepackage{pgfplots}\n  \\usepgfplotslibrary{external}\n  \\tikzexternalize\n");
  pari_fprintf(f, "  \\pgfplotsset{compat=1.16}\n\\begin{document}\n\\tikzsetnextfilename{%s}\n\\begin{tikzpicture}\n  \\begin{axis}", fdata);
  pari_fprintf(f, "[xlabel=C, ylabel=Time,\n");
  pari_fprintf(f, "    xmin=%Pf, xmax=%Pf, ymin=0,\n", gsubgs(Cmin, 1), gaddgs(Cmax, 1));
  pari_fprintf(f, "    scatter/classes={a={mark=o}}, clip mode=individual,]\n");
  pari_fprintf(f, "    \\addplot[scatter, blue, only marks, mark size=0.9]\n      table[x=x,y=y,col sep=space]{%s.dat};\n", fdata);
  pari_fprintf(f, "    \\addplot[red, ultra thick, samples=1000, domain=%Pf:%Pf] (x, %Pf+%Pf*x", Cmin, Cmax, gmael(reg, 1, 1), gmael(reg, 1, 2));
  for(long i=2;i<=2*n;i++) pari_fprintf(f, "*x");
  pari_fprintf(f, ");%%R^2=%Pf\n  \\end{axis}\n\\end{tikzpicture}\n\\end{document}", gel(reg, 2));
  fclose(f);
  set_avma(top);
}

//Returns the time taken to find nelts non-trivial elements
long
enum_timeforNelts(GEN A, GEN p, GEN C, long nelts, GEN R, int type, long prec)
{
  pari_sp top=avma;
  GEN Q=qalg_fdominitialize(A, NULL, NULL, prec);
  if(gequal0(R)){
    GEN area=algfdomarea(A, NULL, 1, prec);
    GEN gamma=dbltor(2.1);
    R=hdiscradius(gpow(area, gamma, prec), prec);
  }
  GEN nf=alg_get_center(A);
  long nfdeg=nf_get_degree(nf), fourn=4*nfdeg;
  GEN nfdecomp=gen_0, nform=gen_0;//nfdecomp used if type=1, nform if type=2
  GEN nformpart=qalg_normform(Q);
  if(type==1) nfdecomp=mat_nfcholesky(nf, nformpart);//qfminim
  else nform=gcopy(nformpart);//condition
  for(long i=1;i<=fourn;i++){
    for(long j=1;j<=fourn;j++){
      gcoeff(nformpart, i, j)=nftrace(nf, gcoeff(nformpart, i, j));//Taking the trace to Q
    }
  }//Tr_{nf/Q}(nrd(elt));
  long tottime=enum_timeforNelts_givendata(Q, p, C, nelts, R, type, nform, nfdecomp, nformpart, prec);
  if(tottime==-1) pari_err_PREC("Precision is too low, please increase");
  set_avma(top);
  return tottime;
}

//enum_timeforNelts given the setup data.
static long
enum_timeforNelts_givendata(GEN Q, GEN p, GEN C, long nelts, GEN R, int type, GEN nform, GEN nfdecomp, GEN nformpart, long prec)
{
  pari_sp top=avma, mid;
  GEN z, elts;
  long found=0, time;
  pari_timer T;
  timer_start(&T);
  if(type==1){//qfminim
    mid=avma;
    while(found<nelts){
      z=randompoint_ud(R, prec);//Random point
      elts=qalg_smallnorm1elts_qfminim(Q, p, C, gen_0, z, 1, nfdecomp, nformpart, prec);
      if(!elts) return gc_long(top, -1);
      if(lg(elts)!=1){
        if(enum_nontrivial(elts)>0) found++;
      }
      set_avma(mid);
    }
    time=timer_delay(&T);
  }
  else{//condition
    mid=avma;
    while(found<nelts){
      z=randompoint_ud(R, prec);//Random point
      elts=qalg_smallnorm1elts_condition(Q, p, C, gen_0, z, 1, nform, nformpart, prec);
      if(lg(elts)!=1){
        if(enum_nontrivial(elts)>0) found++;
      }
      set_avma(mid);
    }
    time=timer_delay(&T);
  }
  set_avma(top);
  return time;
}

//Computes enum_TimeforNelts for ntrials values of C between Cmin and Cmax, and outputs the results to data/fname.txt.
void
enum_timeforNelts_range(GEN A, GEN p, GEN Cmin, GEN Cmax, long ntrials, long nelts, char *fname, int type, int compile, int WSL, long prec)
{
  pari_sp top=avma;
  if(!pari_is_dir("plots/build")){//Checking the directory
    int s=system("mkdir -p plots/build");
    if(s==-1) pari_err(e_MISC, "ERROR CREATING DIRECTORY");
  }
  char *fname_full=pari_sprintf("plots/build/%s.dat", fname);
  FILE *f=fopen(fname_full, "w");
  pari_free(fname_full);
  if(ntrials<=1) ntrials=2;
  GEN Clist=cgetg(ntrials+1, t_VEC);
  GEN C=Cmin;
  GEN blen=gdivgs(gsub(Cmax, Cmin), ntrials-1);
  for(long i=1;i<=ntrials;i++){
    gel(Clist, i)=C;
    C=gadd(C, blen);
  }
  GEN Q=qalg_fdominitialize(A, NULL, NULL, prec);
  GEN area=algfdomarea(A, NULL, 1, prec);
  GEN gamma=dbltor(2.1);
  GEN R=hdiscradius(gpow(area, gamma, prec), prec);
  GEN nf=alg_get_center(A);
  long nfdeg=nf_get_degree(nf), fourn=4*nfdeg;
  GEN nfdecomp=gen_0, nform=gen_0;//nfdecomp used if type=1, nform if type=2
  GEN nformpart=qalg_normform(Q);
  if(type==1) nfdecomp=mat_nfcholesky(nf, nformpart);//qfminim
  else nform=gcopy(nformpart);//condition
  for(long i=1;i<=fourn;i++){
    for(long j=1;j<=fourn;j++){
      gcoeff(nformpart, i, j)=nftrace(nf, gcoeff(nformpart, i, j));//Taking the trace to Q
    }
  }//Tr_{nf/Q}(nrd(elt));
  GEN bestC=qalg_fdombestC(Q, prec);
  int past_bestC=0;
  double bestC_val=(double)0;
  pari_fprintf(f, "C t\n");
  for(long i=1;i<=ntrials;i++){
    C=gel(Clist, i);
    long t=enum_timeforNelts_givendata(Q, p, C, nelts, R, type, nform, nfdecomp, nformpart, prec);
    if(t==-1) pari_err_PREC("Precision is too low, please increase");
    double tins=t/((double)1000);
    pari_fprintf(f, "%Pf %.3f\n", C, tins);
    pari_printf("%Pf %.3f\n", C, tins);
    if(!past_bestC && gcmp(C, bestC)>=0){//We are past the optimal value of C.
      past_bestC=1;
      bestC_val=tins;
    }
  }
  fclose(f);
  if(compile){
    enum_timeforNelts_plot(A, bestC, bestC_val, fname);
    plot_compile(fname, WSL);
  }
  set_avma(top);
}

//Prepares a basic latex plot of the data.
static void
enum_timeforNelts_plot(GEN A, GEN bestC, double bestC_val, char *fname)
{
  pari_sp top=avma;
  if(!pari_is_dir("plots/build")){//Checking the directory
      int s=system("mkdir -p plots/build");
      if(s==-1) pari_err(e_MISC, "ERROR CREATING DIRECTORY");
  }
  char *plotmake=pari_sprintf("plots/build/%s_plotter.tex", fname);
  FILE *f=fopen(plotmake, "w");
  pari_free(plotmake);
  pari_fprintf(f, "\\documentclass{article}\n\\usepackage{pgfplots}\n  \\usepgfplotslibrary{external}\n  \\tikzexternalize\n");
  pari_fprintf(f, "  \\pgfplotsset{compat=1.16}\n\\begin{document}\n\\tikzsetnextfilename{%s}\n\\begin{tikzpicture}\n  \\begin{axis}", fname);
  pari_fprintf(f, "[xlabel=C, ylabel=t (s),\n");
  pari_fprintf(f, "    xmin=0, ymin=0,\n");
  pari_fprintf(f, "    scatter/classes={a={mark=o}}, clip mode=individual,]\n");
  pari_fprintf(f, "    \\addplot[scatter, blue, only marks, mark size=0.9]\n      table[x=C,y=t,col sep=space]{%s.dat};\n", fname);
  pari_fprintf(f, "    \\draw[dashed, red, ultra thick] (%Pf, 0) -- (%Pf, %.3f);\n", bestC, bestC, bestC_val);
  pari_fprintf(f, "  \\end{axis}\n\\end{tikzpicture}\n\\end{document}");
  fclose(f);
  set_avma(top);
}


//3: NUMBER OF ELEMENTS REQUIRED TO GENERATE ALGEBRA


//Initializes and checks the inputs, and computes the fundamental domain
GEN
algfdom_nelts(GEN A, GEN p, GEN CNRdata, int type, long prec)
{
  pari_sp top=avma;
  GEN tol=deftol(prec);
  GEN Q=qalg_fdominitialize(A, NULL, NULL, prec), newA=A, U;
  
  long newprec=prec;
  unsigned int precinc=0;
  for(;;){
    U=qalg_fdom_nelts(Q, p, CNRdata, type, tol, newprec);
    if(U) break;//Success, we have enough precision!
    //If we reach here, we need more precision.
    pari_warn(warner, "Increasing precision");
    newA=algmoreprec(newA, 1, newprec);//Increments the precision by 1.
    newprec++;
    precinc++;
    tol=deftol(newprec);
    Q=qalg_fdominitialize(newA, NULL, NULL, newprec);
  }
  if(precinc) pari_warn(warner, "Precision increased to %d, i.e. \\p%Pd. If U=output, then update the algebra to the correct precision with A=algfdomalg(U), and update the number field with F=algcenter(A). If the original values of a and b in A had denominators, then the new algebra will have cleared them (and hence have different a and b)", newprec, precision00(U, NULL));
  return gerepileupto(top, U);
}

//Computes the number of elements needed to generate the fundamental domain, and returns [nelts, #sides, area].
static GEN
qalg_fdom_nelts(GEN Q, GEN p, GEN CNRdata, int type, GEN tol, long prec)
{
  pari_sp top=avma, mid, bot;
  GEN mats=psltopsu_transmats(p);
  GEN A=qalg_get_alg(Q);
  GEN nf=alg_get_center(A);
  long nfdeg=nf_get_degree(nf);
  GEN area=qalg_fdomarea(Q, 3, prec);//Smallest precision possible.
  GEN areabound=gdivgs(gmulgs(area, 3), 2);//Times 1.5.
  
  GEN C, N, R, epsilon;//Constants used for bounds, can be auto-set or passed in.
  if(gequal0(CNRdata) || gequal0(gel(CNRdata, 1))) C=qalg_fdombestC(Q, prec);
  else C=gel(CNRdata, 1);
  if(gequal0(CNRdata) || gequal0(gel(CNRdata, 2))){//N
    GEN beta=gdivgs(gen_1, 10);
    N=gfloor(gadd(gmul(beta, area), gen_2));//We add 2 to make sure N>=2; errors can happen later otherwise
  }
  else N=gel(CNRdata, 2);
  if(gequal0(CNRdata) || gequal0(gel(CNRdata, 3))){//R
    GEN gamma=dbltor(2.1);//2.1
    R=hdiscradius(gpow(area, gamma, prec), prec);
  }
  else R=gel(CNRdata, 3);
  if(gequal0(CNRdata) || gequal0(gel(CNRdata, 4))) epsilon=gdivgs(gen_1, 6);
  else epsilon=gel(CNRdata, 4);
  
  GEN id=gel(alg_get_basis(A), 1);//The identity  
  long iN;
  GEN points, w;
  GEN U=cgetg(2, t_VEC);
  gel(U, 1)=cgetg(1, t_VEC);//Setting U=[[]], so that the first time normalizedbasis is called, it works

  long fourn=4*nfdeg;//The lg of a normal entry
  GEN nformpart=qalg_normform(Q);
  GEN nfdecomp, nform;//nfdecomp used if type=1, nform if type=2
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

  mid=avma;
  long pass=0, ooend=0, nsidesp1=1;//pass tracks which pass we are on
  GEN oosides=cgetg(1, t_VEC), ang1, ang2;//Initialize, so that lg(oosides)=1 to start.
  long neltsgen=0;

  for(;;){
    pass++;
    if(nsidesp1>1){//We have a partial domain.
      oosides=normalizedboundary_oosides(U);
      ooend=lg(oosides)-1;
    }
    iN=itos(gfloor(N))+1;
    points=cgetg(iN, t_VEC);
    for(long i=1;i<iN;i++){//Random points in ball of radius R
      if(i<=ooend){//Going near infinite side
        ang2=gmael(U, 4, oosides[i]);
        if(oosides[i]==1) ang1=gmael(U, 4, lg(gel(U, 1))-1);//Last side, which is the previous side
        else ang1=gmael(U, 4, oosides[i]-1);
        w=randompoint_udarc(R, ang1, ang2, prec);
      }
      else w=randompoint_ud(R, prec);//Random point
      GEN smallelts;
      if(type==1) smallelts=qalg_smallnorm1elts_qfminim(Q, p, C, gen_0, w, maxelts, nfdecomp, nformpart, prec);
      else smallelts=qalg_smallnorm1elts_condition(Q, p, C, gen_0, w, maxelts, nform, nformpart, prec);
      if(smallelts) gel(points, i)=smallelts;
      else return gc_NULL(top);//Must increase precision I guess.
    }
    points=shallowconcat1(points);
    GEN Unew=normalizedbasis(points, U, mats, id, Q, &qalg_fdomm2rembed, &qalg_fdommul, &qalg_fdominv, &qalg_istriv, tol, prec);
    if(gcmp(gel(Unew, 6), areabound)<0){
      long minnew=1, maxnew=lg(points)-1;
      bot=avma;
      while(minnew+1<maxnew){
        set_avma(bot);
        long cind=(maxnew+minnew)/2;
        GEN newpoints=cgetg(cind+1,t_VEC);
        for(long i=1;i<=cind;i++) gel(newpoints, i)=gel(points, i);
        Unew=normalizedbasis(newpoints, U, mats, id, Q, &qalg_fdomm2rembed, &qalg_fdommul, &qalg_fdominv, &qalg_istriv, tol, prec);
        if(gcmp(gel(Unew, 6), areabound)<0) maxnew=cind;
        else minnew=cind;
      }
      GEN ret=mkvec3(stoi(neltsgen+maxnew), stoi(lg(gel(Unew, 1))-1), area);
      return gerepilecopy(top, ret);
    }
    U=Unew;
    neltsgen=neltsgen+lg(points)-1;
    if(pass>1 && (ooend==0 || nsidesp1==lg(gel(U, 1)))) R=gadd(R, epsilon);//Updating R_n
    nsidesp1=lg(gel(U, 1));//How many sides+1
    if(gc_needed(top, 2)) gerepileall(mid, 3, &U, &N, &R);
  }
}




//3: REGRESSIONS & PLOTS


/*Perform ordinary least squares regression. X is a matrix whose columns are the parameters, and y is a column vector of results. Must include linear term as first variable of X.
The formula is B=Bhat=(X*X^T)^(-1)Xy, for the ordinary least squares regression for y=X^T*B+error (formula differs to Wikipedia due to X here being the transpose of what they define there.
Returns [best fit, R^2]*/
GEN
OLS(GEN X, GEN y, int retrsqr)
{
  pari_sp top=avma;
  if(typ(y)!=t_COL) y=gtocol(y);
  if(lg(y)!=lg(X)) pari_err_TYPE("The inputs must have the same length.", mkvec2(X, y));
  GEN Xy=RgM_RgC_mul(X, y);
  GEN XTX=RgM_multosym(X, shallowtrans(X));//X*X^T, which is symmetric
  GEN fit=RgM_solve(XTX, Xy);//Best fit.
  if(!fit) pari_err(e_MISC, "Could not compute matrix inverse. Error with given data or with precision?");
  if(!retrsqr) return gerepileupto(top, fit);
  GEN rsqr=rsquared(X, y, fit);
  return gerepilecopy(top, mkvec2(fit, rsqr));
}

//Performs OLS where we have one independant variable and assume the intercept is 0 (so y=ax). The formua is now sum(x*y)/sum(x^2).
GEN
OLS_nointercept(GEN X, GEN y, int retrsqr)
{
  pari_sp top=avma;
  GEN xysum=gen_0, xsqrsum=gen_0;
  if(lg(X)!=lg(y)) pari_err_TYPE("The inputs must have the same length.", mkvec2(X, y));
  for(long i=1;i<lg(X);i++){
    xysum=gadd(xysum, gmul(gel(X, i), gel(y, i)));
    xsqrsum=gadd(xsqrsum, gsqr(gel(X, i)));
  }
  GEN fit=gdiv(xysum, xsqrsum);
  if(!retrsqr) return gerepileupto(top, fit);
  long lX=lg(X);
  GEN M=cgetg(lX, t_MAT);
  for(long i=1;i<lX;i++) gel(M, i)=mkcol2(gen_1, gel(X, i));
  GEN rsqr=rsquared(M, y, mkcol2(gen_0, fit));
  return gerepilecopy(top, mkvec2(fit, rsqr));
}

//OLS, where there is only one input variable. This just puts it into a matrix form and calls OLS, and is included for convenience.
GEN
OLS_single(GEN x, GEN y, int retrsqr)
{
  pari_sp top=avma;
  long lgx=lg(x);
  GEN xmat=cgetg(lgx, t_MAT);
  for(long i=1;i<lgx;i++) gel(xmat, i)=mkcol2(gen_1, gel(x, i));
  return gerepileupto(top, OLS(xmat, y, retrsqr));
}

//Given inputs for OLS and the proposed linear fit, this returns the R^2 value of the regression.
GEN
rsquared(GEN X, GEN y, GEN fit)
{
  pari_sp top=avma;
  long n=lg(y)-1;//Number of observations
  GEN predicted=RgV_RgM_mul(shallowtrans(fit), X);//1xn matrix of the fitted values.
  GEN yavg=gen_0;
  for(long i=1;i<=n;i++) yavg=gadd(yavg, gel(y, i));
  yavg=gdivgs(yavg, n);//Average value of y
  GEN sstot=gen_0;
  GEN ssres=gen_0;
  for(long i=1;i<=n;i++){
    sstot=gadd(sstot, gsqr(gsub(gel(y, i), yavg)));
    ssres=gadd(ssres, gsqr(gsub(gel(y, i), gel(predicted, i))));
  }
  return gerepileupto(top, gsubsg(1, gdiv(ssres, sstot)));
}

//Assumes plots/build/fname_plotter.tex points to a file making a plot. This compiles the plot, and views it if WSL=1.
//If there isn't enough memory, something like lualatex -shell-escape NAME_plotter.tex works.
void
plot_compile(char *fname, int WSL)
{
  pari_sp top=avma;
  char *line=pari_sprintf("(cd ./plots/build && pdflatex --interaction=batchmode -shell-escape %s_plotter.tex)", fname);//Build
  int s=system(line);
  if(s==-1) pari_err(e_MISC, "ERROR EXECUTING COMMAND");
  pari_free(line);
  line=pari_sprintf("mv -f ./plots/build/%s.pdf ./plots/", fname);//Move the file
  s=system(line);
  if(s==-1) pari_err(e_MISC, "ERROR EXECUTING COMMAND");
  pari_free(line);
  if(WSL){
    line=pari_sprintf("cmd.exe /C start plots/%s.pdf", fname);//Open the file
    s=system(line);
    if(s==-1) pari_err(e_MISC, "ERROR EXECUTING COMMAND");
    pari_free(line);
  }
  set_avma(top);
}

