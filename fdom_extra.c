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
static GEN algfromnormdisc(GEN F, GEN D, GEN infram);


//SECTION 1: GEOMETRY METHODS


//Writes the set of arcs in the given colour to filename, so that it is ready to be executed by python.
void python_printarcs(GEN arcs, char *filename, int view, char *extrainput, long prec){
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
  avma=top;
}

//Launches the plotviewer with the given inputs.
void python_plotviewer(char *input){
  char *command;
  command=pari_sprintf("cmd.exe /C start py fdviewer.py %s", input);
  int s=system(command);
  if(s==-1) pari_err(e_MISC, "ERROR EXECUTING COMMAND");
  pari_free(command);
}

//Writes the fundamental domain corresponding to U
void python_printfdom(GEN U, char *filename, long prec){
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
	if(i==1) v1=gel(gel(U, 3), lg(arcs)-1);
	else v1=gel(gel(U, 3), i-1);
	v2=gel(gel(U, 3), i);//The two vertices
    pari_fprintf(f, "%lf %lf %lf %lf %lf\n", rtodbl(gtofp(real_i(gel(arc, 1)), prec)), rtodbl(gtofp(imag_i(gel(arc, 1)), prec)), rtodbl(gtofp(gel(arc, 2), prec)), rtodbl(gmul(garg(gsub(v1, gel(gel(gel(U, 2), i), 1)), prec), fact)), rtodbl(gmul(garg(gsub(v2, gel(gel(gel(U, 2), i), 1)), prec), fact)));
  }
  fclose(f);
  avma=top;
}



//SECTION 2: QUATERNIONIC METHODS


//Returns a quaternion algebra over F (of degree n) with |N_{F/Q}(discriminant)|=D and infinite ramification prescribed by infram (a length n vector of 0's/1's), if it exists. If it does not, this returns 0. This is not gerepile suitable, it leaves a dirty stack.
static GEN algfromnormdisc(GEN F, GEN D, GEN infram){
  pari_sp top=avma;
  if(typ(D)!=t_INT || signe(D)!=1) pari_err_TYPE("D should be a positive integer", D);//The absolute norm to Q should be a positive integer
  GEN pfac=Z_factor(D);
  long nfacsp1=lg(gel(pfac, 1));//# prime factors+1
  long nramplaces=nfacsp1-1;
  for(long i=1;i<lg(infram);i++) if(!gequal0(gel(infram, i))) nramplaces++;//Adding the number of oo ramified places
  if(nramplaces%2==1){avma=top;return gen_0;}//Odd number of ramification places, BAD
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
	if(gequal0(gel(pfacideals, i))){avma=top;return gen_0;}//Nope, return 0
	gel(hass, i)=gen_1;//The vector of 1's
  }
  return alginit(F, mkvec3(gen_2, mkvec2(pfacideals, hass), infram), -1, 1);
}

//Returns G[L[1]]*G[L[2]]*...*G[L[n]], where L is a vecsmall or vec
GEN algmulvec(GEN A, GEN G, GEN L){
  pari_sp top=avma;
  GEN Lsmall=gtovecsmall(L);//L in vecsmall
  GEN elt=gel(alg_get_basis(A), 1);//The identity
  pari_CATCH(CATCH_ALL){
	avma=top;
    pari_CATCH_reset();
	pari_err_TYPE("Invalid inputs; perhaps G is not formatted correctly or L has indices that are too large?", mkvec2(G, L));
	return gen_0;
  }
  pari_TRY{
    for(long i=1;i<lg(Lsmall);i++) elt=algmul(A, elt, gel(G, Lsmall[i]));
  }
  pari_ENDCATCH
  return gerepileupto(top, elt);
}

//Returns a quaternion algebra over F (of degree n) with |N_{F/Q}(discriminant)|=D and split at the infinite place place only, if this exists. We also guarantee that a>0. F must be a totally real field.
GEN algshimura(GEN F, GEN D, long place, long maxcomptime, int allowswap){
  pari_sp top=avma;
  if(nf_get_r2(F)>0) return gen_0;//Not totally real!
  long n=nf_get_degree(F);
  if(place<=0 || place>n) return gen_0;//Invalid place!
  GEN infram=cgetg(n+1, t_VEC);
  for(long i=1;i<=n;i++) gel(infram, i)=gen_1;
  gel(infram, place)=gen_0;
  GEN A=algfromnormdisc(F, D, infram);
  if(gequal0(A)){avma=top;return gen_0;}//Nope
  GEN L=alg_get_splittingfield(A), pol=rnf_get_pol(L);//L=K(sqrt(a))
  long varn=rnf_get_varn(L);
  GEN a=gneg(gsubst(pol, varn, gen_0));//Polynomial is of the form x^2-a, so plug in 0 and negate to get a
  
  if(maxcomptime) pari_alarm(maxcomptime);
  pari_CATCH(CATCH_ALL){
	pari_warn(warner, "Time limit or memory exceeded, skipping.");
	avma=top;
	pari_CATCH_reset();
	return gen_0;
  }
  pari_TRY{
    if(gsigne(gsubst(a, nf_get_varn(F), gel(nf_get_roots(F), place)))!=1){//Must swap a, b
	  GEN b=lift(alg_get_b(A));
	  GEN aden=Q_denom(a), bden=Q_denom(b);
      if(!isint1(aden)) a=gmul(a, gsqr(aden));//We need to get rid of the denominator of a
	  if(!isint1(bden)) b=gmul(b, gsqr(bden));//We need to get rid of the denominator of b
	  if(allowswap) A=gerepileupto(top, alginit(F, mkvec2(b, a), -1, 1));//Swapping a, b
	  else{avma=top;A=gen_0;}//We do not allow the swap.
    }
	else A=gerepilecopy(top, A);//No swap required.
  }pari_ENDCATCH
  pari_alarm(0);//Stop the alarm!
  return A;
}

//Returns a quaternion algebra over F (of degree n) with |N_{F/Q}(discriminant)|=D and split at the infinite place place only, if this exists. We also guarantee that a>0. F must be a totally real field. This returns the pair [a, b] that make the algebra. If allowswap=0, we do NOT allow the swapping of the a, b output by alginit, and either return the algebra (if it is valid) or 0 otherwise. Naively swapping a, b when deg(F) is large can be far too slow.
GEN algshimura_ab(GEN F, GEN D, long place, int allowswap){
  pari_sp top=avma;
  if(nf_get_r2(F)>0) return gen_0;//Not totally real!
  long n=nf_get_degree(F);
  if(place<=0 || place>n) return gen_0;//Invalid place!
  GEN infram=cgetg(n+1, t_VEC);
  for(long i=1;i<=n;i++) gel(infram, i)=gen_1;
  gel(infram, place)=gen_0;
  GEN A=algfromnormdisc(F, D, infram);
  if(gequal0(A)){avma=top;return gen_0;}//Nope
  GEN L=alg_get_splittingfield(A), pol=rnf_get_pol(L);//L=K(sqrt(a))
  long varn=rnf_get_varn(L);
  GEN a=gneg(gsubst(pol, varn, gen_0));//Polynomial is of the form x^2-a, so plug in 0 and negate to get a
  GEN b=lift(alg_get_b(A));
  GEN aden=Q_denom(a), bden=Q_denom(b);//When initializing the algebra, PARI insists that a, b have no denominators
  if(!isint1(aden)) a=gmul(a, gsqr(aden));//We need to get rid of the denominator of a
  if(!isint1(bden)) b=gmul(b, gsqr(bden));//We need to get rid of the denominator of b
  if(gsigne(gsubst(a, nf_get_varn(F), gel(nf_get_roots(F), place)))!=1){
    if(allowswap) return gerepilecopy(top, mkvec2(b, a));//Must swap a, b
	avma=top;//If we do not allow the swapping of a, b (recommended if deg(F)>=6), we do not return anything.
	return gen_0;
  }
  return gerepilecopy(top, mkvec2(a, b));//No swap required.
}

//Returns the algebra where a, b are swapped
GEN algswapab(GEN A){
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

//Returns nwant pairs [a, b] for quaternion algebras over the field F that are suitable for fundamental domains. We start at N_F/Q(disc)=Dmin if supplied, and 2 otherwise.
GEN smallalgebras(GEN F, long nwant, GEN Dmin, int allowswap){
  pari_sp top=avma;
  GEN v=cgetg(nwant+1, t_VEC);
  long ind=1;
  GEN D=gceil(gmax(Dmin, gen_2));//Starting D
  while(ind<=nwant){
	GEN A=algshimura_ab(F, D, 1, allowswap);
	if(!gequal0(A)){
	  gel(v, ind)=A;
	  ind++;
	}
	D=gaddgs(D, 1);
  }
  return gerepilecopy(top, v);
}


/*gensmallalgebras(F, nwant)={
  my(Dset, found, A);
  Dset=vector(nwant);
  found=0;
  for(D=2,oo,
    A=algshimura(F, D);
	if(A!=0,
	  found++;
	  Dset[found]=D;
	  if(found==nwant, return(Dset));
	);
  );
}*/

