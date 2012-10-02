/** \file MVLR.cc
 *
 *  `MVLR' is a class implementing the multivariate linear regression
 *  Copyright (C) 2012 Xioaquan (William) Wen
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

using namespace std;

#include "MVLR.h"
#include <math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_cdf.h>
#include <map>

MVLR::MVLR(vector<vector<double> > & Y_in, vector<vector<double> > & Xg_in, vector<vector<double> > & Xc_in, int es ){
  
  s = Y_in.size();
  n = Y_in[0].size();
  p = Xg_in.size();
  q = Xc_in.size();
  
  es_option = es;
  

  // for the intercept
  q++;
  
  //printf("n = %d s = %d p = %d q = %d\n",n,s,p,q);

  Y  = gsl_matrix_calloc(n,s);
  Xg = gsl_matrix_calloc(n,p);
  Xc = gsl_matrix_calloc(n,q);
  

  for(int i=0;i<s;i++){
    for(int j=0;j<n;j++){
      gsl_matrix_set(Y, j,i,Y_in[i][j]);
    }
  }    
  
 


  for(int i=0;i<p;i++){
    for(int j=0;j<n;j++){
      gsl_matrix_set(Xg,j,i,Xg_in[i][j]);
    }
  }


 


  if(q>1){
    for(int i=1;i<q;i++){
      for(int j=0;j<n;j++){
	gsl_matrix_set(Xc,j,i,Xc_in[i-1][j]);
      }
    }
  }
  
  for(int j=0;j<n;j++){
    gsl_matrix_set(Xc,j,0,1.0);
  }


  // default value for IW prior on Sigma H = diag(1e-8), m = s-1
  m = s-1;
  H = gsl_matrix_calloc(s,s);
  for(int i=0;i<s;i++){
    gsl_matrix_set(H,i,i,1e-8);
  }

  hB = Pg = tBc = Kg = Wg = bVi = Vg_inv = Sigma = Sigma_inv = XctXc = XctXc_inv=0;
  compute_common_core();

}


// call before compute_Sigma
void MVLR::set_IW_prior(gsl_matrix *H_in, int m_in){
  if(H!=0)
    gsl_matrix_free(H);
  H = H_in;
  m = m_in;
}
  


void MVLR::set_effect_prior(double phi2, double omg2){
  

  if(Pg!=0){
    gsl_matrix_free(Pg);
  }
  
  Pg = gsl_matrix_calloc(s,s);
  for(int i=0;i<s;i++){
    for(int j=0;j<s;j++){
      double v = omg2;
      if(i==j)
        v+= phi2;
      gsl_matrix_set(Pg,i,j,v);
    }
  }
  

}

  


void MVLR::compute_Sigma1(){
 
  if(Sigma !=0){
    gsl_matrix_free(Sigma);
    Sigma = 0;
  }
  
  if(Sigma_inv !=0){ 
    gsl_matrix_free(Sigma_inv);
    Sigma_inv = 0;
  }
  
  if(hB !=0){
    gsl_matrix_free(hB);
    hB = 0;
  }

  hB = gsl_matrix_calloc(p+q,s);
  Sigma = gsl_matrix_calloc(s,s);


  if(skeleton.size()==0){
    gsl_matrix *X = gsl_matrix_calloc(n,p+q);
    for(int i=0;i<n;i++){
      for(int j=0;j<q;j++){
	gsl_matrix_set(X,i,j,gsl_matrix_get(Xc,i,j));
      }
    }
    
    for(int i=0;i<n;i++){
      for(int j=0;j<p;j++){
	gsl_matrix_set(X,i,q+j,gsl_matrix_get(Xg,i,j));
      }
    }
    
  
    gsl_matrix *XtX = gsl_matrix_calloc(p+q, p+q);
    gsl_blas_dgemm(CblasTrans,CblasNoTrans,1,X,X,0,XtX);
    
  
    // compute inverse of XtX (generalized inverse version)
    gsl_matrix *V = gsl_matrix_calloc(p+q,p+q);
    gsl_vector *S = gsl_vector_calloc(p+q);
    gsl_vector *work = gsl_vector_calloc(q+p);
    gsl_linalg_SV_decomp (XtX, V, S,work);
    
    gsl_matrix *t1 = gsl_matrix_calloc(p+q,p+q);
    for(int i=0;i<p+q;i++){
      double v = gsl_vector_get(S,i);
      if(v>1e-8){
	gsl_matrix_set(t1,i,i,1/v);
      }
    }
    
    
    
    gsl_matrix *t2 = gsl_matrix_calloc(p+q,p+q);
    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,V,t1,0,t2);
    
    gsl_matrix *XtX_inv = gsl_matrix_calloc(p+q,p+q);
    gsl_blas_dgemm(CblasNoTrans,CblasTrans,1,t2,V,0,XtX_inv);
    
    //print_matrix(XtX_inv,p+q,p+q);
    
    // (X'X)^{-1)X'
    gsl_matrix *t3 = gsl_matrix_calloc(p+q,n);
    gsl_blas_dgemm(CblasNoTrans,CblasTrans,1,XtX_inv,X,0,t3);
  
    
    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,t3,Y,0,hB);
    
    gsl_matrix *t4 = gsl_matrix_calloc(n,n);
    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,X,t3,0,t4);
    
    gsl_matrix *t5 = gsl_matrix_calloc(n,n);
    for(int i=0;i<n;i++){
      for(int j=0;j<n;j++){
	double v = -gsl_matrix_get(t4,i,j);
	if(i==j)
	  v += 1;
	gsl_matrix_set(t5,i,j,v);
      }
    }
    
    
    gsl_matrix *t6 = gsl_matrix_calloc(s,n);
    gsl_blas_dgemm(CblasTrans,CblasNoTrans,1,Y,t5,0,t6);
    
    gsl_matrix *t7 = gsl_matrix_calloc(s,s);
    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,t6,Y,0,t7);
    
    //print_matrix(t7,s,s);
  

  
    
  
    for(int i=0;i<s;i++){
      for(int j=0;j<s;j++){	
	double v = (gsl_matrix_get(H,i,j)+gsl_matrix_get(t7,i,j))/(n+m-q-s-1);  	
	gsl_matrix_set(Sigma,i,j,v);
      }
    }
 
    gsl_matrix_free(X);
    gsl_matrix_free(XtX);
    gsl_matrix_free(V);
    gsl_vector_free(S);
    gsl_vector_free(work);
    gsl_matrix_free(t1);
    gsl_matrix_free(t2);
    gsl_matrix_free(XtX_inv);
    gsl_matrix_free(t3);
    gsl_matrix_free(t4);
    gsl_matrix_free(t5);
    gsl_matrix_free(t6);
    gsl_matrix_free(t7);
  }else{
    
    gsl_matrix *Yr = gsl_matrix_calloc(n,s);

    for(int i=0;i<s;i++){   
      gsl_matrix *Ys = gsl_matrix_calloc(n,1);
      for(int j=0;j<n;j++)
	gsl_matrix_set(Ys, j,0,gsl_matrix_get(Y,j,i));
      gsl_matrix *Rs = multiple_regression(Ys, i);
      for(int j=0;j<n;j++)
        gsl_matrix_set(Yr, j,i,gsl_matrix_get(Rs,j,0));   
      
      gsl_matrix_free(Rs);
      gsl_matrix_free(Ys);

    }
    
    gsl_blas_dgemm(CblasTrans,CblasNoTrans,1,Yr,Yr,0,Sigma);
    gsl_matrix_add(Sigma,H);
    gsl_matrix_scale(Sigma,1.0/(n+m-q-s-1));
    
    gsl_matrix_free(Yr);
  }   
  

  invert_Sigma1();



}

gsl_matrix * MVLR::multiple_regression(gsl_matrix * Ys, int sindex){


  
  int count = 0;
  map<int,int> pmap;

  for(int i=0;i<p;i++){
    if(skeleton[i][sindex] == 1){
      pmap[count] = i;
      count++;
    }
  }
  
  int pe = pmap.size();
  gsl_matrix *X = gsl_matrix_calloc(n,pe+q);
  for(int i=0;i<n;i++){
    for(int j=0;j<q;j++){
      gsl_matrix_set(X,i,j,gsl_matrix_get(Xc,i,j));
    }
  }

  for(int j=0;j<pe;j++){
    for(int i=0;i<n;i++){
      gsl_matrix_set(X,i,q+j,gsl_matrix_get(Xg,i,pmap[j]));
    }
  }
  
  
  
  gsl_matrix *XtX = gsl_matrix_calloc(pe+q, pe+q);
  gsl_blas_dgemm(CblasTrans,CblasNoTrans,1,X,X,0,XtX);


  // compute inverse of XtX (generalized inverse version)                                              
  gsl_matrix *V = gsl_matrix_calloc(pe+q,pe+q);
  gsl_vector *S = gsl_vector_calloc(pe+q);
  gsl_vector *work = gsl_vector_calloc(q+pe);
  gsl_linalg_SV_decomp (XtX, V, S,work);

  gsl_matrix *t1 = gsl_matrix_calloc(pe+q,pe+q);
  for(int i=0;i<pe+q;i++){
    double v = gsl_vector_get(S,i);
    if(v>1e-8){
      gsl_matrix_set(t1,i,i,1/v);
    }
  }



  gsl_matrix *t2 = gsl_matrix_calloc(pe+q,pe+q);
  gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,V,t1,0,t2);

  gsl_matrix *XtX_inv = gsl_matrix_calloc(pe+q,pe+q);
  gsl_blas_dgemm(CblasNoTrans,CblasTrans,1,t2,V,0,XtX_inv);


  

  //print_matrix(XtX_inv,p+q,p+q);                                                                     

  // (X'X)^{-1)X'                                                                                      
  gsl_matrix *t3 = gsl_matrix_calloc(pe+q,n);
  gsl_blas_dgemm(CblasNoTrans,CblasTrans,1,XtX_inv,X,0,t3);

  	
  gsl_matrix *sB = gsl_matrix_calloc(q+pe,1);
  gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,t3,Ys,0,sB);

  
  
  for(int i=0;i<q;i++)
    gsl_matrix_set(hB,i,sindex, gsl_matrix_get(sB,i,0));
    
  

  for(int i=0;i<pe;i++)
    gsl_matrix_set(hB,q+pmap[i],sindex, gsl_matrix_get(sB,q+i,0));
  
    

  gsl_matrix *t4 = gsl_matrix_calloc(n,n);
  gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,X,t3,0,t4);

  
  
  gsl_matrix *t5 = gsl_matrix_calloc(n,n);
  for(int i=0;i<n;i++){
    for(int j=0;j<n;j++){
      double v = -gsl_matrix_get(t4,i,j);
      if(i==j)
	v += 1;
      gsl_matrix_set(t5,i,j,v);
    }
  }

  

  gsl_matrix *res = gsl_matrix_calloc(n,1);
  gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,t5,Ys,0,res);

  
  
  

  //print_matrix(t7,s,s);                                                                              

  gsl_matrix_free(X);
  gsl_matrix_free(XtX);
  gsl_matrix_free(V);
  gsl_vector_free(S);
  gsl_vector_free(work);
  gsl_matrix_free(t1);
  gsl_matrix_free(t2);
  gsl_matrix_free(XtX_inv);
  gsl_matrix_free(t3);
  gsl_matrix_free(t4);
  gsl_matrix_free(t5);
  gsl_matrix_free(sB);
  
  return res;

}





void MVLR::invert_Sigma1(){

  // small sample correction for \hat Sigma and inversion

  
  vector<double> factor_v;
  for(int i=0;i<s;i++){
    
    gsl_matrix *bg = gsl_matrix_calloc(p,1);
    
    for(int j=0;j<p;j++){
      gsl_matrix_set(bg,j,0,gsl_matrix_get(hB,q+j,i));
    }
    
    gsl_matrix *t1 = gsl_matrix_calloc(1,p);
    gsl_blas_dgemm(CblasTrans,CblasNoTrans,1,bg,Kg,0,t1);
    gsl_matrix *t2 = gsl_matrix_calloc(1,1);
    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,t1,bg,0,t2);
    
    // hotelling's T2
    double rst = gsl_matrix_get(t2,0,0)/gsl_matrix_get(Sigma,i,i);
    
    // convert to F
    double v1 = p;
    double v2 = n+m-q-s-1;
    
    double F = (v2-v1+1)*rst/(v1*v2);
    // quantile transform it into a Chisq 
    double qv = gsl_cdf_fdist_Q(F, v1, v2-v1+1);
 
    //printf("qv = %f    %f  %f\n",qv,v1, v2);
    double new_F = gsl_cdf_chisq_Qinv (qv, p)/p;
    
    //printf("F:  %f  %f\n",F, new_F);
    if(F<1e-8)
      factor_v.push_back(1);
    else
      factor_v.push_back(new_F/F);
    
    gsl_matrix_free(bg);
    gsl_matrix_free(t1);
    gsl_matrix_free(t2);

  }

  gsl_matrix *t3 = gsl_matrix_calloc(s,s);
  for(int i=0;i<s;i++){
    for(int j=0;j<s;j++){

      double v = gsl_matrix_get(Sigma,i,j);
      if(v<0)
	v=0;
      v = v/sqrt(factor_v[i]*factor_v[j]);
      //if(i!=j)
      //v=0;
      //double v = gsl_matrix_get(Sigma,i,j);
      gsl_matrix_set(t3,i,j,v);
    }
  }
  gsl_matrix_memcpy(Sigma,t3);
  
  /*
  printf("\n");
  print_matrix(Sigma,s,s);
  printf("\n");
  print_matrix(hB,p+q,s);
  printf("\n");
  print_matrix(Kg,p,p);
  printf("\n");
  */


  int ss;
  gsl_permutation * pp = gsl_permutation_alloc(s);
  gsl_linalg_LU_decomp (t3, pp, &ss);
  Sigma_inv = gsl_matrix_calloc(s,s);
  gsl_linalg_LU_invert (t3, pp, Sigma_inv);
  
  log_det_Sigma = gsl_linalg_LU_lndet(t3);


  gsl_permutation_free(pp);
  gsl_matrix_free(t3);


  
}






void MVLR::compute_common_core(){ 


  XctXc = gsl_matrix_calloc(q,q);
  gsl_blas_dgemm(CblasTrans,CblasNoTrans,1,Xc,Xc,0,XctXc);
  XctXc_inv = gsl_matrix_calloc(q,q);
  
  
  if(q==1)
    gsl_matrix_set(XctXc_inv,0,0,1.0/gsl_matrix_get(XctXc,0,0));
  else{
    gsl_matrix *t = gsl_matrix_calloc(q,q);
    gsl_matrix_memcpy(t,XctXc);
    int ss;
    gsl_permutation * pp = gsl_permutation_alloc(q);
    gsl_linalg_LU_decomp (t, pp, &ss);
    gsl_linalg_LU_invert (t, pp, XctXc_inv);
    gsl_permutation_free(pp);
    gsl_matrix_free(t);
  }

  gsl_matrix *t0 = gsl_matrix_calloc(q,n);
  gsl_blas_dgemm(CblasNoTrans,CblasTrans,1,XctXc_inv,Xc,0,t0);
  
  tBc = gsl_matrix_calloc(q,s);
  gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,t0,Y,0,tBc);
  
  gsl_matrix_free(t0);


  // 2. V_g^{-1}
  
  gsl_matrix *t1 = gsl_matrix_calloc(n,q);
  gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,Xc,XctXc_inv,0,t1);

  // t2 = Xc(Xc'Xc)^{-1}Xc'
  gsl_matrix *t2 = gsl_matrix_calloc(n,n);
  gsl_blas_dgemm(CblasNoTrans,CblasTrans,1,t1,Xc,0,t2);

  
  gsl_matrix *XgtXg = gsl_matrix_calloc(p,p);
  gsl_blas_dgemm(CblasTrans,CblasNoTrans,1,Xg,Xg,0,XgtXg);
  
  gsl_matrix *t3 = gsl_matrix_calloc(p,n);
  gsl_blas_dgemm(CblasTrans,CblasNoTrans,1,Xg,t2,0,t3);

  gsl_matrix *t4 = gsl_matrix_calloc(p,p);
  gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,t3,Xg,0,t4);
  
  Kg = gsl_matrix_calloc(p,p);
  
  for(int i=0;i<p;i++){
    for(int j=0;j<p;j++){
      gsl_matrix_set(Kg,i,j, gsl_matrix_get(XgtXg,i,j)-gsl_matrix_get(t4,i,j));
    }
  }
  
  gsl_matrix_free(t1);
  gsl_matrix_free(t2);
  gsl_matrix_free(t3);
  gsl_matrix_free(t4);
  gsl_matrix_free(XgtXg);


  gsl_matrix *t5 = gsl_matrix_calloc(q,n);
  gsl_blas_dgemm(CblasNoTrans,CblasTrans,1,XctXc_inv,Xc,0,t5);
  
  gsl_matrix *t6 = gsl_matrix_calloc(n,n);
  gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,Xc,t5,0,t6);
  
  gsl_matrix *t7 = gsl_matrix_calloc(n,n);
  for(int i=0;i<n;i++){
    for(int j=0;j<n;j++){
      double v = -gsl_matrix_get(t6,i,j);
      if(i==j)
	v += 1;
      gsl_matrix_set(t7,i,j,v);
    }
  }
  
    
  gsl_matrix *t8 = gsl_matrix_calloc(s,n);
  gsl_blas_dgemm(CblasTrans,CblasNoTrans,1,Y,t7,0,t8);
  
  gsl_matrix *t9 = gsl_matrix_calloc(s,s);
  gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,t8,Y,0,t9);
  
  
  
  Sigma0 = gsl_matrix_calloc(s,s);
  
  for(int i=0;i<s;i++){
    for(int j=0;j<s;j++){

      double v = (gsl_matrix_get(H,i,j)+gsl_matrix_get(t9,i,j))/(n+m-q-s-1);      
      gsl_matrix_set(Sigma0,i,j,v);
      
    }
  }
  
  

  gsl_matrix *t10 = gsl_matrix_calloc(s,s);
  gsl_matrix_memcpy(t10,Sigma0);
  
  Sigma0_inv = gsl_matrix_calloc(s,s);

  int ss;
  gsl_permutation * pp = gsl_permutation_alloc(s);
  gsl_linalg_LU_decomp (t10, pp, &ss);
  gsl_linalg_LU_invert (t10, pp, Sigma0_inv);

  log_det_Sigma0 = gsl_linalg_LU_lndet(t10);
 
  
  gsl_permutation_free(pp);
  gsl_matrix_free(t5);
  gsl_matrix_free(t6);
  gsl_matrix_free(t7);
  gsl_matrix_free(t8);
  gsl_matrix_free(t9);
  gsl_matrix_free(t10);
  

  






}

void MVLR::compute_bVi(){

  
  // correct Sigma_inv for small sample size
 
  if(sigma_option==1){
    compute_Sigma1();
  }
  
  if(sigma_option==2 ){
    
    if(Sigma!=0){
      gsl_matrix_free(Sigma);
      gsl_matrix_free(Sigma_inv);
    }
    
    Sigma = gsl_matrix_calloc(s,s);
    Sigma_inv = gsl_matrix_calloc(s,s);
    gsl_matrix_memcpy(Sigma, Sigma0);
    gsl_matrix_memcpy(Sigma_inv, Sigma0_inv);
    
  }

  if(Vg_inv !=0)
    gsl_matrix_free(Vg_inv);
 
  if(bVi != 0)
    gsl_matrix_free(bVi);
  
  Vg_inv = kron(Kg,Sigma_inv, p,s);
  //print_matrix(Vg_inv,p*s,p*s);

  // 1. hbvg'V_g^{-1}
  

  gsl_matrix *t5 = gsl_matrix_calloc(n,s);
  gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,Xc,tBc,0,t5);


  gsl_matrix *t6 = gsl_matrix_calloc(s,n);
  for(int i=0;i<n;i++){
    for(int j=0;j<s; j++){
      double v = gsl_matrix_get(Y,i,j)-gsl_matrix_get(t5,i,j);
      gsl_matrix_set(t6,j,i,v);
    }
  }

  gsl_matrix *t7 = vec(t6,s,n);
  
  gsl_matrix *t8 = kron2(Xg,n,p,Sigma_inv,s,s);
  
  

  
  bVi = gsl_matrix_calloc(1,s*p);
  gsl_blas_dgemm(CblasTrans,CblasNoTrans,1,t7,t8,0,bVi);
  
  
 

  gsl_matrix_free(t5);
  gsl_matrix_free(t6);
  gsl_matrix_free(t7);
  gsl_matrix_free(t8);

}




void MVLR::set_skeleton(vector<vector<int> >& GammaV){
  skeleton = GammaV;
  if(sigma_option==1&&bVi!=0){
    compute_Sigma1();
    compute_bVi();
  }
}


MVLR::~MVLR(){
  
  
  gsl_matrix_free(Y);
  gsl_matrix_free(Xg);
  gsl_matrix_free(Xc);
  gsl_matrix_free(H);

  
  if(hB!=0){
    gsl_matrix_free(hB);
  }
  
  if(Wg!=0){
    gsl_matrix_free(Wg);
  }
  
  if(Pg!=0){
    gsl_matrix_free(Pg);
  }

  if(Sigma!=0){
    gsl_matrix_free(Sigma);
    gsl_matrix_free(Sigma_inv);
  }
  if(Sigma0!=0){
    gsl_matrix_free(Sigma0);
    gsl_matrix_free(Sigma0_inv);
  }
  if(bVi!=0)
    gsl_matrix_free(bVi);
  if(Vg_inv!=0)
    gsl_matrix_free(Vg_inv);  
  
  if(tBc!=0)
    gsl_matrix_free(tBc);
  
  if(Kg!=0)
    gsl_matrix_free(Kg);
  
  if(XctXc !=0){
    gsl_matrix_free(XctXc);
    gsl_matrix_free(XctXc_inv);
  }
  
  
}


double MVLR::compute_log10_ABF(){
   
      
  if(bVi==0)
    compute_bVi();
 
  if(Wg!=0){
    gsl_matrix_free(Wg);
    Wg=0;
  }
  
  // construct Wg from Pg

  if(es_option){
    for(int i=0;i<s;i++){
      for(int j=0;j<s;j++){
	double v1 = sqrt(gsl_matrix_get(Sigma,i,i));
	double v2 = sqrt(gsl_matrix_get(Sigma,j,j));
	double v3 = gsl_matrix_get(Pg,i,j);
	gsl_matrix_set(Pg, i,j,v1*v2*v3);
      }
    }
  }

  gsl_matrix *Ip = gsl_matrix_calloc(p,p);
  for(int i=0;i<p;i++){
    gsl_matrix_set(Ip,i,i,1);
  }
  
  Wg = kron(Ip,Pg,p,s);
  
  gsl_matrix_free(Ip);

  

  // modify Wg according to pre-defined skeleton

  if(skeleton.size()!=0){
    
    gsl_matrix *M = gsl_matrix_calloc(p*s,p*s);
    for(int i=0;i<p;i++){
      for(int j=0;j<s;j++){
	if(skeleton[i][j]==1){
	  gsl_matrix_set(M,i*s+j,i*s+j,1.0);
	}
      }
    }
    

    gsl_matrix *t0 = gsl_matrix_calloc(p*s,p*s);
    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,M,Wg,0,t0);
    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,t0,M,0,Wg);
 
    
    gsl_matrix_free(t0);
    gsl_matrix_free(M);
    
  }

  
  gsl_matrix *t1 = gsl_matrix_calloc(p*s,p*s);
  gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,Vg_inv,Wg,0,t1);
  
    


  for(int i=0;i<p*s;i++){
    gsl_matrix_set(t1,i,i,gsl_matrix_get(t1,i,i)+1);
  }
  
  
  int ss;
  gsl_permutation *pp = gsl_permutation_alloc(p*s);
  gsl_linalg_LU_decomp (t1, pp, &ss);
  gsl_matrix *t2 = gsl_matrix_calloc(p*s,p*s);
  gsl_linalg_LU_invert (t1, pp, t2);
  


  gsl_matrix *t3 = gsl_matrix_calloc(p*s,p*s);
  gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,Wg,t2,0,t3);
  
  
  gsl_matrix *t4 = gsl_matrix_calloc(1,p*s);
  gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,bVi,t3,0,t4);


  gsl_matrix *t5 = gsl_matrix_calloc(1,1);
  gsl_blas_dgemm(CblasNoTrans,CblasTrans,1,t4,bVi,0,t5);
  
  double rst = .5*gsl_matrix_get(t5,0,0);
 

  double log_detVal = gsl_linalg_LU_lndet(t1);  
  rst += -0.5*log_detVal;

  gsl_permutation_free(pp);


  
  

  gsl_matrix_free(t1);
  gsl_matrix_free(t2);
  gsl_matrix_free(t3);
  gsl_matrix_free(t4);
  gsl_matrix_free(t5);


  /*
  if(sigma_option==1){


    gsl_matrix *Hbg = gsl_matrix_calloc(p*s,1);
    for(int i=0;i<p;i++){
      for(int j=0;j<s;j++){      
      gsl_matrix_set(Hbg,j+i*s,0,gsl_matrix_get(hB,q+i,j));
      }
    }
    
    
    printf("\n");
    print_matrix(Hbg,p*s,1);
    printf("\n");
    print_matrix(Vg_inv,p*s,p*s);


    gsl_matrix *t6 = gsl_matrix_calloc(1,p*s);
    gsl_blas_dgemm(CblasTrans,CblasNoTrans,1,Hbg,Vg_inv,0,t6);
    
    gsl_matrix *t7 = gsl_matrix_calloc(1,1);
    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,t6,Hbg,0,t7);
    
  

    double rst2 = -.5*gsl_matrix_get(t7,0,0);//+ 0.5*(n+m-q)*(log_det_Sigma0 - log_det_Sigma);
    
    printf("\nrst2 = %f %f  %f\n", gsl_matrix_get(t7,0,0), log_det_Sigma0, log_det_Sigma);
    
    rst += rst2;
    gsl_matrix_free(t6);
    gsl_matrix_free(t7);
    gsl_matrix_free(Hbg);
    
  }
  */

  return rst/log(10.0);

}






double MVLR::compute_log10_ABF(double phi2, double omg2, vector<vector<int> > & GammaV){
  set_effect_prior(phi2,omg2);
  set_skeleton(GammaV);
  return compute_log10_ABF();
}


double MVLR::compute_log10_ABF(vector<double> &phi2_vec, vector<double> &omg2_vec){
  
  int size = phi2_vec.size();
  vector<double> abf_vec;
  vector<double> wts_vec;
  for(int i=0;i<size;i++){
    set_effect_prior(phi2_vec[i],omg2_vec[i]);
    abf_vec.push_back(compute_log10_ABF());
    wts_vec.push_back(1.0/size);
  }

  return log10_weighted_sum(abf_vec, wts_vec);

}








// utilities


double MVLR::log10_weighted_sum(vector<double> &vec, vector<double> &wts){

    double max = vec[0];
    for(int i=0;i<vec.size();i++){
      if(vec[i]>max)
	max = vec[i];
    }
    double sum = 0;
    for(int i=0;i<vec.size();i++){
      sum += wts[i]*pow(10, (vec[i]-max));
    }

    return (max+log10(sum));
}




gsl_matrix *MVLR::vec (gsl_matrix *M, int a, int b){
  
  gsl_matrix *v = gsl_matrix_calloc(a*b, 1);
  int k = 0;
  for(int j=0;j<b;j++){
    for(int i=0;i<a;i++){
      gsl_matrix_set(v,k++,0,gsl_matrix_get(M,i,j));
    }
  }
  
  return v;
}


gsl_matrix *MVLR::kron (gsl_matrix *M, gsl_matrix *L, int a, int b){
  
  gsl_matrix *R = gsl_matrix_calloc (a*b,a*b);
  for (int i = 0; i < a; i++){
    for (int j = 0; j < a; j++){ 
	  for (int k = 0; k < b; k++){
	      for (int l = 0; l < b; l++){
		  gsl_matrix_set (R, i*b+k,j*b+l, gsl_matrix_get (M, i, j)*gsl_matrix_get (L, k, l));
	      }
	  }
    }
  }
  return R;
}

gsl_matrix *MVLR::kron2(gsl_matrix *A, int nrows, int ncols,gsl_matrix *B, int mrows, int mcols){
  
  gsl_matrix *C = gsl_matrix_calloc(nrows*mrows, ncols*mcols);
  
  int cr = 0;
  int cc = 0;
  for (int i = 0; i<nrows;i++){
    for(int j=0;j<ncols;j++){
      cr = i*mrows;
      cc = j*mcols;
      double a = gsl_matrix_get(A,i,j);
      
      for(int k=0;k<mrows;k++){
	for(int l=0;l<mcols;l++){
	  gsl_matrix_set(C,cr+k,cc+l,a*gsl_matrix_get(B,k,l));
	}
      }
      
    }
    
  }
  
  return C;
}




void MVLR::print_matrix(gsl_matrix *M, int a, int b){
  
  for(int i=0;i<a;i++){
    for(int j=0;j<b;j++){
      //printf("%f  ",gsl_matrix_get(M,i,j));
      printf("%e  ",gsl_matrix_get(M,i,j));
    }
    printf("\n");
  }
  printf("\n");
}

