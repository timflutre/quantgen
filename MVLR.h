/** \file MVLR.h
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

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <vector>


class MVLR {
  
 private:
  
  int n; // sample size
  int s; // subgroup size
  int q; // control size
  int p; // genotype size

  int m; // wishart prior default m = s-1

  int es_option;  // if scale the effect size by noise level

  int sigma_option;


  gsl_matrix *Y; // phenotype matrix nxs
  gsl_matrix *Xg; // genotype matrix nxp
  gsl_matrix *Xc; // control matrix nxq
  gsl_matrix *H; // wishart prior sxs 


  gsl_matrix *XctXc;
  gsl_matrix *XctXc_inv;

  gsl_matrix *Pg;  // effect prior constructor s x s
  gsl_matrix *Wg; // effect prior ps x ps 

  gsl_matrix *Sigma; // residual covariance estimate sxs
  gsl_matrix *Sigma_inv; // Sigma^{-1}

  gsl_matrix *Sigma0; // residual covariance estimate sxs
  gsl_matrix *Sigma0_inv; // Sigma^{-1}

  double log_det_Sigma;
  double log_det_Sigma0;

  gsl_matrix *bVi; // (\hat bv)'(Vg^{-1}) 1 x ps
  gsl_matrix *Vg_inv; //  Var(\hat bv) ps x ps
  gsl_matrix *tBc; // \tilde Bc, pxs
  gsl_matrix *hB;  // \hat B
  gsl_matrix *Kg; // part of  Vg^{-1}  pxp 


  vector<vector<int> > skeleton;

 private:
  

  void compute_Sigma1();
  void invert_Sigma1();
  void compute_Sigma2();
  
  void compute_common_core();
  void compute_bVi();


  gsl_matrix *multiple_regression(gsl_matrix * Ys, int sindex);
  
  double log10_weighted_sum(vector<double> &vec, vector<double> &wts);
  gsl_matrix *vec(gsl_matrix *M, int a, int b);
  gsl_matrix *kron (gsl_matrix *M, gsl_matrix *L, int a, int b);
  gsl_matrix *kron2 (gsl_matrix *M, int mr, int mc, gsl_matrix *L, int lr, int lc);
  void print_matrix(gsl_matrix *M, int a, int b);
  
  
 public:
  
  MVLR(vector<vector<double> > & Y_in, vector<vector<double> > & Xg_in, vector<vector<double> > & Xc_in,int es);
  ~MVLR();
  
  void set_skeleton(vector<vector<int> > & GammaV);
  void set_IW_prior(gsl_matrix *H_in, int m_in);
  void set_effect_prior(double phi2, double omg2);
  void set_Sigma_option(int option){
    if(bVi!=0)
      gsl_matrix_free(bVi);
    bVi = 0;
    sigma_option = option;
  };

  double compute_log10_ABF();
  double compute_log10_ABF(double phi2, double omg2, vector<vector<int> > & GammaV);
  double compute_log10_ABF(vector<double> &phi2_vec, vector<double> &omg2_vec);
  

};
