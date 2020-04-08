/*
IgPhyML: a program that computes maximum likelihood phylogenies under
non-reversible codon models designed for antibody lineages.

Copyright (C) Kenneth B Hoehn. Sept 2016 onward.

built upon

codonPHYML: a program that  computes maximum likelihood phylogenies from
CODON homologous sequences.

Copyright (C) Marcelo Serrano Zanetti. Oct 2010 onward.

built upon

PhyML:  a program that  computes maximum likelihood phylogenies from
DNA or AA homologous sequences.

Copyright (C) Stephane Guindon. Oct 2003 onward.

All parts of the source except where indicated are distributed under
the GNU public licence. See http://www.opensource.org for details.

*/

#ifndef STATS_H
#define STATS_H

#include "utilities.h"
#include "free.h"
#include "lk.h"
#include "optimiz.h"
#include "models.h"
#include "eigen.h"


phydbl *Covariance_Matrix(t_tree *tree);
phydbl *Hessian(t_tree *tree);
void   Recurr_Hessian(t_node *a, t_node *b, int plus_minus, phydbl *inc, phydbl *res, int *is_ok, t_tree *tree);
phydbl stdnormal_inv(phydbl p);
phydbl Uni();
int    Rand_Int(int min, int max);
phydbl Ahrensdietergamma(phydbl alpha);
phydbl Rgamma(phydbl shape, phydbl scale);
phydbl Rexp(phydbl lambda);
phydbl Bico(int n, int k);
phydbl Factln(int n);
phydbl Gammln(phydbl xx);
phydbl Pbinom(int N, int ni, phydbl p);
phydbl LnGamma (phydbl alpha);
phydbl IncompleteGamma(phydbl x, phydbl alpha, phydbl ln_gamma_alpha);
phydbl PointChi2 (phydbl prob, phydbl v);
phydbl Bivariate_Normal_Density(phydbl x, phydbl y, phydbl mux, phydbl muy, phydbl sdx, phydbl sdy, phydbl rho);
phydbl PointNormal (phydbl prob);
int    DiscreteGamma (phydbl freqK[], phydbl rK[],phydbl alfa, phydbl beta, int K, int median);
phydbl Pnorm(phydbl x, phydbl mean, phydbl var);
phydbl Dnorm_Moments(phydbl x, phydbl mean, phydbl var);
phydbl Dnorm(phydbl x, phydbl mean, phydbl sd);
phydbl Pgamma(phydbl x, phydbl shape, phydbl scale);
phydbl Dgamma_Moments(phydbl x, phydbl mean, phydbl var);
phydbl Dgamma(phydbl x, phydbl shape, phydbl scale);
phydbl LnFact(int n);
int    Choose(int n, int k);
phydbl Ppois(phydbl x, phydbl param);
phydbl Dexp(phydbl x, phydbl param);
phydbl Dpois(phydbl x, phydbl param);
phydbl Rand_Normal_Deviate(phydbl mean, phydbl sd);
phydbl Rnorm(phydbl mean, phydbl sd);
phydbl *Rnorm_Multid(phydbl *mu, phydbl *cov, int dim);
phydbl Rnorm_Trunc(phydbl mean, phydbl sd, phydbl min, phydbl max, int *err);
phydbl *Rnorm_Multid_Trunc(phydbl *mean, phydbl *cov, phydbl *min, phydbl *max, int dim);
phydbl *Hessian_Log(t_tree *tree);
void Recurr_Hessian_Log(t_node *a, t_node *d, int plus_minus, phydbl *inc, phydbl *res, int *is_ok, t_tree *tree);
phydbl Log_Det(int *is_ok, t_tree *tree);
phydbl Dnorm_Trunc(phydbl x, phydbl mean, phydbl sd, phydbl lo, phydbl up);
phydbl Normal_Trunc_Mean(phydbl mu, phydbl sd, phydbl min, phydbl max);
phydbl Constraint_Normal_Trunc_Mean(phydbl wanted_mu, phydbl sd, phydbl min, phydbl max);
phydbl Dnorm_Multi(phydbl *x, phydbl *mu, phydbl *cov, int size, int _log);
phydbl Dnorm_Multi_Given_InvCov_Det(phydbl *x, phydbl *mu, phydbl *invcov, phydbl det, int size, int _log);
phydbl Prop_Log_Dnorm_Multi_Given_InvCov_Det(phydbl *x, phydbl *mu, phydbl *invcov, phydbl det, int size);
phydbl Log_Dnorm(phydbl x, phydbl mean, phydbl sd, int *err);
phydbl tt800();
phydbl Pnorm_Ihaka_Derived_From_Cody(phydbl x);
int Matinv(phydbl *x, int n, int m);
phydbl *Matrix_Mult(phydbl *A, phydbl *B, int nra, int nca, int nrb, int ncb);
phydbl *Matrix_Transpose(phydbl *A, int dim);
void Normal_Conditional(phydbl *mu, phydbl *cov, phydbl *a, int n, short int *is_1, int n1, phydbl *cond_mu, phydbl *cond_var);
void Normal_Conditional_Unsorted(phydbl *mu, phydbl *cov, phydbl *a, int n, short int *is_1, int n1, phydbl *cond_mu, phydbl *cond_cov);
phydbl Matrix_Det(phydbl *A, int size, int _log);
void Get_Reg_Coeff(phydbl *mu, phydbl *cov, phydbl *a, int n, short int *is_1, int n1, phydbl *reg_coeff);
phydbl Rnorm_Trunc_Inverse(phydbl mean, phydbl sd, phydbl min, phydbl max, int *error);
phydbl Norm_Trunc_Sd(phydbl mu, phydbl sd, phydbl a, phydbl b);
phydbl Norm_Trunc_Mean(phydbl mu, phydbl sd, phydbl a, phydbl b);
int Norm_Trunc_Mean_Sd(phydbl mu, phydbl sd, phydbl a, phydbl b, phydbl *trunc_mu, phydbl *trunc_sd);
phydbl Log_Dnorm_Trunc(phydbl x, phydbl mean, phydbl sd, phydbl lo, phydbl up, int *err);
phydbl  Pnorm_Marsaglia(phydbl x);


#endif
