#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>
#include <fftw3.h>
#include <assert.h>
#include <complex.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_eigen.h>

#include "gslCompute.h"
#include "MVNsampling.h"
#include "Input.h"
#include "Basis.h"
#include "SampSetUp.h"
#include "Output_t0.h"
#include "SampInit.h"
#include "SampZinter.h"
#include "SampTheta.h"
#include "SampPhi2Inv.h"
#include "SampLoading.h"
#include "SampLatent.h"
#include "SampMu.h"
#include "SampAlpha.h"
#include "SampBeta.h"
#include "SampGamma.h"
#include "Sampweights.h"
#include "SampErr_e.h"
#include "SampErr_zeta.h"
#include "SampErr_eps.h"
#include "testing_training.h"
#include "SampMean.h"
#include "SampOutput.h"

#include <R.h>
#include <Rinternals.h>
//#include <Rmath.h>

SEXP csblf(SEXP outputpath, SEXP seed, SEXP nburnin, SEXP niter){
    gsl_set_error_handler_off();
    //////////////////////////////////////////////////////////
    ///// GSL Random Number Generator Initialization /////
    //////////////////////////////////////////////////////////
    const gsl_rng_type * T;
    gsl_rng * r;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);
    gsl_rng_set(r, asInteger(seed));
    //////////////////////////////////////////////////////////
    ///// Input Data /////
    //////////////////////////////////////////////////////////
    // Input/Output data directory
    char *inpathx, *inpath, *outpath; //, *outputpath;
    inpathx = (char *)calloc(500,sizeof(char));
    inpath = (char *)calloc(500,sizeof(char));
    outpath = (char *)calloc(500,sizeof(char));
    char *outppath = CHAR(asChar(outputpath));
    
    char dataSource[20] = "RealData"; // data from simulation or read data source
    int sim = strcmp(dataSource, "Simulation") == 0 ? 1 : 0;
    strcat(inpathx, outppath);
    strcat(inpathx, "Data/");
    strcat(inpath, inpathx);
    strcat(outpath, outppath);
    strcat(outpath, "Result/");

    // Input data
    struct Inputdata data;
    data = input(inpathx, inpath, sim);
    int L = data.sizes[0]; // number of parcels; restricted to be 1
    int nobs = data.sizes[1]; // number of training observations
    int P = data.sizes[3]; // number of imaging predictors
    int nts = data.sizes[4]; // number of observations for test
    
    ///// Basis Function /////
    // Define basis functions
    struct BasisFunc BF;
    double bandwidth = 1.0/10.0;
    int dd = 6.0;
    BF = genBasis(L, outpath, data, bandwidth, dd);
    int M = BF.M;
    //////////////////////////////////////////////////////////
    ///// MCMC /////
    //////////////////////////////////////////////////////////
    // Number of latent factors
    int K ;
    K = strcmp(dataSource, "Simulation") == 0 ? 20 : 9;
    // Length of chain
    int iter = asInteger(niter); // total iterations
    int burnin = asInteger(nburnin); // iterations after burnin
    // Parameters
    struct Sampling PostSamp;
    PostSamp = setupSamp(M, nobs, nts, L, P, K, BF, data);
    // Initialization
    bool printInit = false;
    //set_initial2(L, nobs, nts, K, P, PostSamp, data, BF, r, outpath, printInit, inpath);
    set_initial(L, nobs, nts, K, P, PostSamp, data, BF, r, outpath, printInit);
    clock_t start, end;
    start = clock();
    int *singular = (int *)calloc(2, sizeof(int));
    int t;
    // Start MCMC
    Rprintf("************ Start MCMC ************\n");
    for(t=1; t<=iter; t++){
        if (t % 50 == 0) {
            Rprintf("**** t=%d *****\n", t);
        }
        // Post sampling of Zinter
        Zinter_samp(nobs, L, K, PostSamp, data, BF, r, singular);
        err2inv_u_samp(L, PostSamp, data, r);
        err2inv_e_samp(nobs, L, PostSamp, data, BF, r);
        // Post sampling of theta
        theta_samp(nobs, L, K, PostSamp, data, BF, r, singular);
        // Post sampling of Phi2Inv
        phi2inv_samp(K, nobs, L, PostSamp, data, BF, r);
        // Post sampling of Loadings and its hyperparameters
        load_samp(K, L, nobs, PostSamp, BF, r);
        // Update latent variables
        latent_samp(L, K, nobs, PostSamp, data, BF, r);
        // Update mu
        mu_samp(L, K, nobs, PostSamp, data, BF, r);
        // Update alpha
        alpha_samp(nobs, K, P, L, PostSamp, data, BF, r);
        // Update beta
        beta_samp_approx(nobs, L, K, PostSamp, data, BF, r);
        // Update Gamma
        gamma_samp(P, nobs, K, L, PostSamp, data, BF, r);
        // Update Error
        err2inv_zeta_samp(nobs, L, K, PostSamp, data, BF, r);
        err2inv_eps_samp(nobs, L, K, PostSamp, data, BF, r);
        ///////// Posterior estimations and predictions /////////
        if(t>(iter-burnin) && t%1==0){
            // Training set
            est_training(nobs, P, L, K, PostSamp, data, BF);
            // test
            if(nts>0){
                est_testing(nts, nobs, P, L, K, PostSamp, data, BF);
            }
            // For posterior mean
            samp_mean(L, K, nobs, nts, P, PostSamp, data, BF);
        }
        // Write posterior samplings
        if(t%50==0){
            output_samp(outpath, t, burnin, nobs, nts, M, K, L, P, singular, PostSamp, data, BF);
        }
        // Write posterior mean
        if(t==iter){
            output_mean(L, K, nobs, nts, P, outpath, PostSamp, data, BF, burnin);
        }
    }
    end = clock();
    double mcmc_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    Rprintf("MCMC: %.2f secs\n", mcmc_time_used);
    ////////////////////////////////////////////////////////////
    ///// Release Memory /////
    ////////////////////////////////////////////////////////////
    Rprintf("Release Memory.\n");
    gsl_rng_free(r);
    free(outpath);
    free(inpath);
    free(inpathx);
    free(singular);
    // Input Data
    free(data.sizes);
    free(data.parcel_len);
    free(data.parcel_len_sum);
    free(data.axes);
    free(data.Z);
    free(data.X);
    int l;
    for(l=0; l<L; l++){
        free(data.Zl[l]);
        free(data.Xl[l]);
    }
    if(nts > 0){
        free(data.Z_test);
        free(data.X_test);
        for(l=0; l<L; l++){
            free(data.Zl_test[l]);
            free(data.Xl_test[l]);
            free(PostSamp.xb_test[l]);
        }
        free(data.Zl_test);
        free(data.Xl_test);
        free(PostSamp.xb_test);
    }
    // Basis Functions
    for(l=0; l<L; l++){
        free(BF.basis[l]);
        free(BF.kernel_loc[l]);
    }
    free(BF.basis);
    free(BF.kernel_loc);
    free(BF.Ml);
    // Post Sampling
    for(l=0; l<L; l++){
        free(PostSamp.basis2[l]);
        free(PostSamp.zb[l]);
        free(PostSamp.zb2[l]);
        free(PostSamp.xb[l]);
        free(PostSamp.rxb[l]);
        
        free(PostSamp.Zinter[l]);
        free(PostSamp.theta[l]);
        free(PostSamp.load[l]);
        free(PostSamp.loadstar[l]);
        free(PostSamp.latent[l]);
        free(PostSamp.latentstar[l]);
        free(PostSamp.latent2[l]);
        
        free(PostSamp.beta[l]);
        free(PostSamp.betastar[l]);
        free(PostSamp.alpha[l]);
        free(PostSamp.alphastar[l]);
        free(PostSamp.xba[l]);
        free(PostSamp.xbar[l]);
        free(PostSamp.mustar[l]);
        free(PostSamp.sumX[l]);
    }
    
    free(PostSamp.basis2);
    free(PostSamp.zb);
    free(PostSamp.zb2);
    free(PostSamp.xb);
    free(PostSamp.rxb);
    free(PostSamp.Zinter);
    free(PostSamp.theta);
    free(PostSamp.load);
    free(PostSamp.loadstar);
    free(PostSamp.latent);
    free(PostSamp.latentstar);
    free(PostSamp.latent2);
    free(PostSamp.Phi2Inv);
    free(PostSamp.beta);
    free(PostSamp.betastar);
    free(PostSamp.alpha);
    free(PostSamp.alphastar);
    free(PostSamp.xba);
    free(PostSamp.xbar);
    free(PostSamp.mustar);
    free(PostSamp.sumX);
    free(PostSamp.err2inv_e);
    free(PostSamp.err2inv_zeta);
    free(PostSamp.err2inv_eps);
    free(PostSamp.err2inv_u);
    
    for(l=0; l<L; l++){
        free(PostSamp.mean_Zinter[l]);
        free(PostSamp.mean_theta[l]);
        free(PostSamp.mean_load[l]);
        free(PostSamp.mean_loadstar[l]);
        free(PostSamp.mean_latent[l]);
        free(PostSamp.mean_latentstar[l]);
        free(PostSamp.mean_beta[l]);
        free(PostSamp.mean_betastar[l]);
        free(PostSamp.mean_alpha[l]);
        free(PostSamp.mean_alphastar[l]);
        free(PostSamp.mean_mu[l]);
    }
    free(PostSamp.mean_Zinter);
    free(PostSamp.mean_theta);
    free(PostSamp.mean_load);
    free(PostSamp.mean_loadstar);
    free(PostSamp.mean_latent);
    free(PostSamp.mean_latentstar);
    free(PostSamp.mean_Phi2Inv);
    free(PostSamp.mean_beta);
    free(PostSamp.mean_betastar);
    free(PostSamp.mean_alpha);
    free(PostSamp.mean_alphastar);
    free(PostSamp.mean_mu);
    free(PostSamp.mean_gamma);
    free(PostSamp.mean_gamma_update);
    free(PostSamp.mean_w);
    free(PostSamp.err2inv_e_mean);
    free(PostSamp.err2inv_zeta_mean);
    free(PostSamp.err2inv_eps_mean);
    free(PostSamp.err2inv_u_mean);
    
    if(nts > 0){
        for(l=0; l<L; l++){
            free(PostSamp.Zinter_test[l]);
            free(PostSamp.latent_test[l]);
            free(PostSamp.theta_test[l]);
            free(PostSamp.outcome_test[l]);
            free(PostSamp.Zinter_test_mean[l]);
            free(PostSamp.latent_test_mean[l]);
            free(PostSamp.theta_test_mean[l]);
            free(PostSamp.outcome_test_mean[l]);
            free(PostSamp.outcome_test_mean2[l]);
        }
        free(PostSamp.Zinter_test);
        free(PostSamp.latent_test);
        free(PostSamp.theta_test);
        free(PostSamp.outcome_test);
        free(PostSamp.Zinter_test_mean);
        free(PostSamp.latent_test_mean);
        free(PostSamp.theta_test_mean);
        free(PostSamp.outcome_test_mean);
        free(PostSamp.outcome_test_mean2);
    }
    
    for(l=0; l<L; l++){
        free(PostSamp.latent_train[l]);
        free(PostSamp.latent_train_err[l]);
        free(PostSamp.theta_train[l]);
        free(PostSamp.theta_train_err[l]);
        free(PostSamp.outcome_train[l]);
        free(PostSamp.latent_train_mean[l]);
        free(PostSamp.theta_train_mean[l]);
        free(PostSamp.outcome_train_mean[l]);
        free(PostSamp.outcome_train_mean2[l]);
    }
    free(PostSamp.latent_train);
    free(PostSamp.latent_train_err);
    free(PostSamp.theta_train);
    free(PostSamp.theta_train_err);
    free(PostSamp.outcome_train);
    free(PostSamp.latent_train_mean);
    free(PostSamp.theta_train_mean);
    free(PostSamp.outcome_train_mean);
    free(PostSamp.outcome_train_mean2);
    
    Rprintf("MCMC complete.\n");
    return R_NilValue;
}
