void mu_samp(const int L, const int K, const int nobs,
              struct Sampling PostSamp, const struct Inputdata data, const struct BasisFunc BF,
              gsl_rng *r){
    gsl_set_error_handler_off();
    int l, i, j, k, q, rr, Ml;
    double err;
    
    /*printf("\n#########\nOld_MuStar\n");
    for(i=0; i<3; i++){
        for(j=0; j<K; j++){
            printf("%.3f ", PostSamp.mustar[0][i*K+j]);
        }
        printf("\n");
    }*/
    
    gsl_matrix * mu_prec = gsl_matrix_alloc(K, K);
    gsl_vector * mu_mean = gsl_vector_alloc(K);
    gsl_matrix * mu_samp = gsl_matrix_alloc(1, K);

    for(l=0; l<L; l++){
        Ml = BF.Ml[l];
        
        // Compute posterior covariance matrix
        for(i=0; i<K; i++){
            err = PostSamp.err2inv_eps[0] / PostSamp.Phi2Inv[i*L+l] + PostSamp.err2inv_mu;
            for(j=0; j<K; j++){
                if(i==j){
                    gsl_matrix_set(mu_prec, i, i, err);
                }else{
                    gsl_matrix_set(mu_prec, i, j, 0.0);
                }
            }
        }
        
        // Cholesky Decomposition
        int status = gsl_linalg_cholesky_decomp(mu_prec);
        if(status){
            Rprintf("Failed Cholesky Decomposition for Mu\n");
            goto stop;
        }
        
        // For i
        for(i=0; i<nobs; i++){
            // meanM
            for(k=0; k<K; k++){
                double tempsum2 = 0.0;
                for(rr=0; rr<Ml; rr++){
                    tempsum2 += (double)PostSamp.alphastar[l][k*Ml+rr] * (double)PostSamp.rxb[l][i*Ml + rr] ;
                }
                tempsum2 = PostSamp.latentstar[l][i*K+k] - tempsum2;
                gsl_vector_set(mu_mean, k, tempsum2 / PostSamp.Phi2Inv[k*L+l]*PostSamp.err2inv_eps[0]);
            }
            // Gibbs Sampling for subject i
            MVNsamp2(mu_samp, 1, mu_mean, mu_prec, K, r);
            for(q=0; q<K; q++){
                PostSamp.mustar[l][i*K+q] = gsl_matrix_get(mu_samp, 0, q);
            }
        }
        
    stop:
        ;
        
    }

    /*printf("\nNew_MuStar\n");
    for(i=0; i<3; i++){
        for(j=0; j<K; j++){
            printf("%.3f ", PostSamp.mustar[0][i*K+j]);
        }
        printf("\n");
    }*/

    gsl_matrix_free(mu_prec);
    gsl_vector_free(mu_mean);
    gsl_matrix_free(mu_samp);
    
}
