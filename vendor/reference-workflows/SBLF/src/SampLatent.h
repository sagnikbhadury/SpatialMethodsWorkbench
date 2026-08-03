//
// Post sampling of latent factors
//

void latent_samp(const int L, const int K, const int nobs, struct Sampling PostSamp, const struct Inputdata data, const struct BasisFunc BF, gsl_rng *r){

    int l, i, j, m, k, q, rr, Ml;
   gsl_set_error_handler_off();
    /*printf("\n#########\nOld_LatentStar\n");
    for(i=0; i<3; i++){
        for(j=0; j<K; j++){
            printf("%.3f ", PostSamp.latentstar[0][i*K+j]);
        }
        printf("\n");
    }*/

    gsl_matrix * latent_prec = gsl_matrix_alloc(K, K);
    gsl_vector * latent_mean = gsl_vector_alloc(K);
    gsl_matrix * latent_samp = gsl_matrix_alloc(1, K);
    

    for(l=0; l<L; l++){
        Ml =  BF.Ml[l];
        
        // CovInv
        for(i=0; i<K; i++){
            for(j=0; j<K; j++){
                double tempsum = 0.0;
                for(m=0; m<Ml; m++){
                    tempsum += (double)PostSamp.loadstar[l][m*K+i] * (double)PostSamp.loadstar[l][m*K+j]* (double)PostSamp.err2inv_zeta[0];
                }
                if(i==j){
                    tempsum += PostSamp.err2inv_eps[0]/PostSamp.Phi2Inv[j*L+l];
                }
                gsl_matrix_set(latent_prec, i, j, tempsum);
            }
        }
        
        // Cholesky Decomposition
        int status = gsl_linalg_cholesky_decomp(latent_prec);
        if(status){
            Rprintf("Failed Cholesky Decomposition for Latent\n");
            goto stop;
        }
      
        // For each i
        for(i=0; i<nobs; i++){
            
            // meanM
            for(k=0; k<K; k++){
                double tempsum = 0.0;
                for(j=0; j<Ml; j++){
                    tempsum += (double)PostSamp.theta[l][i*Ml+j] * (double)PostSamp.loadstar[l][j*K+k] *(double)PostSamp.err2inv_zeta[0];
                }
                double tempsum2 = 0.0;
                for(rr=0; rr<Ml; rr++){
                    tempsum2 += (double) PostSamp.alphastar[l][k*Ml+rr] * (double)PostSamp.rxb[l][i*Ml+rr];
                }
                tempsum2 += PostSamp.mustar[l][i*K+k];
                tempsum +=  tempsum2 / PostSamp.Phi2Inv[k*L+l]*PostSamp.err2inv_eps[0];
                gsl_vector_set(latent_mean, k, tempsum);
            }
            // Gibbs Sampling for subject i
            MVNsamp2(latent_samp, 1, latent_mean, latent_prec, K, r);
            for(q=0; q<K; q++){
                PostSamp.latentstar[l][i*K+q] = gsl_matrix_get(latent_samp, 0, q);
            }
        }
        
    stop:
        ;
        
    }

    // Transfer back to latent
    for(l=0; l<L; l++){
        Ml = BF.Ml[l];
        for(k=0; k<K; k++){
            double diag = PostSamp.loadstar[l][k*K + k];
            double sign;
            if(diag >=0){
                sign = 1.0;
            }else{
                sign = -1.0;
            }
            double phi = 1.0/sqrt((double)PostSamp.Phi2Inv[k*L+l]);
            for(i=0; i<nobs; i++){
                PostSamp.latent[l][i*K+k] = (PostSamp.latentstar[l][i*K+k] - PostSamp.mustar[l][i*K+k])* sign * phi;
            }
        }
    }
    
    /*printf("\nNew_LatentStar\n");
    for(i=0; i<3; i++){
        for(j=0; j<K; j++){
            printf("%.3f ", PostSamp.latentstar[0][i*K+j]);
        }
        printf("\n");
    }*/

    gsl_matrix_free(latent_prec);
    gsl_vector_free(latent_mean);
    gsl_matrix_free(latent_samp);
    
}
