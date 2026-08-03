//
// Post sampling of loading matrix
//
void load_samp(int K, int L, int nobs,
               struct Sampling PostSamp, const struct BasisFunc BF,
               gsl_rng *r){
   
    int l, m, i, j, k, Ml;
    gsl_set_error_handler_off();
    /*printf("\n#########\nOld_LoadingStar\n");
    for(i=0; i<5; i++){
        for(j=0; j<K; j++){
            printf("%.3f ", PostSamp.loadstar[0][i*K+j]);
        }
        printf("\n");
    }*/
    
    gsl_matrix * load_prec;
    gsl_vector * load_mean;
    gsl_matrix * load_samp;
    
    load_prec = gsl_matrix_alloc(K, K);
    load_mean = gsl_vector_alloc(K);
    load_samp = gsl_matrix_alloc(1, K);
    
    // latent^T*latent
    for(l=0; l<L; l++){
        for(k=0; k<K; k++){
            for(j=0; j<K; j++){
                double temp = 0.0;
                for(i=0; i<nobs; i++) {
                    temp += (double)PostSamp.latentstar[l][i*K+k] * (double)PostSamp.latentstar[l][i*K+j];
                }
                PostSamp.latent2[l][k*K+j] = temp;
            }
        }
    }

    // Gibbs Sampling
    for(l=0; l<L; l++){
        Ml = BF.Ml[l];
        
        for(m=0; m<Ml; m++){
            
            if(m<(K-1)){
                if(m==0){
                    double vr = 1.0 / (PostSamp.latent2[l][0]*PostSamp.err2inv_zeta[0] + PostSamp.err2inv_load);
                    double mn = 0;
                    for(i=0; i<nobs; i++){
                        mn += (double)PostSamp.latentstar[l][i*K] * (double)PostSamp.theta[l][i*Ml];
                    }
                    
                    mn = mn * (double)PostSamp.err2inv_zeta[0] * vr;
                    PostSamp.loadstar[l][0] = gsl_ran_gaussian(r, sqrt(vr)) + mn;
                    
                    for(k=1; k<K; k++){
                        PostSamp.loadstar[l][k] = 0.0;
                    }
                }else{
                    int Km = m+1;
                    gsl_matrix * load_prec_m = gsl_matrix_alloc(Km, Km);
                    gsl_vector * load_mean_m = gsl_vector_alloc(Km);
                    gsl_matrix * load_samp_m = gsl_matrix_alloc(1, Km);
                    // covInv
                    for(k=0; k<Km; k++){
                        for(j=0; j<Km; j++){
                            double temp;
                            temp = PostSamp.latent2[l][k*K+j] * PostSamp.err2inv_zeta[0];
                            if(k==j){
                                temp += PostSamp.err2inv_load;
                            }
                            gsl_matrix_set(load_prec_m, k, j, temp);
                        }
                    }
                    // meanM
                    for(k=0; k<Km; k++){
                        double temp1 = 0.0;
                        for(i=0; i<nobs; i++){
                            temp1 += (double)PostSamp.latentstar[l][i*K+k] * (double)PostSamp.theta[l][i*Ml+m];
                        }
                        temp1 *= (double)PostSamp.err2inv_zeta[0];
                        gsl_vector_set(load_mean_m, k, temp1);
                    }
                    // Gibb Sampler
                    MVNsamp(load_samp_m, 1, load_mean_m, load_prec_m, Km, r);
                    for(k=0; k<K; k++){
                        if(k<Km){
                            PostSamp.loadstar[l][m*K+k] = gsl_matrix_get(load_samp_m, 0, k);
                        }else{
                            PostSamp.loadstar[l][m*K+k] = 0.0f;
                        }
                    }
                    // Free
                    gsl_matrix_free(load_prec_m);
                    gsl_vector_free(load_mean_m);
                    gsl_matrix_free(load_samp_m);
                }
            }else{ // m>=(K-1)
                // covInv
                double sig = PostSamp.err2inv_zeta[0];
                for(k=0; k<K; k++){
                    for(j=0; j<K; j++){
                        double temp;
                        temp = PostSamp.latent2[l][k*K+j] * sig;
                        if(k==j){
                            temp += PostSamp.err2inv_load;
                        }
                        gsl_matrix_set(load_prec, k, j, temp);
                    }
                }
                // Posterior mean vector
                for(k=0; k<K; k++){
                    double temp1 = 0.0;
                    for(i=0; i<nobs; i++){
                        temp1 += (double)PostSamp.latentstar[l][i*K+k] * (double)PostSamp.theta[l][i*Ml+m];
                    }
                    temp1 *= (double)PostSamp.err2inv_zeta[0];
                    gsl_vector_set(load_mean, k, temp1);
                }

                // Gibbs Sampling
                MVNsamp(load_samp, 1, load_mean, load_prec, K, r);
                for(k=0; k<K; k++){
                    PostSamp.loadstar[l][m*K+k] = gsl_matrix_get(load_samp, 0, k);
                }
            }
        }
    }
    
    // Transfer to loading matrix
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
            double phi = sqrt((double)PostSamp.Phi2Inv[k*L+l]);
            for(m=0; m<Ml; m++){
                PostSamp.load[l][m*K+k] = sign * phi * PostSamp.loadstar[l][m*K+k];
            }
        }
    }
    
    /*printf("\nNew_LoadingStar\n");
    for(i=0; i<5; i++){
        for(j=0; j<K; j++){
            printf("%.3f ", PostSamp.loadstar[0][i*K+j]);
        }
        printf("\n");
    }*/
    
    gsl_vector_free(load_mean);
    gsl_matrix_free(load_samp);
    gsl_matrix_free(load_prec);
}


