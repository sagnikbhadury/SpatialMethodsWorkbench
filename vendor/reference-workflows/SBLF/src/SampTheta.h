//
// Post sampling of theta
//

void theta_samp(int nobs, int L, int K, struct Sampling PostSamp, struct Inputdata data, const struct BasisFunc BF,
                gsl_rng *r, int * singular){
    gsl_set_error_handler_off();
    int i, j, l, k, m, rr, Ml, loc;
    double tmp;
    
    /*printf("\n########\nOld_Theta:\n");
    for(i=0; i<3; i++){
        for(j=0; j<3; j++){
            printf("%.3f ", PostSamp.theta[0][i*BF.Ml[0]+j]);
        }
        printf("\n");
    }*/
    
    
    // Gibbs Sampling
    for(l=0; l<L; l++){
        Ml = BF.Ml[l];
        
        gsl_matrix * theta_prec = gsl_matrix_alloc(Ml, Ml);
        gsl_vector * theta_mean = gsl_vector_alloc(Ml);
        gsl_matrix * theta_samp = gsl_matrix_alloc(1, Ml);

        // compute precision matrix
        for(i=0; i<Ml; i++){
            for(j=0; j<Ml; j++){
                loc = i*Ml + j;
                if(i==j){
                    tmp = PostSamp.basis2[l][loc] * PostSamp.err2inv_e[0] + PostSamp.err2inv_zeta[0];
                }else{
                    tmp = PostSamp.basis2[l][loc] * PostSamp.err2inv_e[0];
                }
                gsl_matrix_set(theta_prec, i, j, tmp);
            }
        }
        
        // Cholesky decomposition
        int status = gsl_linalg_cholesky_decomp(theta_prec);
        if(status){
            Rprintf("Failed Cholesky Decomposition for Theta\n");
            goto stop;
        }
        
        for(k=0; k<nobs; k++){
            // Compute posterior mean of theta
            for(m=0; m<Ml; m++){
                double tempsum = 0.0;
                for(rr=0; rr<K; rr++){
                    tempsum += (double)PostSamp.latentstar[l][k*K+rr] * (double)PostSamp.loadstar[l][m*K+rr];
                }
                tmp = tempsum*PostSamp.err2inv_zeta[0] + (double)PostSamp.zb2[l][k*Ml+m]*PostSamp.err2inv_e[0];
                gsl_vector_set(theta_mean, m, tmp);
            }

            // Sampling
            MVNsamp2(theta_samp, 1, theta_mean, theta_prec, Ml, r);
            for(m=0; m<Ml; m++){
                PostSamp.theta[l][k*Ml+m] = gsl_matrix_get(theta_samp, 0, m);
            }

        }

    stop:
        ;
        
        gsl_vector_free(theta_mean);
        gsl_matrix_free(theta_samp);
        gsl_matrix_free(theta_prec);
        
    }
    
    /*printf("\nNew_Theta:\n");
    for(i=0; i<3; i++){
        for(j=0; j<3; j++){
            printf("%.3f ", PostSamp.theta[0][i*BF.Ml[0]+j]);
        }
        printf("\n");
    }*/
    
}



