//
// Post sampling of alpha
//

void alpha_samp(const int nobs, const int K, const int P, const int L,
               struct Sampling PostSamp, const struct Inputdata data, const struct BasisFunc BF, gsl_rng *r){
    
    gsl_set_error_handler_off();
    int i, j, l, k, m, p, Ml;
    
    /*printf("\n#########\nOld_AlphaStar\n");
    for(i=0; i<K; i++){
        for(j=0; j<5; j++){
            printf("%.3f ", PostSamp.alphastar[0][i*BF.Ml[0]+j]);
        }
        printf("\n");
     }*/

    double temp;
    double tempsum;

    for(l=0; l<L; l++){
        Ml = BF.Ml[l];
        float * tempmm = (float *)calloc(Ml*Ml, sizeof(float));
        gsl_vector * alpha_mean = gsl_vector_alloc(Ml);
        gsl_matrix * alpha_samp = gsl_matrix_alloc(1, Ml);
        gsl_matrix * alpha_prec = gsl_matrix_alloc(Ml, Ml);

        // rxb^2
        for(m=0; m<Ml; m++){
            for(j=0; j<Ml; j++){
                tempmm[m*Ml+j] = 0.0;
            }
        }
        for(i=0; i<nobs; i++){
            for(m=0; m<Ml; m++){
                for(j=0; j<Ml; j++){
                    tempmm[m*Ml+j] += PostSamp.rxb[l][i*Ml+m] * PostSamp.rxb[l][i*Ml+j];
                }
            }
        }

        for(k=0; k<K; k++){
            
           // covInv
            for(m=0; m<Ml; m++){
                for(j=0; j<Ml; j++){
                    temp =  tempmm[m*Ml+j] * PostSamp.err2inv_eps[0] / PostSamp.Phi2Inv[k*L+l];
                    if(m==j){
                        temp += PostSamp.err2inv_alpha;
                    }
                    gsl_matrix_set(alpha_prec, m, j, temp);
                }
            }

            // Cholesky Decomposation
            int status = gsl_linalg_cholesky_decomp(alpha_prec);
            if(status){
                Rprintf("Failed Cholesky Decomposition for Alpha\n");
                goto stop;
            }
            
            // meanM
            for(i=0; i<Ml; i++){
                tempsum = 0.0;
                for(j=0; j<nobs; j++){
                    tempsum += (PostSamp.latentstar[l][j*K+k] - PostSamp.mustar[l][j*K+k]) * PostSamp.rxb[l][j*Ml+i];
                }
                
                tempsum *= PostSamp.err2inv_eps[0]/ PostSamp.Phi2Inv[k*L+l];
                gsl_vector_set(alpha_mean, i, tempsum);
            }
            
            // Gibbs sampling
            MVNsamp2(alpha_samp, 1, alpha_mean, alpha_prec, Ml, r);
            for(j=0; j<Ml; j++){
                PostSamp.alphastar[l][k*Ml+j] = gsl_matrix_get(alpha_samp, 0, j);
            }
            
        stop:
            ;
        }
        
        free(tempmm);
        gsl_vector_free(alpha_mean);
        gsl_matrix_free(alpha_samp);
        gsl_matrix_free(alpha_prec);
    }

    // Transfer back to alpha
    for(k=0; k<K; k++){
        for(l=0; l<L; l++){
            Ml = BF.Ml[l];
            double diag = PostSamp.loadstar[l][k*K + k];
            double sign;
            if(diag >=0){
                sign = 1.0;
            }else{
                sign = -1.0;
            }
            double phi = 1.0/sqrt((double)PostSamp.Phi2Inv[k*L+l]);
            
            for(i=0; i<Ml; i++){
                PostSamp.alpha[l][k*Ml+i] = sign * phi * PostSamp.alphastar[l][k*Ml+i];
            }
        }
    }

    /*printf("\nNew_AlphaStar\n");
    for(i=0; i<K; i++){
        for(j=0; j<5; j++){
            printf("%.3f ", PostSamp.alphastar[0][i*BF.Ml[0]+j]);
        }
        printf("\n");
    }*/
    

    // Update: BETA * X_i,p,l
    for(i=0; i<nobs; i++){
        for(p=0; p<P; p++){
            for(l=0; l<L; l++){
                Ml = BF.Ml[l];
                for(k=0; k<K; k++){
                    tempsum = 0.0;
                    for(m=0; m<Ml; m++){
                        tempsum += PostSamp.xb[l][i*P*Ml+p*Ml+m] * PostSamp.alphastar[l][k*Ml+m];
                    }
                    PostSamp.xba[l][i*P*K+p*K+k] = tempsum;
                }
            }
        }
    }
    
}

