//
// Post sampling of gamma
//

void gamma_samp(const int P, const int nobs, const int K, const int L, struct Sampling PostSamp, const struct Inputdata data, const struct BasisFunc BF, gsl_rng *r){
    gsl_set_error_handler_off();
    int i, k, p, l, m, Ml;
    double tempsum, tempsum1, temp, tau1, tau2, tau3, tau4, p0, p1;
    
    /*printf("\n#########\nOld_Gamma\n");
    for(i=0; i<10; i++){
        for(l=0; l<L; l++){
            printf("%.3f ", PostSamp.gamma[i*L+l]);
        }
    }
    printf("\n");
    */
    
    // Tau_ikpl
    float * tau = (float *)calloc(nobs*K*L*P, sizeof(float));
    for(i=0; i<nobs; i++){
        for(l=0; l<L; l++){
            for(k=0; k<K; k++){
                for(p=0; p<P; p++){
                    tempsum = PostSamp.latentstar[l][i*K+k] - PostSamp.mustar[l][i*K+k] - PostSamp.xbar[l][i*K+k]+ PostSamp.xba[l][i*P*K+p*K+k]*PostSamp.gamma[p*L+l];
                    tau[i*L*K*P + p*K*L + l*K +k ] = tempsum;
                }
            }
        }
    }
    
    for(l=0; l<L; l++){
        for(p=0; p<P; p++){
            
            // term_0
            tempsum = 0.0;
            tempsum1 = 0.0;
            for(i=0; i<nobs; i++){
                for(k=0; k<K; k++){
                    tau1 = tau[i*L*K*P+p*K*L+l*K+k];
                    tau2 = tau1 * tau1;
                    tau3 = tau1 - PostSamp.xba[l][i*P*K+p*K+k];
                    tau4 = tau3 * tau3;
                    temp = PostSamp.err2inv_eps[l] / PostSamp.Phi2Inv[k*L+l];
                    tempsum += temp * tau2;
                    tempsum1 += temp * tau4;
                }
            }
            
            p0 = -0.5*tempsum+ log(1.0 - (double)PostSamp.w[l]);
            p1 = -0.5*tempsum1+ log((double)PostSamp.w[l]);
            p1 = 1.0 / (1 + exp(p0-p1));
            
            PostSamp.gamma[p*L+l] = gsl_ran_bernoulli(r, p1);

        }
    }
    
    /*printf("\nNew_Gamma\n");
    for(i=0; i<10; i++){
        for(l=0; l<L; l++){
            printf("%.3f ", PostSamp.gamma[i*L+l]);
        }
    }
    printf("\n");*/
    
    // Update gamma * X * Basis
    for(i=0; i<nobs; i++){
        for(l=0; l<L; l++){
            Ml = BF.Ml[l];
            for(m=0; m<Ml; m++){
                tempsum = 0.0;
                for(p=0; p<P; p++){
                    tempsum += PostSamp.gamma[p*L+l] * PostSamp.xb[l][i*P*Ml+p*Ml+m];
                }
                PostSamp.rxb[l][i*Ml+m] = tempsum;
            }
        }
    }
    
    // Beta * X * Gamma
    for(i=0; i<nobs; i++){
        for(l=0; l<L; l++){
            for(k=0; k<K; k++){
                tempsum = 0.0;
                for(p=0; p<P; p++){
                    tempsum += PostSamp.xba[l][i*P*K+p*K+k] * PostSamp.gamma[p*L+l];
                }
                PostSamp.xbar[l][i*K+k] = tempsum;
            }
        }
    }
    
    free(tau);
}


