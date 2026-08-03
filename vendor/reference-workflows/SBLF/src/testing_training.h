//
// To compute posterior estimations and predictions
//

void est_testing(const int nts, const int nobs, const int P, const int L, const int K,
                 struct Sampling PostSamp, const struct Inputdata data, const struct BasisFunc BF);
void est_training(const int nobs, const int P, const int L, const int K,
                  struct Sampling PostSamp, const struct Inputdata data, const struct BasisFunc BF);


//////////////////////////////////////////////////////////////////////////////////////////////////
void est_testing(const int nts, const int nobs, const int P, const int L, const int K,
                 struct Sampling PostSamp, const struct Inputdata data, const struct BasisFunc BF){
    
    int i, l, k, p, m, Ml, parcel_len;
    double tempsum, tempsum1, tempsum2;
    float ** xba = (float **)malloc(L * sizeof(float *));
    for(l=0; l<L; l++){
        xba[l] = (float *)malloc(nts*K*P * sizeof(float));
    }
    
    // X * beta = X * basis * alpha
    for(i=0; i<nts; i++){
        for(p=0; p<P; p++){
            for(l=0; l<L; l++){
                Ml = BF.Ml[l];
                for(k=0; k<K; k++){
                    tempsum = 0.0;
                    for(m=0; m<Ml; m++){
                        tempsum += PostSamp.xb_test[l][i*P*Ml+p*Ml+m] * PostSamp.alpha[l][k*Ml+m];
                    }
                    xba[l][i*P*K+p*K+k] = tempsum;
                }
            }
        }
    }
    
    // latent factors; E(latent) = \sum_p x*beta*gamma
    for(i=0; i<nts; i++){
        for(l=0; l<L; l++){
            Ml = BF.Ml[l];
            for(k=0; k<K; k++){
                // mean part
                tempsum = 0.0;
                for(p=0; p<P; p++){
                    tempsum += xba[l][i*P*K+p*K+k] * PostSamp.gamma[p*L+l];
                }
                PostSamp.latent_test[l][i*K+k] = tempsum;
            }
        }
    }
    
    // theta; E(theta) = load * latent
    for(i=0; i<nts; i++){
        for(l=0; l<L; l++){
            Ml = BF.Ml[l];
            for(m=0; m<Ml; m++){
                
                // mean
                tempsum = 0.0;
                tempsum1 = 0.0;
                tempsum2 = 0.0;
                for(k=0; k<K; k++){
                    tempsum += PostSamp.latent_test[l][i*K+k] * PostSamp.load[l][m*K+k];
                }
                PostSamp.theta_test[l][i*Ml+m] = tempsum;
            }
        }
    }
    
    // outcome; E(outcome) = theta * basis
    for(i=0; i<nts; i++){
        for(l=0; l<L; l++){
            Ml = BF.Ml[l];
            parcel_len = data.parcel_len[l];
            
            for(p=0; p<parcel_len; p++){
                tempsum = 0.0;
                tempsum1 = 0.0;
                tempsum2 = 0.0;
                for(m=0; m<Ml; m++){
                    tempsum += PostSamp.theta_test[l][i*Ml+m] * BF.basis[l][p*Ml+m];
                }
                PostSamp.outcome_test[l][i*parcel_len + p] = tempsum;
                // Intercept
                PostSamp.Zinter_test[l][p] = PostSamp.Zinter[l][p];
                PostSamp.outcome_test[l][i*parcel_len + p] += PostSamp.Zinter_test[l][p];
            }
        }
    }
    for(l=0; l<L; l++){
        free(xba[l]);
    }
    free(xba);
}




void est_training(const int nobs, const int P, const int L, const int K,
                  struct Sampling PostSamp, const struct Inputdata data, const struct BasisFunc BF){
    
    
    int i, l, m, p, k, Ml, parcel_len;
    double tempsum;
    float ** xba = (float **)malloc(L * sizeof(float *));
    for(l=0; l<L; l++){
        xba[l] = (float *)malloc(nobs*K*P * sizeof(float));
    }
    

    // latent
    for(i=0; i<nobs; i++){
        for(p=0; p<P; p++){
            for(l=0; l<L; l++){
                Ml = BF.Ml[l];
                for(k=0; k<K; k++){
                    tempsum = 0.0;
                    for(m=0; m<Ml; m++){
                        tempsum += PostSamp.xb[l][i*P*Ml+p*Ml+m] * PostSamp.alpha[l][k*Ml+m];
                    }
                    xba[l][i*P*K+p*K+k] = tempsum;
                }
            }
        }
    }
    for(i=0; i<nobs; i++){
        for(l=0; l<L; l++){
            for(k=0; k<K; k++){
                tempsum = 0.0;
                for(p=0; p<P; p++){
                    tempsum += xba[l][i*P*K+p*K+k] * PostSamp.gamma[p*L+l];
                }
                PostSamp.latent_train[l][i*K+k] = tempsum;
                PostSamp.latent_train_err[l][i*K+k] = PostSamp.latent[l][i*K+k] - PostSamp.latent_train[l][i*K+k];
            }
        }
    }
  
    // theta
    for(i=0; i<nobs; i++){
        for(l=0; l<L; l++){
            Ml = BF.Ml[l];
            for(m=0; m<Ml; m++){
                tempsum = 0.0;
                for(k=0; k<K; k++){
                    tempsum += PostSamp.latent_train[l][i*K+k] * PostSamp.load[l][m*K+k];
                }
                PostSamp.theta_train[l][i*Ml+m] = tempsum;
                PostSamp.theta_train_err[l][i*Ml+m] = PostSamp.theta[l][i*Ml+m] - PostSamp.theta_train[l][i*Ml+m];
            }
        }
    }
    
    // outcome
    for(i=0; i<nobs; i++){
        for(l=0; l<L; l++){
            Ml = BF.Ml[l];
            parcel_len = data.parcel_len[l];
            
            for(p=0; p<parcel_len; p++){
                tempsum = 0.0;
                for(m=0; m<Ml; m++){
                    tempsum += PostSamp.theta[l][i*Ml+m] * BF.basis[l][p*Ml+m];
                }
                PostSamp.outcome_train[l][i*parcel_len + p] = PostSamp.Zinter[l][p] + tempsum;
            }
        }
    }
    
    for(l=0; l<L; l++){
        free(xba[l]);
    }
    free(xba);

}

