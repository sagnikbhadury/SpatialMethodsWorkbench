//////////////////////////////////////////////////////////
///// Structure for saveing post samplings /////
//////////////////////////////////////////////////////////

struct Sampling{
    
    float ** basis2;
    float ** zb;
    float ** zb2;
    float ** xb;
    float ** xb_test;
    float ** rxb;
    
    // Intercept for outcomes
    float ** Zinter;
    float ** mean_Zinter;

    // Theta
    float ** theta;
    float ** mean_theta;
    
    // Loading
    float ** load;
    float ** loadstar;
    float ** mean_load;
    float ** mean_loadstar;
    
    // Latent
    float ** latent;
    float ** latentstar;
    float ** latent2;
    float ** mean_latent;
    float ** mean_latentstar;
    
    // Phi
    float * Phi2Inv;
    float * mean_Phi2Inv;
    float Phi2Inv_a;
    float Phi2Inv_b;
    
    // Beta
    float ** beta;
    float ** betastar;
    float ** mean_beta;
    float ** mean_betastar;
    
    // Alpha
    float ** alpha;
    float ** alphastar;
    float ** mean_alpha;
    float ** mean_alphastar;
    
    // XBA, XBAR
    float ** xba; // X*Baisi*Aalpha
    float ** xbar; // X*Baisi*Aalpha*Gamma
    
    // Mu
    float ** mustar;
    float ** mean_mu;
    
    // Summarized image
    float ** sumX;
    
    // Predictor parameters
    float * gamma;
    float * mean_gamma;
    float * mean_gamma_update;
    float * w;
    float * mean_w;
    float w_a;
    float w_b;
    
    // Errors
    float * err2inv_e;
    float * err2inv_e_mean;
    float * err2inv_zeta;
    float * err2inv_zeta_mean;
    float * err2inv_eps;
    float * err2inv_eps_mean;
    float * err2inv_u;
    float * err2inv_u_mean;
    
    // prior of errors
    float err2inv_e_a;
    float err2inv_e_b;
    float err2inv_zeta_a;
    float err2inv_zeta_b;
    float err2inv_eps_a;
    float err2inv_eps_b;
    float err2inv_alpha;
    float err2inv_load;
    float err2inv_mu;
    float err2inv_u_a;
    float err2inv_u_b;
    
    // Testing Set
    float ** Zinter_test;
    float ** latent_test;
    float ** theta_test;
    float ** outcome_test;
    float ** Zinter_test_mean;
    float ** latent_test_mean;
    float ** theta_test_mean;
    float ** outcome_test_mean;
    float ** outcome_test_mean2;
    
    // Training Set
    float ** latent_train;
    float ** theta_train;
    float ** outcome_train;
    float ** latent_train_mean;
    float ** theta_train_mean;
    float ** outcome_train_mean;
    float ** outcome_train_mean2;
    float ** latent_train_err;
    float ** theta_train_err;
};


struct Sampling setupSamp(const int M, const int nobs, const int nts, const int L, const int P, const int K, struct BasisFunc BF, struct Inputdata data){
    
    struct Sampling PostSamp;
    int l, Ml, parcel_len;
    
    ///// Allocate Memory /////
    // data
    PostSamp.basis2 = (float **)malloc(L * sizeof(float *));
    PostSamp.zb = (float **)malloc(L * sizeof(float *));
    PostSamp.zb2 = (float **)malloc(L * sizeof(float *));
    PostSamp.xb = (float **)malloc(L * sizeof(float *));
    PostSamp.rxb = (float **)malloc(L * sizeof(float *));
    for(l=0; l<L; l++){
        Ml = BF.Ml[l];
        PostSamp.basis2[l] = (float *)malloc(Ml * Ml * sizeof(float));
        PostSamp.zb[l] = (float *)malloc(nobs * Ml * sizeof(float));
        PostSamp.zb2[l] = (float *)malloc(nobs * Ml * sizeof(float));
        PostSamp.xb[l] = (float *)malloc(nobs * P * Ml * sizeof(float));
        PostSamp.rxb[l] = (float *)malloc(nobs * Ml * sizeof(float));
    }
    if(nts > 0){
        PostSamp.xb_test = (float **)malloc(L * sizeof(float *));
        for(l=0; l<L; l++){
            Ml = BF.Ml[l];
            PostSamp.xb_test[l] = (float *)malloc(nts*P*Ml* sizeof(float));
        }
    }
    
    
    
    
    // Parameter
    PostSamp.Zinter = (float **)malloc(L * sizeof(float *));
    PostSamp.theta = (float **)malloc(L * sizeof(float *));
    PostSamp.load = (float **)malloc(L * sizeof(float *));
    PostSamp.loadstar = (float **)malloc(L * sizeof(float *));
    PostSamp.latent = (float **)malloc(L * sizeof(float *));
    PostSamp.latentstar = (float **)malloc(L * sizeof(float *));
    PostSamp.latent2 = (float **)malloc(L * sizeof(float *));
    PostSamp.beta = (float **)malloc(L * sizeof(float *));
    PostSamp.betastar = (float **)malloc(L * sizeof(float *));
    PostSamp.alpha = (float **)malloc(L * sizeof(float *));
    PostSamp.alphastar = (float **)malloc(L * sizeof(float *));
    PostSamp.xba = (float **)malloc(L * sizeof(float *));
    PostSamp.xbar = (float **)malloc(L * sizeof(float *));
    PostSamp.mustar = (float **)malloc(L * sizeof(float *));
    PostSamp.sumX = (float **)malloc(L * sizeof(float *));
    for(l=0; l<L; l++){
        Ml = BF.Ml[l];
        parcel_len = data.parcel_len[l];
        PostSamp.Zinter[l] = (float *)malloc(parcel_len * sizeof(float));
        PostSamp.theta[l] = (float *)malloc(Ml * nobs * sizeof(float));
        PostSamp.load[l] = (float *)malloc(Ml * K * sizeof(float));
        PostSamp.loadstar[l] = (float *)malloc(Ml * K * sizeof(float));
        PostSamp.latent[l] = (float *)malloc(nobs * K * sizeof(float));
        PostSamp.latentstar[l] = (float *)malloc(nobs * K * sizeof(float));
        PostSamp.latent2[l] = (float *)malloc(K * K * sizeof(float));
        PostSamp.beta[l] = (float *)malloc(parcel_len * K * sizeof(float));
        PostSamp.betastar[l] = (float *)malloc(parcel_len * K * sizeof(float));
        PostSamp.alpha[l] = (float *)malloc(Ml * K * sizeof(float));
        PostSamp.alphastar[l] = (float *)malloc(Ml * K * sizeof(float));
        PostSamp.xba[l] = (float *)malloc(nobs*K*P * sizeof(float));
        PostSamp.xbar[l] = (float *)malloc(nobs*K * sizeof(float));
        PostSamp.mustar[l] = (float *)malloc(nobs*K * sizeof(float));
        PostSamp.sumX[l] = (float *)malloc(nobs*parcel_len * sizeof(float));
    }
    PostSamp.Phi2Inv = (float *)calloc(K*L, sizeof(float));
    PostSamp.gamma = (float *)calloc(P*L, sizeof(float));
    PostSamp.w = (float *)calloc(L, sizeof(float));
    PostSamp.err2inv_e = (float *)calloc(1, sizeof(float));
    PostSamp.err2inv_zeta = (float *)calloc(1, sizeof(float));
    PostSamp.err2inv_eps = (float *)calloc(1, sizeof(float));
    PostSamp.err2inv_u = (float *)calloc(L, sizeof(float));
    
    
    // Posterior Mean
    PostSamp.mean_Zinter = (float **)malloc(L * sizeof(float *));
    PostSamp.mean_theta = (float **)malloc(L * sizeof(float *));
    PostSamp.mean_load = (float **)malloc(L * sizeof(float *));
    PostSamp.mean_loadstar = (float **)malloc(L * sizeof(float *));
    PostSamp.mean_latent = (float **)malloc(L * sizeof(float *));
    PostSamp.mean_latentstar = (float **)malloc(L * sizeof(float *));
    PostSamp.mean_beta = (float **)malloc(L * sizeof(float *));
    PostSamp.mean_betastar = (float **)malloc(L * sizeof(float *));
    PostSamp.mean_alpha = (float **)malloc(L * sizeof(float *));
    PostSamp.mean_alphastar = (float **)malloc(L * sizeof(float *));
    PostSamp.mean_mu = (float **)malloc(L * sizeof(float *));
    for(l=0; l<L; l++){
        Ml = BF.Ml[l];
        parcel_len = data.parcel_len[l];
        PostSamp.mean_Zinter[l] = (float *)malloc(parcel_len * sizeof(float));
        PostSamp.mean_theta[l] = (float *)malloc(Ml * nobs * sizeof(float));
        PostSamp.mean_load[l] = (float *)malloc(Ml * K * sizeof(float));
        PostSamp.mean_loadstar[l] = (float *)malloc(Ml * K * sizeof(float));
        PostSamp.mean_latent[l] = (float *)malloc(nobs * K * sizeof(float));
        PostSamp.mean_latentstar[l] = (float *)malloc(nobs * K * sizeof(float));
        PostSamp.mean_beta[l] = (float *)malloc(parcel_len * K * sizeof(float));
        PostSamp.mean_betastar[l] = (float *)malloc(parcel_len * K * sizeof(float));
        PostSamp.mean_alpha[l] = (float *)malloc(Ml * K * sizeof(float));
        PostSamp.mean_alphastar[l] = (float *)malloc(Ml * K * sizeof(float));
        PostSamp.mean_mu[l] = (float *)malloc(nobs * K * sizeof(float));
    }
    PostSamp.mean_Phi2Inv = (float *)calloc(K*L, sizeof(float));
    PostSamp.mean_gamma = (float *)calloc(P*L, sizeof(float));
    PostSamp.mean_gamma_update = (float *)calloc(P*L, sizeof(float));
    PostSamp.mean_w = (float *)calloc(L, sizeof(float));
    PostSamp.err2inv_e_mean = (float *)calloc(1, sizeof(float));
    PostSamp.err2inv_zeta_mean = (float *)calloc(1, sizeof(float));
    PostSamp.err2inv_eps_mean = (float *)calloc(1, sizeof(float));
    PostSamp.err2inv_u_mean = (float *)calloc(L, sizeof(float));
    
    ///// Priors
    PostSamp.err2inv_alpha = 1.0f;
    PostSamp.err2inv_load = 1.0f;
    PostSamp.err2inv_mu = 1.0f;
    PostSamp.Phi2Inv_a = 1.0f;
    PostSamp.Phi2Inv_b = 1.0f;
    PostSamp.err2inv_e_a = 1.0f;
    PostSamp.err2inv_e_b = 1.0f;
    PostSamp.err2inv_zeta_a = 1.0f;
    PostSamp.err2inv_zeta_b = 1.0f;
    PostSamp.err2inv_eps_a = 1.0f;
    PostSamp.err2inv_eps_b = 1.0f;
    PostSamp.err2inv_u_a = 1.0f;
    PostSamp.err2inv_u_b = 1.0f;
    PostSamp.w_a = 1.0f;
    PostSamp.w_b = 1.0f;
    
    
    ///// Testing
    if(nts > 0){
        PostSamp.Zinter_test = (float **)malloc(L * sizeof(float *));
        PostSamp.latent_test = (float **)malloc(L * sizeof(float *));
        PostSamp.theta_test = (float **)malloc(L * sizeof(float *));
        PostSamp.outcome_test = (float **)malloc(L * sizeof(float *));
        PostSamp.Zinter_test_mean = (float **)malloc(L * sizeof(float *));
        PostSamp.latent_test_mean = (float **)malloc(L * sizeof(float *));
        PostSamp.theta_test_mean = (float **)malloc(L * sizeof(float *));
        PostSamp.outcome_test_mean = (float **)malloc(L * sizeof(float *));
        PostSamp.outcome_test_mean2 = (float **)malloc(L * sizeof(float *));
        for(l=0; l<L; l++){
            Ml = BF.Ml[l];
            parcel_len = data.parcel_len[l];
            PostSamp.Zinter_test[l] = (float *)malloc(parcel_len * sizeof(float));
            PostSamp.latent_test[l] = (float *)malloc(nts * K * sizeof(float));
            PostSamp.theta_test[l] = (float *)malloc(Ml * nts * sizeof(float));
            PostSamp.outcome_test[l] = (float *)malloc(parcel_len * nts * sizeof(float));
            PostSamp.Zinter_test_mean[l] = (float *)malloc(parcel_len * sizeof(float));
            PostSamp.latent_test_mean[l] = (float *)malloc(nts * K * sizeof(float));
            PostSamp.theta_test_mean[l] = (float *)malloc(Ml * nts * sizeof(float));
            PostSamp.outcome_test_mean[l] = (float *)malloc(parcel_len * nts * sizeof(float));
            PostSamp.outcome_test_mean2[l] = (float *)malloc(parcel_len * nts * sizeof(float));
        }
    }
    
    
    ///// Training
    PostSamp.latent_train = (float **)malloc(L * sizeof(float *));
    PostSamp.theta_train = (float **)malloc(L * sizeof(float *));
    PostSamp.outcome_train = (float **)malloc(L * sizeof(float *));
    PostSamp.latent_train_mean = (float **)malloc(L * sizeof(float *));
    PostSamp.theta_train_mean = (float **)malloc(L * sizeof(float *));
    PostSamp.outcome_train_mean = (float **)malloc(L * sizeof(float *));
    PostSamp.outcome_train_mean2 = (float **)malloc(L * sizeof(float *));
    PostSamp.latent_train_err = (float **)malloc(L * sizeof(float *));
    PostSamp.theta_train_err = (float **)malloc(L * sizeof(float *));
    for(l=0; l<L; l++){
        Ml = BF.Ml[l];
        parcel_len = data.parcel_len[l];
        PostSamp.latent_train[l] = (float *)malloc(nobs * K * sizeof(float));
        PostSamp.theta_train[l] = (float *)malloc(Ml * nobs * sizeof(float));
        PostSamp.outcome_train[l] = (float *)malloc(parcel_len * nobs * sizeof(float));
        PostSamp.latent_train_mean[l] = (float *)malloc(nobs * K * sizeof(float));
        PostSamp.theta_train_mean[l] = (float *)malloc(Ml * nobs * sizeof(float));
        PostSamp.outcome_train_mean[l] = (float *)malloc(parcel_len * nobs * sizeof(float));
        PostSamp.outcome_train_mean2[l] = (float *)malloc(parcel_len * nobs * sizeof(float));
        PostSamp.latent_train_err[l] = (float *)malloc(nobs * K * sizeof(float));
        PostSamp.theta_train_err[l] = (float *)malloc(Ml * nobs * sizeof(float));
    }
    
    // Compute basis^T*basis
    int i, j, k, p;
    double tempsum;
    
    for(l=0; l<L; l++){
        //printf("l=%d\n", l);
        Ml = BF.Ml[l];
        parcel_len = data.parcel_len[l];
        
        for(i=0; i<Ml; i++){
            for(j=0; j<Ml; j++){
                tempsum = 0.0;
                for(k=0; k<parcel_len; k++){
                    tempsum += (double)BF.basis[l][k*Ml+i] * (double)BF.basis[l][k*Ml+j];
                }
                PostSamp.basis2[l][i*Ml+j] = tempsum;
            }
        }
    }
    
    
    // Compute z_parcel*basis
    for(l=0; l<L; l++){
        Ml = BF.Ml[l];
        parcel_len = data.parcel_len[l];
        for(i=0; i<nobs; i++){
            for(j=0; j<Ml; j++){
                tempsum=0.0;
                for(k=0; k<parcel_len; k++){
                    tempsum += (double)data.Zl[l][i*parcel_len+k] * (double)BF.basis[l][k*Ml+j];
                }
                PostSamp.zb[l][i*Ml+j] = tempsum;
            }
        }
    }
    
    // Compute basis * X_{ip}
    for(i=0; i<nobs; i++){

        for(l=0; l<L; l++){
            Ml = BF.Ml[l];
            parcel_len = data.parcel_len[l];
            
            for(p=0; p<P; p++){
                for(k=0; k<Ml; k++){
                    tempsum = 0.0;
                    for(j=0; j<parcel_len; j++){
                        tempsum += (double)data.Xl[l][i*parcel_len*P+p*parcel_len+j] * (double)BF.basis[l][j*Ml+k];
                    }
                    PostSamp.xb[l][i*P*Ml + p*Ml+ k] = tempsum;
                    
                }
            }
        }
    }
    
    // Compute basis * X_test
    if(nts > 0){
        for(i=0; i<nts; i++){
            for(l=0; l<L; l++){
                Ml = BF.Ml[l];
                parcel_len = data.parcel_len[l];
                
                for(p=0; p<P; p++){
                    for(k=0; k<Ml; k++){
                        
                        tempsum = 0.0;
                        for(j=0; j<parcel_len; j++){
                            tempsum += (double)data.Xl_test[l][i*parcel_len*P+p*parcel_len+j] * (double)BF.basis[l][j*Ml+k];
                        }
                        PostSamp.xb_test[l][i*P*Ml + p*Ml+ k] = tempsum;
                    }
                }
            }
        }
    }
    

    return PostSamp;
}
