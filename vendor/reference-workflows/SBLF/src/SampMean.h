//
// Compute and write posterior mean estimations
//

void samp_mean_gamma(int P, int L, struct Sampling PostSamp);
void samp_mean(int L, int K, int nobs, int nts, int P, struct Sampling PostSamp, const struct Inputdata data, const struct BasisFunc BF);
void output_mean(const int L, const int K, const int nobs, const int nts, const int P, const char * outpath, struct Sampling PostSamp, const struct Inputdata data, const struct BasisFunc BF, const int iter);

///////////////////////////////////////////////////////////////////////////////////////////////////
void samp_mean_gamma(int P, int L, struct Sampling PostSamp){
    int p;
    for(p=0; p<L*P; p++){
        PostSamp.mean_gamma_update[p] += PostSamp.gamma[p];
    }
}


void samp_mean(int L, int K, int nobs, int nts, int P,
               struct Sampling PostSamp, const struct Inputdata data, const struct BasisFunc BF){
    
    int i, k, m, p, l, Ml, parcel_len;
    double temp;
    
    // outcome
    if(nts > 0){
        for(i=0; i<nts; i++){
            for(l=0; l<L; l++){
                Ml = BF.Ml[l];
                parcel_len = data.parcel_len[l];
                for(p=0; p<parcel_len; p++){
                    temp =  PostSamp.outcome_test[l][i*parcel_len + p];
                    PostSamp.outcome_test_mean[l][i*parcel_len + p] += temp;
                    PostSamp.outcome_test_mean2[l][i*parcel_len + p] += temp * temp;
                }
            }
        }
    }
    for(i=0; i<nobs; i++){
        for(l=0; l<L; l++){
            Ml = BF.Ml[l];
            parcel_len = data.parcel_len[l];
            for(p=0; p<parcel_len; p++){
                temp = PostSamp.outcome_train[l][i*parcel_len + p];
                PostSamp.outcome_train_mean[l][i*parcel_len + p] += temp;
                PostSamp.outcome_train_mean2[l][i*parcel_len + p] += temp * temp;
            }
        }
    }
    
    // Intercept in submodel 1
    for(l=0; l<L; l++){
        parcel_len = data.parcel_len[l];
        for(p=0; p<parcel_len; p++){
            PostSamp.mean_Zinter[l][p] += PostSamp.Zinter[l][p];
        }
    }
    if(nts > 0){
        for(l=0; l<L; l++){
            parcel_len = data.parcel_len[l];
            for(p=0; p<parcel_len; p++){
                PostSamp.Zinter_test_mean[l][p] += PostSamp.Zinter_test[l][p];
            }
        }
    }
    
    // theta
    for(l=0; l<L; l++){
        Ml = BF.Ml[l];
        for(i=0; i<nobs; i++){
            for(m=0; m<Ml; m++){
                PostSamp.mean_theta[l][i*Ml+m] += PostSamp.theta[l][i*Ml+m] ;
                PostSamp.theta_train_mean[l][i*Ml+m] += PostSamp.theta_train[l][i*Ml+m];
            }
        }
    }
    if(nts > 0){
        for(l=0; l<L; l++){
            Ml = BF.Ml[l];
            for(i=0; i<nts; i++){
                for(m=0; m<Ml; m++){
                    PostSamp.theta_test_mean[l][i*Ml+m] += PostSamp.theta_test[l][i*Ml+m];
                }
            }
        }
    }
    
    // Phi2Inv
    for(k=0; k<K*L; k++){
        PostSamp.mean_Phi2Inv[k] += PostSamp.Phi2Inv[k];
    }
    
    // Loading matrix
    for(l=0; l<L; l++){
        Ml = BF.Ml[l];
        for(m=0; m<Ml; m++){
            for(k=0; k<K; k++){
                PostSamp.mean_load[l][m*K+k] += PostSamp.load[l][m*K+k];
                PostSamp.mean_loadstar[l][m*K+k] += PostSamp.loadstar[l][m*K+k];
            }
        }
        
    }

    // Latent & mu
    for(i=0; i<nobs; i++){
        for(l=0; l<L; l++){
            for(k=0; k<K; k++){
                PostSamp.mean_latentstar[l][i*K+k] += PostSamp.latentstar[l][i*K+k];
                PostSamp.mean_latent[l][i*K+k] += PostSamp.latent[l][i*K+k];
                PostSamp.mean_mu[l][i*K+k] += PostSamp.mustar[l][i*K+k];
                PostSamp.latent_train_mean[l][i*K+k] += PostSamp.latent_train[l][i*K+k];
            }
        }
    }
    if(nts > 0){
        for(i=0; i<nts; i++){
            for(l=0; l<L; l++){
                for(k=0; k<K; k++){
                    PostSamp.latent_test_mean[l][i*K+k] += PostSamp.latent_test[l][i*K+k];
                }
            }
        }
    }

    // beta
    for(l=0; l<L; l++){
        Ml = BF.Ml[l];
        parcel_len = data.parcel_len[l];
        for(p=0; p<parcel_len; p++){
            for(k=0; k<K; k++){
                PostSamp.mean_beta[l][p*K+k] += PostSamp.beta[l][p*K+k];
                PostSamp.mean_betastar[l][p*K+k] += PostSamp.betastar[l][p*K+k];
            }
        }
    }
    
    // alpha
    for(l=0; l<L; l++){
        Ml = BF.Ml[l];
        for(k=0; k<K; k++){
            for(m=0; m<Ml; m++){
                PostSamp.mean_alpha[l][k*Ml+m] += PostSamp.alpha[l][k*Ml+m];
                PostSamp.mean_alphastar[l][k*Ml+m] += PostSamp.alphastar[l][k*Ml+m];
            }
        }
    }

    // gamma, weight and omega2inv
    for(p=0; p<L*P; p++){
        PostSamp.mean_gamma[p] += PostSamp.gamma[p];
    }
    for(p=0; p<L; p++){
        PostSamp.mean_w[p] += PostSamp.w[p];
    }
    
    // Errors
    PostSamp.err2inv_e_mean[0] += PostSamp.err2inv_e[0];
    PostSamp.err2inv_zeta_mean[0] += PostSamp.err2inv_zeta[0];
    PostSamp.err2inv_eps_mean[0] += PostSamp.err2inv_eps[0];
    
}


void output_mean(const int L, const int K, const int nobs, const int nts, const int P, const char * outpath, struct Sampling PostSamp, const struct Inputdata data, const struct BasisFunc BF, const int iter){
    int i, j, k, l, Ml, parcel_len;
    float niter = (float)iter;
    // outcome_test
    if(nts > 0){
        char *outts = (char *)calloc(500,sizeof(char));
        strcat(outts, outpath);
        strcat(outts, "PostMean_Out_test.txt");
        FILE * foutts = fopen(outts, "wb");
        if(foutts == NULL){
            Rprintf("Cannot open file for PostMean_Out_test.\n");
        }
        for(i=0; i<nts; i++){
            for(l=0; l<L; l++){
                Ml = BF.Ml[l];
                parcel_len = data.parcel_len[l];
                for(j=0; j<parcel_len; j++){
                    fprintf(foutts, "%f ", PostSamp.outcome_test_mean[l][i*parcel_len+j]/niter);
                }
            }
            fprintf(foutts, "\n");
        }
        fprintf(foutts, "\n");
        fclose(foutts);
        free(outts);

        char *outts2 = (char *)calloc(500,sizeof(char));
        strcat(outts2, outpath);
        strcat(outts2, "PostMean2_Out_test.txt");
        FILE * foutts2 = fopen(outts2, "wb");
        if(foutts2 == NULL){
            Rprintf("Cannot open file for PostMean2_Out_test.\n");
        }
        for(i=0; i<nts; i++){
            for(l=0; l<L; l++){
                Ml = BF.Ml[l];
                parcel_len = data.parcel_len[l];
                for(j=0; j<parcel_len; j++){
                    fprintf(foutts2, "%f ", PostSamp.outcome_test_mean2[l][i*parcel_len+j]/niter);
                }
            }
            fprintf(foutts2, "\n");
        }
        fprintf(foutts2, "\n");
        fclose(foutts2);
        free(outts2);
        

    }
    
    // outcome_train
    char *outtr = (char *)calloc(500,sizeof(char));
    strcat(outtr, outpath);
    strcat(outtr, "PostMean_Out_train.txt");
    FILE * fouttr = fopen(outtr, "wb");
    if(fouttr == NULL){
        Rprintf("Cannot open file for PostMean_Out_train.\n");
    }
    for(i=0; i<nobs; i++){
        for(l=0; l<L; l++){
            Ml = BF.Ml[l];
            parcel_len = data.parcel_len[l];
            for(j=0; j<parcel_len; j++){
                fprintf(fouttr, "%f ", PostSamp.outcome_train_mean[l][i*parcel_len+j]/niter);
            }
        }
        fprintf(fouttr, "\n");
    }
    fprintf(fouttr, "\n");
    fclose(fouttr);
    free(outtr);
    
    char *outtr2 = (char *)calloc(500,sizeof(char));
    strcat(outtr2, outpath);
    strcat(outtr2, "PostMean2_Out_train.txt");
    FILE * fouttr2 = fopen(outtr2, "wb");
    if(fouttr2 == NULL){
        Rprintf("Cannot open file for PostMean2_Out_train.\n");
    }
    for(i=0; i<nobs; i++){
        for(l=0; l<L; l++){
            Ml = BF.Ml[l];
            parcel_len = data.parcel_len[l];
            for(j=0; j<parcel_len; j++){
                fprintf(fouttr2, "%f ", PostSamp.outcome_train_mean2[l][i*parcel_len+j]/niter);
            }
        }
        fprintf(fouttr2, "\n");
    }
    fprintf(fouttr2, "\n");
    fclose(fouttr2);
    free(outtr2);
    
    // Intercept_test
    if(nts > 0){
        char *ozinter = (char *)calloc(500,sizeof(char));
        strcat(ozinter, outpath);
        strcat(ozinter, "PostMean_Zinter_test.txt");
        FILE * fzinter = fopen(ozinter, "wb");
        if(fzinter == NULL){
            Rprintf("Cannot open file for PostMean_Zinter_test.\n");
        }
        for(l=0; l<L; l++){
            parcel_len = data.parcel_len[l];
            for(j=0; j<parcel_len; j++){
                fprintf(fzinter, "%f ", PostSamp.Zinter_test_mean[l][j]/niter);
            }
            fprintf(fzinter, "\n");
        }
        fprintf(fzinter, "\n");
        fclose(fzinter);
        free(ozinter);
    }
    
    
    // Intercept_train
    char *ozintertr = (char *)calloc(500,sizeof(char));
    strcat(ozintertr, outpath);
    strcat(ozintertr, "PostMean_Zinter.txt");
    FILE * fzintertr = fopen(ozintertr, "wb");
    if(fzintertr == NULL){
        Rprintf("Cannot open file for PostMean_Zinter.\n");
    }
    for(l=0; l<L; l++){
        parcel_len = data.parcel_len[l];
        for(j=0; j<parcel_len; j++){
            fprintf(fzintertr, "%f ", PostSamp.mean_Zinter[l][j]/niter);
        }
        fprintf(fzintertr, "\n");
    }
    fprintf(fzintertr, "\n");
    fclose(fzintertr);
    free(ozintertr);
    
    
    // Theta
    char *mptheta = (char *)calloc(500,sizeof(char));
    strcat(mptheta, outpath);
    strcat(mptheta, "PostMean_Theta.txt");
    FILE * mftheta = fopen(mptheta, "wb");
    if(mftheta == NULL){
        Rprintf("Cannot open file for PostMean_Theta.\n");
    }
    for(i=0; i<nobs; i++){
        for(l=0; l<L; l++){
            Ml = BF.Ml[l];
            for(j=0; j<Ml; j++){
                fprintf(mftheta, "%f ", PostSamp.mean_theta[l][i*Ml+j]/niter);
            }
        }
        fprintf(mftheta, "\n");
    }
    fprintf(mftheta, "\n");
    fclose(mftheta);
    free(mptheta);
    
    
    // Theta_test
    if(nts > 0){
        char *mpthetats = (char *)calloc(500,sizeof(char));
        strcat(mpthetats, outpath);
        strcat(mpthetats, "PostMean_Theta_test.txt");
        FILE * mfthetats = fopen(mpthetats, "wb");
        if(mfthetats == NULL){
            Rprintf("Cannot open file for PostMean_Theta_test.\n");
        }
        for(i=0; i<nts; i++){
            for(l=0; l<L; l++){
                Ml = BF.Ml[l];
                for(j=0; j<Ml; j++){
                    fprintf(mfthetats, "%f ", PostSamp.theta_test_mean[l][i*Ml+j]/niter);
                }
            }
            fprintf(mfthetats, "\n");
        }
        fprintf(mfthetats, "\n");
        fclose(mfthetats);
        free(mpthetats);
    }
    
    
    // Theta_train
    char *mpthetatr = (char *)calloc(500,sizeof(char));
    strcat(mpthetatr, outpath);
    strcat(mpthetatr, "PostMean_Theta_train.txt");
    FILE * mfthetatr = fopen(mpthetatr, "wb");
    if(mfthetatr == NULL){
        Rprintf("Cannot open file for PostMean_Theta_train.\n");
    }
    for(i=0; i<nobs; i++){
        for(l=0; l<L; l++){
            Ml = BF.Ml[l];
            for(j=0; j<Ml; j++){
                fprintf(mfthetatr, "%f ", PostSamp.theta_train_mean[l][i*Ml+j]/niter);
            }
        }
        fprintf(mfthetatr, "\n");
    }
    fprintf(mfthetatr, "\n");
    fclose(mfthetatr);
    free(mpthetatr);
    
    // Loading matrix
    char *mpload = (char *)calloc(500,sizeof(char));
    strcat(mpload, outpath);
    strcat(mpload, "PostMean_Loading.txt");
    FILE * mfload = fopen(mpload, "wb");
    if(mfload == NULL){
        Rprintf("Cannot open file for PostMean_Load.\n");
    }
    for(l=0; l<L; l++){
        Ml = BF.Ml[l];
        for(i=0; i<Ml; i++){
            for(j=0; j<K; j++){
                fprintf(mfload, "%f ", PostSamp.mean_load[l][i*K+j]/niter);
            }
            fprintf(mfload, "\n");
        }
    }
    fprintf(mfload, "\n");
    fclose(mfload);
    free(mpload);

    // LoadingSTAR matrix
    char *mploads = (char *)calloc(500,sizeof(char));
    strcat(mploads, outpath);
    strcat(mploads, "PostMean_LoadingStar.txt");
    FILE * mfloads = fopen(mploads, "wb");
    if(mfloads == NULL){
        Rprintf("Cannot open file for PostMean_LoadStar.\n");
    }
    for(l=0; l<L; l++){
        Ml = BF.Ml[l];
        for(i=0; i<Ml; i++){
            for(j=0; j<K; j++){
                fprintf(mfloads, "%f ", PostSamp.mean_loadstar[l][i*K+j]/niter);
            }
            fprintf(mfloads, "\n");
        }
    }
    fprintf(mfloads, "\n");
    fclose(mfloads);
    free(mploads);
    
    // Latent
    char *mlatent = (char *)calloc(500,sizeof(char));
    strcat(mlatent, outpath);
    strcat(mlatent, "PostMean_Latent.txt");
    FILE * mflatent = fopen(mlatent, "wb");
    if(mflatent == NULL){
        Rprintf("Cannot open file for PostMean_Latent.\n");
    }
    for(i=0; i<nobs; i++){
        for(l=0; l<L; l++){
            Ml = BF.Ml[l];
            for(k=0; k<K; k++){
                fprintf(mflatent, "%f ", PostSamp.mean_latent[l][i*K+k]/niter);
            }
        }
        fprintf(mflatent, "\n");
    }
    fprintf(mflatent, "\n");
    fclose(mflatent);
    free(mlatent);
    
    // LatentStar
    char *mlatents = (char *)calloc(500,sizeof(char));
    strcat(mlatents, outpath);
    strcat(mlatents, "PostMean_LatentStar.txt");
    FILE * mflatents = fopen(mlatents, "wb");
    if(mflatents == NULL){
        Rprintf("Cannot open file for PostMean_LatentStar.\n");
    }
    
    for(i=0; i<nobs; i++){
        for(l=0; l<L; l++){
            Ml = BF.Ml[l];
            for(k=0; k<K; k++){
                fprintf(mflatents, "%f ", PostSamp.mean_latentstar[l][i*K+k]/niter);
            }
        }
        fprintf(mflatents, "\n");
    }
    fprintf(mflatents, "\n");
    fclose(mflatents);
    free(mlatents);
    
    // Latent_test
    if(nts > 0){
        char *mlatentts = (char *)calloc(500,sizeof(char));
        strcat(mlatentts, outpath);
        strcat(mlatentts, "PostMean_Latent_test.txt");
        FILE * mflatentts = fopen(mlatentts, "wb");
        if(mflatentts == NULL){
            Rprintf("Cannot open file for PostMean_Latent_test.\n");
        }
        for(i=0; i<nts; i++){
            for(l=0; l<L; l++){
                Ml = BF.Ml[l];
                for(k=0; k<K; k++){
                    fprintf(mflatentts, "%f ", PostSamp.latent_test_mean[l][i*K+k]/niter);
                }
            }
            fprintf(mflatentts, "\n");
        }
        fprintf(mflatentts, "\n");
        fclose(mflatentts);
        free(mlatentts);
    }
    
    
    // Latent_train
    char *mlatenttr = (char *)calloc(500,sizeof(char));
    strcat(mlatenttr, outpath);
    strcat(mlatenttr, "PostMean_Latent_train.txt");
    FILE * mflatenttr = fopen(mlatenttr, "wb");
    if(mflatenttr == NULL){
        Rprintf("Cannot open file for PostMean_Latent_train.\n");
    }
    for(i=0; i<nobs; i++){
        for(l=0; l<L; l++){
            Ml = BF.Ml[l];
            for(k=0; k<K; k++){
                fprintf(mflatenttr, "%f ", PostSamp.latent_train_mean[l][i*K+k]/niter);
            }
        }
        fprintf(mflatenttr, "\n");
    }
    fprintf(mflatenttr, "\n");
    fclose(mflatenttr);
    free(mlatenttr);
    
    
    // Phi2Inv
    char *mphi = (char *)calloc(500,sizeof(char));
    strcat(mphi, outpath);
    strcat(mphi, "PostMean_Phi.txt");
    FILE * mfphi = fopen(mphi, "wb");
    if(mfphi == NULL){
        Rprintf("Cannot open file for PostMean_Phi.\n");
    }
    for(j=0; j<K*L; j++){
        fprintf(mfphi, "%f ", PostSamp.mean_Phi2Inv[j]/niter);
    }
    fprintf(mfphi, "\n");
    fclose(mfphi);
    free(mphi);
    
    // Beta
    char *mbetal = (char *)calloc(500,sizeof(char));
    strcat(mbetal, outpath);
    strcat(mbetal, "PostMean_Beta.txt");
    FILE * mfbetal = fopen(mbetal, "wb");
    if(mfbetal == NULL){
        Rprintf("Cannot open file for PostMean_Beta.\n");
    }
    for(l=0; l<L; l++){
        Ml = BF.Ml[l];
        parcel_len = data.parcel_len[l];
        for(i=0; i<parcel_len; i++){
            for(j=0; j<K; j++){
                fprintf(mfbetal, "%f ", PostSamp.mean_beta[l][i*K+j]/niter);
            }
            fprintf(mfbetal, "\n");
        }
    }
    fprintf(mfbetal, "\n");
    fclose(mfbetal);
    free(mbetal);
    
    // BetaStar
    char *mbetals = (char *)calloc(500,sizeof(char));
    strcat(mbetals, outpath);
    strcat(mbetals, "PostMean_BetaStar.txt");
    FILE * mfbetals = fopen(mbetals, "wb");
    if(mfbetals == NULL){
        Rprintf("Cannot open file for PostMean_BetaStar.\n");
    }
    for(l=0; l<L; l++){
        Ml = BF.Ml[l];
        parcel_len = data.parcel_len[l];
        for(i=0; i<parcel_len; i++){
            for(j=0; j<K; j++){
                fprintf(mfbetals, "%f ", PostSamp.mean_betastar[l][i*K+j]/niter);
            }
            fprintf(mfbetals, "\n");
        }
    }
    fprintf(mfbetals, "\n");
    fclose(mfbetals);
    free(mbetals);
    
    
    
    // alpha values
    char *pa = (char *)calloc(500,sizeof(char));
    strcat(pa, outpath);
    strcat(pa, "PostMean_Alpha.txt");
    FILE * fa = fopen(pa, "wb");
    if(fa == NULL){
        Rprintf("Cannot open file for PostMean_Alpha.\n");
    }
    for(i=0; i<K; i++){
        for(l=0; l<L; l++){
            Ml = BF.Ml[l];
            for(k=0; k<Ml; k++){
                fprintf(fa, "%f ", PostSamp.mean_alpha[l][i*Ml+k]/niter);
            }
        }
        fprintf(fa, "\n");
    }
    fclose(fa);
    free(pa);
    
    // AlphaStar
    char *pas = (char *)calloc(500,sizeof(char));
    strcat(pas, outpath);
    strcat(pas, "PostMean_AlphaStar.txt");
    FILE * fas = fopen(pas, "wb");
    if(fas == NULL){
        Rprintf("Cannot open file for PostMean_AlphaStar.\n");
    }
    for(i=0; i<K; i++){
        for(l=0; l<L; l++){
            Ml = BF.Ml[l];
            for(k=0; k<Ml; k++){
                fprintf(fas, "%f ", PostSamp.mean_alphastar[l][i*Ml+k]/niter);
            }
        }
        fprintf(fas, "\n");
    }
    fclose(fas);
    free(pas);
    
    // Gamma
    char *mgamma = (char *)calloc(500,sizeof(char));
    strcat(mgamma, outpath);
    strcat(mgamma, "PostMean_Gamma.txt");
    FILE * mfgamma = fopen(mgamma, "wb");
    if(mfgamma == NULL){
        Rprintf("Cannot open file for PostMean_Gamma.\n");
    }
    for(j=0; j<P; j++){
        for(l=0; l<L; l++){
            fprintf(mfgamma, "%f ", PostSamp.mean_gamma[j*L+l]/niter);
        }
        fprintf(mfgamma, "\n");
    }
    fclose(mfgamma);
    free(mgamma);
    
    // Weight
    char *mw = (char *)calloc(500,sizeof(char));
    strcat(mw, outpath);
    strcat(mw, "PostMean_Weight.txt");
    FILE * mfw = fopen(mw, "wb");
    if(mfw == NULL){
        Rprintf("Cannot open file for PostMean_Weight.\n");
    }
    for(l=0; l<L; l++){
        fprintf(mfw, "%f ", PostSamp.mean_w[l]/niter);
    }
    fprintf(mfw, "\n");
    fclose(mfw);
    free(mw);
    
    // Mustar
    char *pmu = (char *)calloc(500,sizeof(char));
    strcat(pmu, outpath);
    strcat(pmu, "PostMean_MuStar.txt");
    FILE * fmu = fopen(pmu, "wb");
    if(fmu == NULL){
        Rprintf("Cannot open file for PostMean_MuStar.\n");
    }
    for(i=0; i<nobs; i++){
        for(l=0; l<L; l++){
            Ml = BF.Ml[l];
            for(k=0; k<K; k++){
                fprintf(fmu, "%f ", PostSamp.mean_mu[l][i*K+k]/niter);
            }
        }
        fprintf(fmu, "\n");
    }
    fprintf(fmu, "\n");
    fclose(fmu);
    free(pmu);

    // Errors: err2inv_e
    // Err2inv_e
    char * pe = (char*)calloc(500, sizeof(char));
    strcat(pe, outpath);
    strcat(pe, "PostMean_Err2inv_e.txt");
    FILE * fe = fopen(pe, "wb");
    if(fe == NULL){
        Rprintf("Cannot open file for PostMean_Err2inv_e.\n");
    }
    fprintf(fe, "%f\n", PostSamp.err2inv_e_mean[0]/niter);
    fprintf(fe, "\n");
    fclose(fe);
    free(pe);
    
    //err2inv_zeta
    char * pezeta = (char*)calloc(500, sizeof(char));
    strcat(pezeta, outpath);
    strcat(pezeta, "PostMean_Err2inv_zeta.txt");
    FILE * fezeta = fopen(pezeta, "wb");
    if(fezeta == NULL){
        Rprintf("Cannot open file for PostMean-Err2inv_zeta.\n");
    }
    fprintf(fezeta, "%f\n", PostSamp.err2inv_zeta_mean[0]/niter);
    fprintf(fezeta, "\n");
    fclose(fezeta);
    free(pezeta);
    
    //err2inv_eps
    char * peps = (char*)calloc(500, sizeof(char));
    strcat(peps, outpath);
    strcat(peps, "PostMean_Err2inv_eps.txt");
    FILE * feps = fopen(peps, "wb");
    if(feps == NULL){
        Rprintf("Cannot open file for PostMean-Err2inv_eps.\n");
    }
    fprintf(feps, "%f\n", PostSamp.err2inv_eps_mean[0]/niter);
    fclose(feps);
    free(peps);

}
