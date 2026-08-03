void output_samp(const char * outpath, const int t, const int burnin, const int nobs, const int nts, const int M, const int K, const int L, const int P, const int *singular, struct Sampling PostSamp, const struct Inputdata data, const struct BasisFunc BF){
    
    int i, j, k, m, l, Ml, parcel_len;
    
    // Training Outcome
    char * pztr = (char *)calloc(500,sizeof(char));
    strcat(pztr, outpath);
    strcat(pztr, "PostOut_train.txt");
    FILE * fztr = fopen(pztr, "a+");
    if(fztr == NULL){
        Rprintf("Cannot open file for PostOut_train.\n");
    }
    for(i=0; i<nobs; i++){
        
        for(l=0; l<L; l++){
            parcel_len = data.parcel_len[l];
            for(j=0; j<parcel_len; j++){
                fprintf(fztr, "%f ", PostSamp.outcome_train[l][i*parcel_len+j]);
            }
        }
    }
    fprintf(fztr, "\n");
    fclose(fztr);
    free(pztr);
    
    // Test Outcome
    if(nts > 0){
        char * pzts = (char *)calloc(500,sizeof(char));
        strcat(pzts, outpath);
        strcat(pzts, "PostOut_test.txt");
        FILE * fzts = fopen(pzts, "a+");
        if(fzts == NULL){
            Rprintf("Cannot open file for PostOut_test.\n");
        }
        for(i=0; i<nts; i++){
            for(l=0; l<L; l++){
                parcel_len = data.parcel_len[l];
                for(j=0; j<parcel_len; j++){
                    fprintf(fzts, "%f ", PostSamp.outcome_test[l][i*parcel_len+j]);
                }
            }
        }
        fprintf(fzts, "\n");
        fclose(fzts);
        free(pzts);

    }
    
    
    // Intercept term
    char *pzInter = (char *)calloc(500, sizeof(char));
    strcat(pzInter, outpath);
    strcat(pzInter, "PostZInter.txt");
    FILE * fzInter = fopen(pzInter, "a+");
    if(fzInter == NULL){
        Rprintf("Cannot open file for PostZInter.\n");
    }
    for(l=0; l<L; l++) {
        parcel_len = data.parcel_len[l];
        for(i=0; i<parcel_len; i++) {
            fprintf(fzInter, "%f ", PostSamp.Zinter[l][i]);
        }
        fprintf(fzInter, "\n");
    }
    fprintf(fzInter, "\n");
    fclose(fzInter);
    free(pzInter);
    
    
    // theta
    char *ptheta = (char *)calloc(500,sizeof(char));
    strcat(ptheta, outpath);
    strcat(ptheta, "PostTheta.txt");
    FILE * ftheta = fopen(ptheta, "a+");
    if(ftheta == NULL){
        Rprintf("Cannot open file for PostTheta.\n");
    }
    for(i=0; i<nobs; i++){
        for(l=0; l<L; l++){
            Ml = BF.Ml[l];
            for(j=0; j<Ml; j++){
                fprintf(ftheta, "%f ", PostSamp.theta[l][i*Ml+j]);
            }
        }
    }
    fprintf(ftheta, "\n");
    fclose(ftheta);
    free(ptheta);
    
    // error term related to theta
    char *pthetaerr = (char *)calloc(500,sizeof(char));
    strcat(pthetaerr, outpath);
    strcat(pthetaerr, "PostThetaErr.txt");
    FILE * fthetaerr = fopen(pthetaerr, "a+");
    if(fthetaerr == NULL){
        Rprintf("Cannot open file for PostThetaErr.\n");
    }
    for(i=0; i<nobs; i++){
        for(l=0; l<L; l++){
            Ml = BF.Ml[l];
            for(j=0; j<Ml; j++){
                fprintf(fthetaerr, "%f ", PostSamp.theta_train_err[l][i*Ml+j]);
            }
        }
    }
    fprintf(fthetaerr, "\n");
    fclose(fthetaerr);
    free(pthetaerr);
    
    
    // Theta_test
    if(nts > 0){
        char *pthetats = (char *)calloc(500,sizeof(char));
        strcat(pthetats, outpath);
        strcat(pthetats, "PostTheta_test.txt");
        FILE * fthetats = fopen(pthetats, "a+");
        if(fthetats == NULL){
            Rprintf("Cannot open file for PostTheta_test.\n");
        }
        for(i=0; i<nts; i++){
            for(l=0; l<L; l++){
                Ml = BF.Ml[l];
                for(j=0; j<Ml; j++){
                    fprintf(fthetats, "%f ", PostSamp.theta_test[l][i*Ml+j]);
                }
            }
        }
        fprintf(fthetats, "\n");
        fclose(fthetats);
        free(pthetats);
    }
    
    
    // loading matrix
    char *pload = (char *)calloc(500,sizeof(char));
    strcat(pload, outpath);
    strcat(pload, "PostLoading.txt");
    FILE * fload = fopen(pload, "a+");
    if(fload == NULL){
        Rprintf("Cannot open file for PostLoading.\n");
    }
    
    for(l=0; l<L; l++){
        Ml = BF.Ml[l];
        for(i=0; i<Ml; i++){
            for(j=0; j<K; j++){
                fprintf(fload, "%f ", PostSamp.load[l][i*K+j]);
            }
        }
    }
    
    fprintf(fload, "\n");
    fclose(fload);
    free(pload);
    
    // Loading Star
    char *ploads = (char *)calloc(500,sizeof(char));
    strcat(ploads, outpath);
    strcat(ploads , "PostLoadingStar.txt");
    FILE * floads = fopen(ploads, "a+");
    if(floads == NULL){
        Rprintf("Cannot open file for PostLoadingStar.\n");
    }
    for(l=0; l<L; l++){
        Ml = BF.Ml[l];
        for(m=0; m<Ml; m++){
            for(k=0; k<K; k++){
                fprintf(floads, "%f ", PostSamp.loadstar[l][m*K+k]);
            }
        }
    }
    fprintf(floads, "\n");
    fclose(floads);
    free(ploads);
    
    // Latent values
    char *platent = (char *)calloc(500,sizeof(char));
    strcat(platent, outpath);
    strcat(platent, "PostLatent.txt");
    FILE * flatent = fopen(platent, "a+");
    if(flatent == NULL){
        Rprintf("Cannot open file for PostLatent.\n");
    }
    for(i=0; i<nobs; i++){
        for(l=0; l<L; l++){
            Ml = BF.Ml[l];
            for(k=0; k<K; k++){
                fprintf(flatent, "%f ", PostSamp.latent[l][i*K+k]);
            }
        }
    }
    fprintf(flatent, "\n");
    fclose(flatent);
    free(platent);
    
    // LatentStar values
    char *platents = (char *)calloc(500,sizeof(char));
    strcat(platents, outpath);
    strcat(platents, "PostLatentStar.txt");
    FILE * flatents = fopen(platents, "a+");
    if(flatents == NULL){
        Rprintf("Cannot open file for PostLatentStar.\n");
    }
    for(i=0; i<nobs; i++){
        for(l=0; l<L; l++){
            Ml = BF.Ml[l];
            for(k=0; k<K; k++){
                fprintf(flatents, "%f ", PostSamp.latentstar[l][i*K+k]);
            }
        }
    }
    fprintf(flatents, "\n");
    fclose(flatents);
    free(platents);
    
    // Errors related to Latent values
    char *platenterr = (char *)calloc(500,sizeof(char));
    strcat(platenterr, outpath);
    strcat(platenterr, "PostLatentErr.txt");
    FILE * flatenterr = fopen(platenterr, "a+");
    if(flatenterr == NULL){
        Rprintf("Cannot open file for PostLatentErr.\n");
    }
    for(i=0; i<nobs; i++){
        for(l=0; l<L; l++){
            Ml = BF.Ml[l];
            for(k=0; k<K; k++){
                fprintf(flatenterr, "%f ", PostSamp.latent_train_err[l][i*K+k]);
            }
        }
    }
    fprintf(flatenterr, "\n");
    fclose(flatenterr);
    free(platenterr);
    
    // Latent_test
    if(nts > 0){
        char *platentts = (char *)calloc(500,sizeof(char));
        strcat(platentts, outpath);
        strcat(platentts, "PostLatent_test.txt");
        FILE * flatentts = fopen(platentts, "a+");
        if(flatentts == NULL){
            Rprintf("Cannot open file for PostLatent_test.\n");
        }
        for(i=0; i<nts; i++){
            for(l=0; l<L; l++){
                Ml = BF.Ml[l];
                for(k=0; k<K; k++){
                    fprintf(flatentts, "%f ", PostSamp.latent_test[l][i*K+k]);
                }
            }
        }
        fprintf(flatentts, "\n");
        fclose(flatentts);
        free(platentts);
    }
    
    
    
    // Phi
    char *pphi = (char *)calloc(500,sizeof(char));
    strcat(pphi, outpath);
    strcat(pphi, "PostPhi2Inv.txt");
    FILE * fphi = fopen(pphi, "a+");
    if(fphi == NULL){
        Rprintf("Cannot open file for PostPhi2Inv.\n");
    }
    for(k=0; k<K*L; k++){
        fprintf(fphi, "%f ", PostSamp.Phi2Inv[k]);
    }
    fprintf(fphi, "\n");
    fclose(fphi);
    free(pphi);
    

    // Beta
    char *pb = (char *)calloc(500,sizeof(char));
    strcat(pb, outpath);
    strcat(pb, "PostBeta.txt");
    FILE * fb = fopen(pb, "a+");
    if(fb == NULL){
        Rprintf("Cannot open file for PostBeta.\n");
    }
    for(l=0; l<L; l++){
        parcel_len = data.parcel_len[l];
        for(i=0; i<parcel_len; i++){
            for(j=0; j<K; j++){
                fprintf(fb, "%f ", PostSamp.beta[l][i*K+j]);
            }
        }
    }
    fprintf(fb, "\n");
    fclose(fb);
    free(pb);
    
    // BetaSTAR
    char *pbs = (char *)calloc(500,sizeof(char));
    strcat(pbs, outpath);
    strcat(pbs, "PostBetaStar.txt");
    FILE * fbs = fopen(pbs, "a+");
    if(fbs == NULL){
        Rprintf("Cannot open file for PostBetaStar.\n");
    }
    for(l=0; l<L; l++){
        parcel_len = data.parcel_len[l];
        for(i=0; i<parcel_len; i++){
            for(j=0; j<K; j++){
                fprintf(fbs, "%f ", PostSamp.betastar[l][i*K+j]);
            }
        }
    }
    fprintf(fbs, "\n");
    fclose(fbs);
    free(pbs);
    
    // alpha values
    char *pa = (char *)calloc(500,sizeof(char));
    strcat(pa, outpath);
    strcat(pa, "PostAlpha.txt");
    FILE * fa = fopen(pa, "a+");
    if(fa == NULL){
        Rprintf("Cannot open file for PostAlpha.\n");
    }
    for(i=0; i<K; i++){
        for(l=0; l<L;l++){
            Ml = BF.Ml[l];
            for(k=0; k<Ml; k++){
                fprintf(fa, "%f ", PostSamp.alpha[l][i*Ml+k]);
            }
        }
    }
    fprintf(fa, "\n");
    fclose(fa);
    free(pa);
    
    // AlphaSTAR
    char *pas = (char *)calloc(500,sizeof(char));
    strcat(pas, outpath);
    strcat(pas, "PostAlphaStar.txt");
    FILE * fas = fopen(pas, "a+");
    if(fas == NULL){
        Rprintf("Cannot open file for PostAlphaStar.\n");
    }
    for(i=0; i<K; i++){
        for(l=0; l<L;l++){
            Ml = BF.Ml[l];
            for(k=0; k<Ml; k++){
                fprintf(fas, "%f ", PostSamp.alphastar[l][i*Ml+k]);
            }
        }
    }
    fprintf(fas, "\n");
    fclose(fas);
    free(pas);
    
    
    // Gamma
    char *pgamma = (char *)calloc(500,sizeof(char));
    strcat(pgamma, outpath);
    strcat(pgamma, "PostGamma.txt");
    FILE * fgamma = fopen(pgamma, "a+");
    if(fgamma == NULL){
        Rprintf("Cannot open file for PostGamma.\n");
    }
    for(k=0; k<P*L; k++){
        fprintf(fgamma, "%f ", PostSamp.gamma[k]);
    }
    fprintf(fgamma, "\n");
    fclose(fgamma);
    free(pgamma);
    
    // Weight
    char *pw = (char *)calloc(500,sizeof(char));
    strcat(pw, outpath);
    strcat(pw, "PostWeight.txt");
    FILE * fw = fopen(pw, "a+");
    if(fw== NULL){
        Rprintf("Cannot open file for PostWeight.\n");
    }
    for(k=0; k<L; k++){
        fprintf(fw, "%f ", PostSamp.w[k]);
    }
    fprintf(fw, "\n");
    fclose(fw);
    free(pw);
    
    // Err2inv_e
    char * pe = (char*)calloc(500, sizeof(char));
    strcat(pe, outpath);
    strcat(pe, "PostErr2inv_e.txt");
    FILE * fe = fopen(pe, "a+");
    if(fe == NULL){
        Rprintf("Cannot open file for PostErr2inv_e.\n");
    }
    fprintf(fe, "%f\n", PostSamp.err2inv_e[0]);
    fprintf(fe, "\n");
    fclose(fe);
    free(pe);
    
    // Err2inv_eps
    char * peps = (char*)calloc(500, sizeof(char));
    strcat(peps, outpath);
    strcat(peps, "PostErr2inv_eps.txt");
    FILE * feps = fopen(peps, "a+");
    if(feps == NULL){
        Rprintf("Cannot open file for PostErr2inv_eps.\n");
    }
    fprintf(feps, "%f\n", PostSamp.err2inv_eps[0]);
    fclose(feps);
    free(peps);
    
    // Err2inv_e
    char * pzeta = (char*)calloc(500, sizeof(char));
    strcat(pzeta, outpath);
    strcat(pzeta, "PostErr2inv_zeta.txt");
    FILE * fzeta = fopen(pzeta, "a+");
    if(fzeta == NULL){
        Rprintf("Cannot open file for PostErr2inv_zeta.\n");
    }
    fprintf(fzeta, "%f\n", PostSamp.err2inv_zeta[0]);
    fprintf(fzeta, "\n");
    fclose(fzeta);
    free(pzeta);
    
    // mustar
    char *pmu = (char *)calloc(500,sizeof(char));
    strcat(pmu, outpath);
    strcat(pmu, "PostMu.txt");
    FILE * fmu = fopen(pmu, "a+");
    if(fmu == NULL){
        Rprintf("Cannot open file for PostMuStar.\n");
    }
    for(i=0; i<nobs; i++){
        
        for(l=0; l<L; l++){
            Ml = BF.Ml[l];
            for(k=0; k<K; k++){
                fprintf(fmu, "%f ", PostSamp.mustar[l][i*K+k]);
            }
        }
    }
    fprintf(fmu, "\n");
    fclose(fmu);
    free(pmu);
    
}

