void out_t0(struct Sampling PostSamp, const struct Inputdata data, const struct BasisFunc BF,
            const int L, const int nobs, const int nts,
            const int K, const int P,
            gsl_rng *r, const char *outpath){
    gsl_set_error_handler_off();
    //////////////////////////////////////////////////////////
    //////////////////////// Output //////////////////////////
    //////////////////////////////////////////////////////////
    
    int i, m, k, p, j, l, Ml, parcel_len;

    // Output Initial Training Outcome
    char * pztr = (char *)calloc(500,sizeof(char));
    strcat(pztr, outpath);
    strcat(pztr, "PostOut_train.txt");
    FILE * fztr = fopen(pztr, "wb");
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
    
    // Output Initial Testing Outcome
    if(nts > 0){
        char * pzts = (char *)calloc(500,sizeof(char));
        strcat(pzts, outpath);
        strcat(pzts, "PostOut_test.txt");
        FILE * fzts = fopen(pzts, "wb");
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
    FILE * fzInter = fopen(pzInter, "wb");
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
    
    // Output Initial Theta
    char *ptheta = (char *)calloc(500,sizeof(char));
    strcat(ptheta, outpath);
    strcat(ptheta, "PostTheta.txt");
    FILE * ftheta = fopen(ptheta, "wb");
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
    
    // Initial theta related error terms
    char *pthetaerr = (char *)calloc(500,sizeof(char));
    strcat(pthetaerr, outpath);
    strcat(pthetaerr, "PostThetaErr.txt");
    FILE * fthetaerr = fopen(pthetaerr, "wb");
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
    
    // Output Initial Theta _ testing
    if(nts > 0){
        char *pthetats = (char *)calloc(500,sizeof(char));
        strcat(pthetats, outpath);
        strcat(pthetats, "PostTheta_test.txt");
        FILE * fthetats = fopen(pthetats, "wb");
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
    
    
    
    // Output Initial Loading Star
    char *ploads = (char *)calloc(500,sizeof(char));
    strcat(ploads, outpath);
    strcat(ploads, "PostLoadingStar.txt");
    FILE * floads = fopen(ploads, "wb");
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
 
    // Output Initial Loading
    char *pload = (char *)calloc(500,sizeof(char));
    strcat(pload, outpath);
    strcat(pload , "PostLoading.txt");
    FILE * fload = fopen(pload, "wb");
    if(fload == NULL){
        Rprintf("Cannot open file for PostLoading.\n");
    }
    for(l=0; l<L; l++){
        Ml = BF.Ml[l];
        for(m=0; m<Ml; m++){
            for(k=0; k<K; k++){
                fprintf(fload, "%f ", PostSamp.load[l][m*K+k]);
            }
        }
    }
    fprintf(fload, "\n");
    fclose(fload);
    free(pload);
    
    
    
    // Initial Latent values
    char *platent = (char *)calloc(500,sizeof(char));
    strcat(platent, outpath);
    strcat(platent, "PostLatent.txt");
    FILE * flatent = fopen(platent, "wb");
    if(flatent == NULL){
        Rprintf("Cannot open file for PostLatent.\n");
    }
    for(i=0; i<nobs; i++){
        for(l=0; l<L; l++){
            for(k=0; k<K; k++){
                fprintf(flatent, "%f ", PostSamp.latent[l][i*K+k]);
            }
        }
    }
    fprintf(flatent, "\n");
    fclose(flatent);
    free(platent);

    

    
    // Initial LatentStar values
    char *platents = (char *)calloc(500,sizeof(char));
    strcat(platents, outpath);
    strcat(platents, "PostLatentStar.txt");
    FILE * flatents = fopen(platents, "wb");
    if(flatents == NULL){
        Rprintf("Cannot open file for PostLatentStar.\n");
    }
    for(i=0; i<nobs; i++){
        for(l=0; l<L; l++){
            for(k=0; k<K; k++){
                fprintf(flatents, "%f ", PostSamp.latentstar[l][i*K+k]);
            }
        }
    }
    fprintf(flatents, "\n");
    fclose(flatents);
    free(platents);
    
    // Initial latent variable related error terms
    char *platenterr = (char *)calloc(500,sizeof(char));
    strcat(platenterr, outpath);
    strcat(platenterr, "PostLatentErr.txt");
    FILE * flatenterr = fopen(platenterr, "wb");
    if(flatenterr == NULL){
        Rprintf("Cannot open file for PostLatentErr.\n");
    }
    for(i=0; i<nobs; i++){
        for(l=0; l<L; l++){
            for(k=0; k<K; k++){
                fprintf(flatenterr, "%f ", PostSamp.latent_train_err[l][i*K+k]);
            }
        }
    }
    fprintf(flatenterr, "\n");
    fclose(flatenterr);
    free(platenterr);
    
    
    // Initial testing Latent values
    if(nts > 0){
        char *platentts = (char *)calloc(500,sizeof(char));
        strcat(platentts, outpath);
        strcat(platentts, "PostLatent_test.txt");
        FILE * flatentts = fopen(platentts, "wb");
        if(flatentts == NULL){
            Rprintf("Cannot open file for PostLatent_test.\n");
        }
        for(i=0; i<nts; i++){
            for(l=0; l<L; l++){
                for(k=0; k<K; k++){
                    fprintf(flatentts, "%f ", PostSamp.latent_test[l][i*K+k]);
                }
            }
        }
        fprintf(flatentts, "\n");
        fclose(flatentts);
        free(platentts);
        
    }
    
    
    
    // Initial Phi2Inv
    char *pphi = (char *)calloc(500,sizeof(char));
    strcat(pphi, outpath);
    strcat(pphi, "PostPhi2Inv.txt");
    FILE * fphi = fopen(pphi, "wb");
    if(fphi == NULL){
        Rprintf("Cannot open file for PostPhi2Inv.\n");
    }
    for(k=0; k<K*L; k++){
        fprintf(fphi, "%f ", PostSamp.Phi2Inv[k]);
    }
    fprintf(fphi, "\n");
    fclose(fphi);
    free(pphi);
    
    // Initial beta values
    char *pbeta = (char *)calloc(500,sizeof(char));
    strcat(pbeta, outpath);
    strcat(pbeta, "PostBeta.txt");
    FILE * fbeta = fopen(pbeta, "wb");
    if(fbeta == NULL){
        Rprintf("Cannot open file for PostBeta.\n");
    }
    for(l=0; l<L; l++){
        parcel_len = data.parcel_len[l];
        for(i=0; i<parcel_len; i++){
            for(k=0; k<K; k++){
                fprintf(fbeta, "%f ", PostSamp.beta[l][i*K+k]);
            }
        }
    }
    fprintf(fbeta, "\n");
    fclose(fbeta);
    free(pbeta);
    
    // BetaStar
    char *pbetas = (char *)calloc(500,sizeof(char));
    strcat(pbetas, outpath);
    strcat(pbetas, "PostBetaStar.txt");
    FILE * fbetas = fopen(pbetas, "wb");
    if(fbetas == NULL){
        Rprintf("Cannot open file for PostBetaStar.\n");
    }
    for(l=0; l<L; l++){
        parcel_len = data.parcel_len[l];
        for(i=0; i<parcel_len; i++){
            for(k=0; k<K; k++){
                fprintf(fbetas, "%f ", PostSamp.betastar[l][i*K+k]);
            }
        }
    }
    fprintf(fbetas, "\n");
    fclose(fbetas);
    free(pbetas);
    
    // Initial alpha values
    char *pa = (char *)calloc(500,sizeof(char));
    strcat(pa, outpath);
    strcat(pa, "PostAlpha.txt");
    FILE * fa = fopen(pa, "wb");
    if(fa == NULL){
        Rprintf("Cannot open file for PostAlpha.\n");
    }
    for(i=0; i<K; i++){
        for(l=0; l<L; l++){
            Ml = BF.Ml[l];
            for(k=0; k<Ml; k++){
                fprintf(fa, "%f ", PostSamp.alpha[l][i*Ml+k]);
            }
        }
    }
    fprintf(fa, "\n");
    fclose(fa);
    free(pa);
    
    // AlphaStar
    char *pas = (char *)calloc(500,sizeof(char));
    strcat(pas, outpath);
    strcat(pas, "PostAlphaStar.txt");
    FILE * fas = fopen(pas, "wb");
    if(fas == NULL){
        Rprintf("Cannot open file for PostAlphaStar.\n");
    }
    for(i=0; i<K; i++){
        for(l=0; l<L; l++){
            Ml = BF.Ml[l];
            for(k=0; k<Ml; k++){
                fprintf(fas, "%f ", PostSamp.alphastar[l][i*Ml+k]);
            }
        }
    }
    fprintf(fas, "\n");
    fclose(fas);
    free(pas);
    
    
    // mustar
    char *pmu = (char *)calloc(500,sizeof(char));
    strcat(pmu, outpath);
    strcat(pmu, "PostMu.txt");
    FILE * fmu = fopen(pmu, "wb");
    if(fmu == NULL){
        Rprintf("Cannot open file for PostMuStar.\n");
    }
    for(i=0; i<nobs; i++){
        for(l=0; l<L; l++){
            for(k=0; k<K; k++){
                fprintf(fmu, "%f ", PostSamp.mustar[l][i*K+k]);
            }
        }
    }
    fprintf(fmu, "\n");
    fclose(fmu);
    free(pmu);
    
    // Initial Gamma
    char *pgamma = (char *)calloc(500,sizeof(char));
    strcat(pgamma, outpath);
    strcat(pgamma, "PostGamma.txt");
    FILE * fgamma = fopen(pgamma, "wb");
    if(fgamma == NULL){
        Rprintf("Cannot open file for PostGamma.\n");
    }
    for(p=0; p<P*L; p++){
        fprintf(fgamma, "%f ", PostSamp.gamma[p]);
    }
    fprintf(fgamma, "\n");
    fclose(fgamma);
    free(pgamma);
    
    // Initial Weight
    char *pw = (char *)calloc(500,sizeof(char));
    strcat(pw, outpath);
    strcat(pw, "PostWeight.txt");
    FILE * fw = fopen(pw, "wb");
    if(fw == NULL){
        Rprintf("Cannot open file for PostWeight.\n");
    }
    for(p=0; p<L; p++){
        fprintf(fw, "%f ", PostSamp.w[p]);
    }
    fprintf(fw, "\n");
    fclose(fw);
    free(pw);
    
    //Initial err2inv_e
    char * pe = (char*)calloc(500, sizeof(char));
    strcat(pe, outpath);
    strcat(pe, "PostErr2inv_e.txt");
    FILE * fe = fopen(pe, "wb");
    if(fe == NULL){
        Rprintf("Cannot open file for PostErr2inv_e.\n");
    }
    fprintf(fe, "%f\n", PostSamp.err2inv_e[0]);
    fprintf(fe, "\n");
    fclose(fe);
    free(pe);
    
    //Initial err2inv_zeta
    char * pz = (char*)calloc(500, sizeof(char));
    strcat(pz, outpath);
    strcat(pz, "PostErr2inv_zeta.txt");
    FILE * fz = fopen(pz, "wb");
    if(fz == NULL){
        Rprintf("Cannot open file for PostErr2inv_zeta.\n");
    }
    fprintf(fz, "%f\n", PostSamp.err2inv_zeta[0]);
    fprintf(fz, "\n");
    fclose(fz);
    free(pz);
    
    //Initial err2inv_eps
    char * peps = (char*)calloc(500, sizeof(char));
    strcat(peps, outpath);
    strcat(peps, "PostErr2inv_eps.txt");
    FILE * feps = fopen(peps, "wb");
    if(feps == NULL){
        Rprintf("Cannot open file for PostErr2inv_eps.\n");
    }
    fprintf(feps, "%f ", PostSamp.err2inv_eps[0]);
    fprintf(feps, "\n");
    fclose(feps);
    free(peps);

}


void print_t0(struct Sampling PostSamp, const struct Inputdata data, const struct BasisFunc BF,
              const int P, const int K, const int L){
    int p, k, m, l, Ml, parcel_len;
    Ml = BF.Ml[0];
    parcel_len = data.parcel_len[0];
    
    // Print
    Rprintf("\nInitial Gamma\n");
    for(p=0; p<3; p++){
        for(l=0; l<L; l++){
            Rprintf("%.6f ", PostSamp.gamma[p*L+l]);
        }
        Rprintf("\n");
    }
    
    Rprintf("\nInitial Alpha\n");
    for(k=0; k<K; k++){
        for(m=0; m<5; m++){
            Rprintf("%.6f ", PostSamp.alpha[0][k*Ml+m]);
        }
        Rprintf("\n");
    }
    
    Rprintf("\nInitial Beta\n");
    for(k=0; k<5; k++){
        for(m=0; m<K; m++){
            Rprintf("%.6f ", PostSamp.beta[0][k*K+m]);
        }
        Rprintf("\n");
    }
    
    Rprintf("\nInitial BetaStar\n");
    for(k=0; k<5; k++){
        for(m=0; m<K; m++){
            Rprintf("%.6f ", PostSamp.betastar[0][k*K+m]);
        }
        Rprintf("\n");
    }
    
    
    Rprintf("\nInitial SumX\n");
    for(k=0; k<5; k++){
        for(m=0; m<5; m++){
            Rprintf("%.6f ", PostSamp.sumX[0][k*parcel_len+m]);
        }
        Rprintf("\n");
    }
    
    Rprintf("\nInitial X\n");
    for(k=0; k<5; k++){
        for(m=0; m<5; m++){
            Rprintf("%.6f ", data.Xl[0][k*parcel_len*P+m]);
        }
        Rprintf("\n");
    }
    
    Rprintf("\nInitial Latent\n");
    for(k=0; k<5; k++){
        for(m=0; m<K; m++){
            Rprintf("%.6f ", PostSamp.latentstar[0][k*K+m]);
        }
        Rprintf("\n");
    }

    Rprintf("\nInitial Mustar\n");
    for(k=0; k<5; k++){
        for(m=0; m<K; m++){
            Rprintf("%.6f ", PostSamp.mustar[0][k*K+m]);
        }
        Rprintf("\n");
    }
    
    
    Rprintf("\nInitial Load\n");
    for(k=0; k<5; k++){
        for(m=0; m<K; m++){
            Rprintf("%.6f ", PostSamp.load[0][k*K+m]);
        }
        Rprintf("\n");
    }
    
    Rprintf("\nInitial Theta\n");
    for(k=0; k<5; k++){
        for(m=0; m<5; m++){
            Rprintf("%.6f ", PostSamp.theta[0][k*Ml+m]);
        }
        Rprintf("\n");
    }
    
    Rprintf("\nInitial Phi\n");
    for(k=0; k<K; k++){
        for(l=0; l<L; l++){
            Rprintf("%.6f ", PostSamp.Phi2Inv[k*L+l]);
        }
        Rprintf("\n");
    }
    
    // Initial sigma
    Rprintf("\nInitial sigma^2\n");
    Rprintf("e=%.6f, zeta=%.6f, eps=%.6f\n",
           PostSamp.err2inv_e[0],PostSamp.err2inv_zeta[0], PostSamp.err2inv_eps[0]);

}
