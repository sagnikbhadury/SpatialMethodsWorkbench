void init_theta_rand(const int nobs, const int L, struct Sampling PostSamp, const struct Inputdata data, const struct BasisFunc BF);
void init_latent_rand(struct Sampling PostSamp, const struct BasisFunc BF, const struct Inputdata data,
                      const int K, const int L, const int nobs, gsl_rng*r);
void init_load_rand(struct Sampling PostSamp, struct BasisFunc BF, const int L, const int K, const int nobs);
void init_phi_rand(struct Sampling PostSamp, const int K,const int L, gsl_rng*r);
void init_alpha_rand(struct Sampling PostSamp, const struct Inputdata data, const struct BasisFunc BF,
                     const int K, const int L, gsl_rng *r);
void init_beta_rand(struct Sampling PostSamp, const struct BasisFunc BF, const struct Inputdata data,
                    const int K, const int L);
void init_w_rand(struct Sampling PostSamp, const int L);
void init_gamma_rand(struct Sampling PostSamp, int P,  int L, gsl_rng *r);
void init_sumx(struct Sampling PostSamp, const struct Inputdata data, const struct BasisFunc BF,
               const int nobs, const int L, const int P );
void init_mu_rand(struct Sampling PostSamp, const int L, const int nobs, const int K, gsl_rng*r);
void init_zInter(const int nobs, const int L, struct Sampling PostSamp,
                 const struct Inputdata data, const struct BasisFunc BF);
void set_initial(const int L, const int nobs, const int nts, const int K, const int P, struct Sampling PostSamp, const struct Inputdata data, const struct BasisFunc BF, gsl_rng *r,  const char *outpath, const bool printInit);

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void set_initial2(const int L, const int nobs, const int nts, const int K, const int P, struct Sampling PostSamp, const struct Inputdata data, const struct BasisFunc BF, gsl_rng *r,  const char *outpath, const bool printInit, const char *inpath){
    gsl_set_error_handler_off();
    int i, j;
    int image_len = 1024;
    int Ml = BF.Ml[0];
    
    // Errors
    PostSamp.err2inv_e[0] = 0.1;//gsl_ran_flat(r, 1.0, 10.);
    PostSamp.err2inv_zeta[0] = 0.2;//gsl_ran_flat(r, 1.0, 10.);
    PostSamp.err2inv_eps[0] =  0.2;//gsl_ran_flat(r, 1.0, 10.);
    PostSamp.err2inv_e_mean[0] = 0.0f;
    PostSamp.err2inv_zeta_mean[0] = 0.0f;
    PostSamp.err2inv_eps_mean[0] = 0.0f;
    int l;
    for(l=0; l<L; l++){
        PostSamp.err2inv_u[l] = 0.0f;
        PostSamp.err2inv_u_mean[l] = 0.0f;
    }
    
    // Theta
    char *stheta = (char *)calloc(500,sizeof(char));
    strcat(stheta, inpath);
    strcat(stheta, "simtheta.txt");
    FILE *rtheta = fopen(stheta, "r");
    if(rtheta == NULL){
        Rprintf("Cannot open file for theta.txt\n");
    }
    for(i=0; i<nobs; i++){
        for(j=0; j<Ml; j++){
            if(fscanf(rtheta, "%f", &PostSamp.theta[0][i*Ml+j]) !=1 ){
                break;
            }
        }
        if (fscanf(rtheta, "\n") != 0){
            break;
        }
    }
    fclose(rtheta);
    free(stheta);
    
    for(i=0; i<nobs; i++){
        for(j=0; j<Ml; j++){
            PostSamp.mean_theta[0][i*Ml+j] = 0.0f;
        }
    }
    
    /*for(i=0; i<5; i++){
        for(j=0; j<5; j++){
            printf("%.3f ", PostSamp.theta[0][i*Ml+j]);
        }
        printf("\n");
    }*/
    
    // phi
    for(i=0; i<K*L; i++){
        PostSamp.Phi2Inv[i] = 1.0;
        PostSamp.mean_Phi2Inv[i] = 0.0f;
    }
    
    // Loading
    char *sload = (char *)calloc(500,sizeof(char));
    strcat(sload, inpath);
    strcat(sload, "simloading.txt");
    FILE *rload = fopen(sload, "r");
    if(rload == NULL){
        Rprintf("Cannot open file for load.txt\n");
    }
    for(i=0; i<Ml; i++){
        for(j=0; j<K; j++){
            if(fscanf(rload, "%f", &PostSamp.load[0][i*K+j]) !=1 ){
                break;
            }
        }
        if (fscanf(rload, "\n") != 0){
            break;
        }
    }
    fclose(rload);
    free(sload);
    
    for(i=0; i<Ml; i++){
        for(j=0; j<K; j++){
            PostSamp.mean_load[0][i*K+j] = 0.0f;
        }
    }
    
    /*for(i=0; i<10; i++){
        for(j=0; j<K; j++){
            printf("%.3f ", PostSamp.load[0][i*K+j]);
        }
        printf("\n");
    }*/
    
    // LoadingStar
    char *sloads = (char *)calloc(500,sizeof(char));
    strcat(sloads, inpath);
    strcat(sloads, "simloadingStar.txt");
    FILE *rloads = fopen(sloads, "r");
    if(rload == NULL){
        Rprintf("Cannot open file for loadStar.txt\n");
    }
    for(i=0; i<Ml; i++){
        for(j=0; j<K; j++){
            if(fscanf(rloads, "%f", &PostSamp.loadstar[0][i*K+j]) !=1 ){
                break;
            }
        }
        if (fscanf(rloads, "\n") != 0){
            break;
        }
    }
    fclose(rloads);
    free(sloads);
    
    for(i=0; i<Ml; i++){
        for(j=0; j<K; j++){
            PostSamp.mean_loadstar[0][i*K+j] = 0.0f;
        }
    }
    

    // latent
    char *slatent = (char *)calloc(500,sizeof(char));
    strcat(slatent, inpath);
    strcat(slatent, "simlatent.txt");
    FILE *rlatent = fopen(slatent, "r");
    if(rlatent == NULL){
        Rprintf("Cannot open file for latent.txt\n");
    }
    for(i=0; i<nobs; i++){
        for(j=0; j<K; j++){
            if(fscanf(rlatent, "%f", &PostSamp.latent[0][i*K+j]) !=1 ){
                break;
            }
        }
        if (fscanf(rlatent, "\n") != 0){
            break;
        }
    }
    fclose(rlatent);
    free(slatent);
    
    for(i=0; i<nobs; i++){
        for(j=0; j<K; j++){
            PostSamp.mean_latent[0][i*K+j] = 0.0f;
        }
    }

    // latentstar
    char *slatents = (char *)calloc(500,sizeof(char));
    strcat(slatents, inpath);
    strcat(slatents, "simlatentStar.txt");
    FILE *rlatents = fopen(slatents, "r");
    if(rlatents == NULL){
        Rprintf("Cannot open file for latentS.txt\n");
    }
    for(i=0; i<nobs; i++){
        for(j=0; j<K; j++){
            if(fscanf(rlatents, "%f", &PostSamp.latentstar[0][i*K+j]) !=1 ){
                break;
            }
        }
        if (fscanf(rlatents, "\n") != 0){
            break;
        }
    }
    fclose(rlatents);
    free(slatents);
    
    for(i=0; i<nobs; i++){
        for(j=0; j<K; j++){
            PostSamp.mean_latentstar[0][i*K+j] = 0.0f;
        }
    }
    
    //mustar
    char *smu = (char *)calloc(500,sizeof(char));
    strcat(smu, inpath);
    strcat(smu, "simmustar.txt");
    FILE *rmu = fopen(smu, "r");
    if(rmu == NULL){
        Rprintf("Cannot open file for mu.txt\n");
    }
    for(i=0; i<nobs; i++){
        for(j=0; j<K; j++){
            if(fscanf(rmu, "%f", &PostSamp.mustar[0][i*K+j]) !=1 ){
                break;
            }
        }
        if (fscanf(rmu, "\n") != 0){
            break;
        }
    }
    fclose(rmu);
    free(smu);
    
    for(i=0; i<nobs; i++){
        for(j=0; j<K; j++){
            PostSamp.mean_mu[0][i*K+j] = 0.0f;
        }
    }
    
    // Alpha
    char *sa = (char *)calloc(500,sizeof(char));
    strcat(sa, inpath);
    strcat(sa, "simalpha.txt");
    FILE *ra = fopen(sa, "r");
    if(ra == NULL){
        Rprintf("Cannot open file for alpha.txt\n");
    }
    for(i=0; i<K; i++){
        for(j=0; j<Ml; j++){
            if(fscanf(ra, "%f", &PostSamp.alpha[0][i*Ml+j]) !=1 ){
                break;
            }
        }
        if (fscanf(ra, "\n") != 0){
            break;
        }
    }
    fclose(ra);
    free(sa);
    
    for(i=0; i<K; i++){
        for(j=0; j<Ml; j++){
            PostSamp.mean_alpha[0][i*Ml+j] = 0.0f;
            PostSamp.mean_alphastar[0][i*Ml+j] = 0.0f;
        }
    }
    
    // Alphastar
    char *sas = (char *)calloc(500,sizeof(char));
    strcat(sas, inpath);
    strcat(sas, "simalphaStar.txt");
    FILE *ras = fopen(sas, "r");
    if(ras == NULL){
        Rprintf("Cannot open file for alphas.txt\n");
    }
    for(i=0; i<K; i++){
        for(j=0; j<Ml; j++){
            if(fscanf(ras, "%f", &PostSamp.alphastar[0][i*Ml+j]) !=1 ){
                break;
            }
        }
        if (fscanf(ras, "\n") != 0){
            break;
        }
    }
    fclose(ras);
    free(sas);
    
    
    // Beta
    char *sb = (char *)calloc(500,sizeof(char));
    strcat(sb, inpath);
    strcat(sb, "simbeta.txt");
    FILE *rb = fopen(sb, "r");
    if(rb == NULL){
        Rprintf("Cannot open file for beta.txt\n");
    }
    for(i=0; i<1024; i++){
        for(j=0; j<K; j++){
            if(fscanf(rb, "%f", &PostSamp.beta[0][i*K+j]) !=1 ){
                break;
            }
        }
        if (fscanf(rb, "\n") != 0){
            break;
        }
    }
    fclose(rb);
    free(sb);

    for(i=0; i<image_len; i++){
        for(j=0; j<K; j++){
            PostSamp.mean_beta[0][i*K+j] = 0.0f;
            PostSamp.mean_betastar[0][i*K+j] = 0.0f;
        }
    }
    
    // betastar
    char *sbs = (char *)calloc(500,sizeof(char));
    strcat(sbs, inpath);
    strcat(sbs, "simbetaStar.txt");
    FILE *rbs = fopen(sbs, "r");
    if(rbs == NULL){
        Rprintf("Cannot open file for betas.txt\n");
    }
    for(i=0; i<1024; i++){
        for(j=0; j<K; j++){
            if(fscanf(rbs, "%f", &PostSamp.betastar[0][i*K+j]) !=1 ){
                break;
            }
        }
        if (fscanf(rbs, "\n") != 0){
            break;
        }
    }
    fclose(rbs);
    free(sbs);
    
    // gamma
    for(i=0; i<20; i++){
        if(i < 5){
            PostSamp.gamma[i] = 1.0;
        }else{
            PostSamp.gamma[i] = 0.0f;
        }
        PostSamp.mean_gamma[i] = 0.0f;
    }

    // w
    PostSamp.w[0]=0.5;
    PostSamp.mean_w[0] = 0.0f;
    
    // summarized_x
    init_sumx(PostSamp, data, BF, nobs, L, P);
    
    // Compute \sum_p gamma*basis*X
    int m, p, k;
    double tempsum;
    for(i=0; i<nobs; i++){
        for(l=0; l<L; l++){
            Ml = BF.Ml[l];
            for(m=0; m<Ml; m++){
                tempsum = 0.0;
                for(p=0; p<P; p++){
                    tempsum += PostSamp.gamma[p*L+l] * PostSamp.xb[l][i*P*Ml + p*Ml+ m];
                }
                PostSamp.rxb[l][i*Ml+m] = tempsum;
            }
        }
    }
    
    
    // Compute xba
    for(i=0; i<nobs; i++){
        for(p=0; p<P; p++){
            for(l=0; l<L; l++){
                Ml = BF.Ml[l];
                for(k=0; k<K; k++){
                    
                    tempsum = 0.0;
                    for(m=0; m<Ml; m++){
                        tempsum += PostSamp.xb[l][i*P*Ml + p*Ml +m] * PostSamp.alphastar[l][k*Ml+m];
                    }
                    PostSamp.xba[l][i*P*K+p*K+k] = tempsum;
                }
            }
        }
    }
    
    // Compute xbar
    for(i=0; i<nobs; i++){
        for(l=0; l<L; l++){
            Ml = BF.Ml[l];
            for(k=0; k<K; k++){
                tempsum = 0.0;
                for(p=0; p<P; p++){
                    tempsum += PostSamp.xba[l][i*P*K+p*K+k] * PostSamp.gamma[p*L+l];
                }
                PostSamp.xbar[l][i*K+k] = tempsum;
            }
        }
    }
    
    // compute zb2= (Z-Zinter)*basis
    int parcel_len;
    for(l=0; l<L; l++){
        Ml = BF.Ml[l];
        parcel_len = data.parcel_len[l];
        for(i=0; i<nobs; i++){
            for(j=0; j<Ml; j++){
                tempsum=0.0;
                for(k=0; k<parcel_len; k++){
                    tempsum += (double)PostSamp.Zinter[l][k] * (double)BF.basis[l][k*Ml+j];
                }
                PostSamp.zb2[l][i*Ml+j] = PostSamp.zb[l][i*Ml+j] - tempsum;
            }
        }
    }
    
    // Training
    for(i=0; i<nobs; i++){
        for(l=0; l<L; l++){
            Ml = BF.Ml[l];
            parcel_len = data.parcel_len[l];
            for(p=0; p<parcel_len; p++){
                PostSamp.outcome_train_mean[l][i*parcel_len + p] = 0.0f;
                PostSamp.outcome_train_mean2[l][i*parcel_len + p] = 0.0f;
            }
        }
    }
    for(l=0; l<L; l++){
        Ml = BF.Ml[l];
        for(i=0; i<nobs; i++){
            for(m=0; m<Ml; m++){
                PostSamp.theta_train_mean[l][i*Ml+m] = 0.0f;
            }
        }
    }
    for(i=0; i<nobs; i++){
        for(l=0; l<L; l++){
            for(k=0; k<K; k++){
                PostSamp.latent_train_mean[l][i*K+k] = 0.0f;
            }
        }
    }
    
    // Printout
    if(printInit == true){
        Rprintf("************ Initial Values of Model Parameters ************\n");
        print_t0(PostSamp, data, BF, P, K, L);
    }
    
    
    // Output
    out_t0(PostSamp, data, BF, L, nobs, nts, K, P, r, outpath);
}
    
    
    
void set_initial(const int L, const int nobs, const int nts, const int K, const int P, struct Sampling PostSamp, const struct Inputdata data, const struct BasisFunc BF, gsl_rng *r,  const char *outpath, const bool printInit){
    
    gsl_set_error_handler_off();
    // Errors
    PostSamp.err2inv_e[0] = gsl_ran_flat(r, 1.0, 10.0);
    PostSamp.err2inv_zeta[0] = gsl_ran_flat(r, 1.0, 10.0);
    PostSamp.err2inv_eps[0] =  gsl_ran_flat(r, 1.0, 10.0);
    PostSamp.err2inv_e_mean[0] = 0.0f;
    PostSamp.err2inv_zeta_mean[0] = 0.0f;
    PostSamp.err2inv_eps_mean[0] = 0.0f;
    int l;
    for(l=0; l<L; l++){
        PostSamp.err2inv_u[l] = gsl_ran_flat(r, 1.0, 10.0);
        PostSamp.err2inv_u_mean[l] = 0.0f;
    }
    
    // Theta
    //init_theta_rand(nobs, L, PostSamp, data, BF);
    
    // Phi2Inv
    init_phi_rand(PostSamp, K, L, r);

    // Alpha
    init_alpha_rand(PostSamp, data, BF, K, L, r);
    
    // Beta
    init_beta_rand(PostSamp, BF, data, K, L);
    
    // Weight
    init_w_rand(PostSamp, L);
    
    // Gamma
    init_gamma_rand(PostSamp, P, L, r);
    
    // Summarized X
    init_sumx(PostSamp, data, BF, nobs, L, P);
    
    // Mu
    init_mu_rand(PostSamp, L, nobs, K, r);
    
    // Latent factors
    init_latent_rand(PostSamp, BF, data, K, L, nobs, r);
    
    // Loading matrix
    init_load_rand(PostSamp, BF, L, K, nobs);
    
    // Compute \sum_p gamma*basis*X
    int i, j, m, p, k, Ml;
    double tempsum;
    for(i=0; i<nobs; i++){
        for(l=0; l<L; l++){
            Ml = BF.Ml[l];
            for(m=0; m<Ml; m++){
                tempsum = 0.0;
                for(p=0; p<P; p++){
                    tempsum += PostSamp.gamma[p*L+l] * PostSamp.xb[l][i*P*Ml + p*Ml+ m];
                }
                PostSamp.rxb[l][i*Ml+m] = tempsum;
            }
        }
    }

    
    // Compute xba
    for(i=0; i<nobs; i++){
        for(p=0; p<P; p++){
            for(l=0; l<L; l++){
                Ml = BF.Ml[l];
                for(k=0; k<K; k++){
                    
                    tempsum = 0.0;
                    for(m=0; m<Ml; m++){
                        tempsum += PostSamp.xb[l][i*P*Ml + p*Ml +m] * PostSamp.alphastar[l][k*Ml+m];
                    }
                    PostSamp.xba[l][i*P*K+p*K+k] = tempsum;
                }
            }
        }
    }
    
    // Compute xbar
    for(i=0; i<nobs; i++){
        for(l=0; l<L; l++){
            Ml = BF.Ml[l];
            for(k=0; k<K; k++){
                tempsum = 0.0;
                for(p=0; p<P; p++){
                    tempsum += PostSamp.xba[l][i*P*K+p*K+k] * PostSamp.gamma[p*L+l];
                }
                PostSamp.xbar[l][i*K+k] = tempsum;
            }
        }
    }

    // compute zb2= (Z-Zinter)*basis
    int parcel_len;
    for(l=0; l<L; l++){
        Ml = BF.Ml[l];
        parcel_len = data.parcel_len[l];
        for(i=0; i<nobs; i++){
            for(j=0; j<Ml; j++){
                tempsum=0.0;
                for(k=0; k<parcel_len; k++){
                    tempsum += (double)PostSamp.Zinter[l][k] * (double)BF.basis[l][k*Ml+j];
                }
                PostSamp.zb2[l][i*Ml+j] = PostSamp.zb[l][i*Ml+j] - tempsum;
            }
        }
    }
    
    // Training
    for(i=0; i<nobs; i++){
        for(l=0; l<L; l++){
            Ml = BF.Ml[l];
            parcel_len = data.parcel_len[l];
            for(p=0; p<parcel_len; p++){
                PostSamp.outcome_train_mean[l][i*parcel_len + p] = 0.0f;
                PostSamp.outcome_train_mean2[l][i*parcel_len + p] = 0.0f;
            }
        }
    }
    for(l=0; l<L; l++){
        Ml = BF.Ml[l];
        for(i=0; i<nobs; i++){
            for(m=0; m<Ml; m++){
                PostSamp.theta_train_mean[l][i*Ml+m] = 0.0f;
            }
        }
    }
    for(i=0; i<nobs; i++){
        for(l=0; l<L; l++){
            for(k=0; k<K; k++){
                PostSamp.latent_train_mean[l][i*K+k] = 0.0f;
            }
        }
    }    
    
    // Printout
    if(printInit == true){
        Rprintf("************ Initial Values of Model Parameters ************\n");
        print_t0(PostSamp, data, BF, P, K, L);
    }
    
    
    // Output
    out_t0(PostSamp, data, BF, L, nobs, nts, K, P, r, outpath);
    
}

void init_zInter(const int nobs, const int L, struct Sampling PostSamp,
                 const struct Inputdata data, const struct BasisFunc BF){
    int i, j, l;
    gsl_set_error_handler_off();
    for (l=0; l<L; l++) {
        int parcel_len = data.parcel_len[l];
        for(j=0; j<parcel_len; j++) {
            double tempsum = 0.0;
            for(i=0; i<nobs; i++){
                tempsum += (double)data.Zl[l][i*parcel_len+j];
            }
            PostSamp.Zinter[l][j] = (float) (tempsum / (double)nobs);
        }

    }
}


void init_theta_rand(const int nobs, const int L, struct Sampling PostSamp, const struct Inputdata data, const struct BasisFunc BF){
    gsl_set_error_handler_off();
    int i, j, k, l, Ml;
    double tempsum;
    
    // estimates of theta_k^l (k=1,...,nobs; l=1,...,L)
    for(l=0; l<L; l++){
        Ml = BF.Ml[l];
        gsl_matrix * btb = gsl_matrix_alloc(Ml, Ml);
        gsl_matrix * btb_inv = gsl_matrix_alloc(Ml, Ml);
        
        for(i=0; i<Ml; i++){
            for(j=0; j<Ml; j++){
                tempsum = PostSamp.basis2[l][i*Ml+j];
                gsl_matrix_set(btb, i, j, tempsum);
                gsl_matrix_set(btb_inv, i, j, 0.0f);
            }
        }
        // inverse of btb
        gsl_matrix_inv(btb, btb_inv);
        
        for(k=0; k<nobs; k++){
            for(i=0; i<Ml; i++){
                tempsum = 0.0;
                for(j=0; j<Ml; j++){
                    tempsum += gsl_matrix_get(btb_inv, i, j) * (double)PostSamp.zb[l][k*Ml+j];
                }
                PostSamp.theta[l][k*Ml+i] = tempsum;
                PostSamp.mean_theta[l][k*Ml+i] = 0.0f;
            }
        }
        gsl_matrix_free(btb);
        gsl_matrix_free(btb_inv);
    }
}


void init_latent_rand(struct Sampling PostSamp, const struct BasisFunc BF, const struct Inputdata data,
                      const int K, const int L, const int nobs, gsl_rng*r){
    int l, i, k, v, Ml, parcel_len;
    gsl_set_error_handler_off();
    // Initial Latent values
    double err, tempsum;
    for(l=0; l<L; l++){
        Ml = BF.Ml[l];
        parcel_len = data.parcel_len[l];
        for(i=0; i<nobs; i++){
            for(k=0; k<K; k++){
                tempsum = 0.0;
                for(v=0; v<parcel_len; v++){
                    tempsum += (double)PostSamp.sumX[l] [i*parcel_len+v] * (double)PostSamp.beta[l][v*K+k];
                }
                //err = gsl_ran_gaussian(r, sqrt(1.0/PostSamp.err2inv_eps[0]));
                PostSamp.latent[l][i*K+k] = gsl_ran_gaussian(r, 2.0);//(float)tempsum + err;
                PostSamp.latentstar[l][i*K+k] = PostSamp.latent[l][i*K+k] * sqrtf(PostSamp.Phi2Inv[k*L+l]) + PostSamp.mustar[l][i*K+k];
                PostSamp.mean_latent[l][i*K+k] = 0.0f;
                PostSamp.mean_latentstar[l][i*K+k] = 0.0f;
            }
        }
        
        for(i=0; i<K; i++){
            for(k=0; k<K; k++){
                tempsum = 0.0;
                for(v=0; v<nobs; v++){
                    tempsum += PostSamp.latent[l][v*K+i] * PostSamp.latent[l][v*K+k];
                }
                PostSamp.latent2[l][i*K+k] = tempsum;
            }
        }
        
    }
}


void init_load_rand(struct Sampling PostSamp, struct BasisFunc BF, const int L, const int K, const int nobs){
    gsl_set_error_handler_off();
    gsl_matrix * ete = gsl_matrix_alloc(K, K);
    gsl_matrix * ete_inv = gsl_matrix_alloc(K, K);
    
    int l, i, j, k, m, Ml;
    double tempsum = 0.0;

    for(l=0; l<L; l++){
        Ml = BF.Ml[l];
        float * kmmat = (float *)calloc(K*Ml, sizeof(float));
        for(k=0; k<K; k++){
            for(j=0; j<K; j++){
                gsl_matrix_set(ete, k, j, PostSamp.latent2[l][k*K+j]);
            }
        }
        
        gsl_matrix_inv(ete, ete_inv);
        Ml = BF.Ml[l];
        for(k=0; k<K; k++){
            for(m=0; m<Ml; m++){
                tempsum = 0.0;
                for(i=0; i<nobs; i++){
                    tempsum += PostSamp.latent[l][i*K+k] * PostSamp.theta[l][i*Ml+m];
                }
                kmmat[k*Ml+m] = tempsum;
            }
        }
        
        for(k=0; k<K; k++){
            for(m=0; m<Ml; m++){
                tempsum = 0.0;
                for(i=0; i<K; i++){
                    tempsum += gsl_matrix_get(ete_inv, k, i) * kmmat[i*Ml+m];
                }
                PostSamp.load[l][m*K+k]=tempsum;
                PostSamp.loadstar[l][m*K+k]=tempsum;
                PostSamp.mean_loadstar[l][m*K+k] = 0.0f;
                PostSamp.mean_load[l][m*K+k] = 0.0f;
            }
        }
        free(kmmat);
    }
    gsl_matrix_free(ete);
    gsl_matrix_free(ete_inv);
}






void init_phi_rand(struct Sampling PostSamp, const int K,const int L, gsl_rng*r){
    int i;
    gsl_set_error_handler_off();
    for(i=0; i<K*L; i++){
        PostSamp.Phi2Inv[i] = gsl_ran_flat(r, 0.5, 2.0);
        PostSamp.mean_Phi2Inv[i] = 0.0f;
    }
}

void init_alpha_rand(struct Sampling PostSamp, const struct Inputdata data, const struct BasisFunc BF,
                     const int K, const int L, gsl_rng *r){
    gsl_set_error_handler_off();
    float s;
    // Initial alpha values
    
    int k, l, m, Ml;
    for(k=0; k<K; k++){
        for(l=0; l<L; l++){
            Ml = BF.Ml[l];
            // Sign of loading matrix
            s = PostSamp.loadstar[l][k*K+k];
            if(s >= 0){
                s= 1.0;
            }else{
                s = -1.0;
            }
            // alpha
            for(m=0; m<Ml; m++){
                PostSamp.alpha[l][k*Ml+m] = gsl_ran_gaussian(r, 1.0/sqrt(PostSamp.err2inv_alpha));
                PostSamp.alphastar[l][k*Ml+m] = PostSamp.alpha[l][k*Ml+m] * sqrt(PostSamp.Phi2Inv[k*L+l]) * s;
                PostSamp.mean_alpha[l][k*Ml+m] = 0.0f;
                PostSamp.mean_alphastar[l][k*Ml+m] = 0.0f;
            }
        }
    }
}



void init_beta_rand(struct Sampling PostSamp, const struct BasisFunc BF, const struct Inputdata data,
                    const int K, const int L){
    int k, l, i, v, Ml, parcel_len;
    float s;
    // Initial beta values
    double tempsum;
    for(k=0; k<K; k++){
        for(l=0; l<L; l++){
            Ml = BF.Ml[l];
            parcel_len = data.parcel_len[l];
            
            // Sign of loading matrix
            s = PostSamp.loadstar[l][k*K+k];
            if(s >= 0){
                s= 1.0;
            }else{
                s = -1.0;
            }
            
            for(v=0; v<parcel_len; v++){
                tempsum = 0.0;
                for(i=0; i<Ml; i++){
                    tempsum += PostSamp.alpha[l][k*Ml+i] * BF.basis[l][v*Ml+i];
                }
                PostSamp.beta[l][v*K+k] = tempsum;
                PostSamp.betastar[l][v*K+k] = PostSamp.beta[l][v*K+k] * PostSamp.Phi2Inv[k*L+l] * s;
                PostSamp.mean_beta[l][v*K+k] = 0.0f;
                PostSamp.mean_betastar[l][v*K+k] = 0.0f;
            }
        }
    }
}


void init_w_rand(struct Sampling PostSamp, const int L){
    int l;
    for(l=0; l<L; l++){
        PostSamp.w[l]=0.5;
        PostSamp.mean_w[l] = 0.0f;
    }
}


void init_gamma_rand(struct Sampling PostSamp, int P,  int L, gsl_rng *r){
    int p;
    gsl_set_error_handler_off();
    for(p=0; p<P*L; p++){
        PostSamp.gamma[p]= gsl_ran_bernoulli(r, .5);
        PostSamp.mean_gamma[p] = 0.0f;
        PostSamp.mean_gamma_update[p] = 0.0f;
    }
}



void init_sumx(struct Sampling PostSamp, const struct Inputdata data, const struct BasisFunc BF, const int nobs, const int L, const int P ){
    int i, l, j, p, Ml, parcel_len;
    double tempx;
    
    for(i=0; i<nobs; i++){
        for(l=0; l<L; l++){
            Ml = BF.Ml[l];
            parcel_len = data.parcel_len[l];
            for(j=0; j<parcel_len; j++){
                tempx = 0.0;
                for(p=0; p<P; p++){
                    tempx += PostSamp.gamma[p*L+l] * data.Xl[l][i*parcel_len*P + p*parcel_len + j];
                }
                PostSamp.sumX[l][i*parcel_len + j] = tempx;
            }
        }
    }
}


void init_mu_rand(struct Sampling PostSamp, const int L, const int nobs, const int K, gsl_rng*r){
    int i, k, l;
    gsl_set_error_handler_off();
    for(l=0; l<L; l++){
        for(i=0; i<nobs; i++){
            for(k=0; k<K; k++){
                PostSamp.mustar[l][i*K+k] = gsl_ran_gaussian(r, 1);
                PostSamp.mean_mu[l][i*K+k] = 0.0f;
            }
        }
    }
}
