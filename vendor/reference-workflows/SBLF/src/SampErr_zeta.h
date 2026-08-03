//
// Post sampling of the error term in submodel 2
//


void err2inv_zeta_samp(const int nobs, const int L, const int K,
                       struct Sampling PostSamp, const struct Inputdata data, const struct BasisFunc BF, const gsl_rng *r){
    gsl_set_error_handler_off();
    int i, j, l,m;
    
    //printf("\n########\nOld_Err2inv_zeta:");
    //printf("%.6f\n", PostSamp.err2inv_zeta[0]);
    
    double post_b = PostSamp.err2inv_zeta_b;
    double temp, tempsum2, tempsum, templlh;
    int Ml, M;
    M = 0;
    for(l=0; l<L; l++){
        M += BF.Ml[l];
    }
    
    templlh = 0.0;
    for(l=0; l<L; l++){
        Ml = BF.Ml[l];
        
        for(i=0; i<nobs; i++){
            //loadstar * latentstar
            tempsum2 = 0.0;
            for(j=0; j<Ml; j++){
                tempsum = 0.0;
                for(m=0; m<K; m++){
                    tempsum += PostSamp.loadstar[l][j*K+m] * PostSamp.latentstar[l][i*K+m];
                }
                temp = PostSamp.theta[l][i*Ml+j] - tempsum;
                tempsum2 += temp * temp;
            }
            post_b += 0.5 * tempsum2;
            templlh += 0.5 * tempsum2 * PostSamp.err2inv_zeta[0] ;
        }
    }
   
    // sampling err_zeta^(-2)
    double post_a = (double)PostSamp.err2inv_zeta_a + (double)(M * nobs) * 0.5;
    PostSamp.err2inv_zeta[0] =  gsl_ran_gamma(r, post_a, 1.0/post_b);
    
    //printf("New_Err2inv_zeta: ");
    //printf("%.6f\n", PostSamp.err2inv_zeta[0]);
}



