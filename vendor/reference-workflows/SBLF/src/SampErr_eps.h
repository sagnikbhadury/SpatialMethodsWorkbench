//
// Post sampling of the error term in submodel 3
//



void err2inv_eps_samp(const int nobs, const int L, const int K,
                       struct Sampling PostSamp, const struct Inputdata data, const struct BasisFunc BF, const gsl_rng *r){
    gsl_set_error_handler_off();
    int i, j, l, m;
    
    //printf("\n########\nOld_Err2inv_eps: ");
    //printf("%.6f\n", PostSamp.err2inv_eps[0]);
    
    double post_a = (double)PostSamp.err2inv_eps_a + (double)(K * nobs*L) * 0.5;
    double post_b = (double)PostSamp.err2inv_eps_b;
    double temp, tempsum2, tempsum, templlh;

    int Ml;
    templlh = 0.0;
    for(l=0; l<L; l++){
        Ml = BF.Ml[l];
        
        for(i=0; i<nobs; i++){
            //loadstar * latentstar
            tempsum2 = 0.0;
            for(j=0; j<K; j++){
                tempsum = 0.0;
                for(m=0; m<Ml; m++){
                    tempsum += PostSamp.alphastar[l][j*Ml+m] * PostSamp.rxb[l][i*Ml+m];
                }
                temp = PostSamp.latentstar[l][i*K+j] - PostSamp.mustar[l][i*K+j]- tempsum;
                tempsum2 += temp * temp / PostSamp.Phi2Inv[j*L+l];
            }
            post_b += 0.5 * tempsum2;
            templlh += 0.5 * tempsum2 *PostSamp.err2inv_eps[0];
        }
    }
    PostSamp.err2inv_eps[0] =  gsl_ran_gamma(r, post_a, 1.0/post_b);
    
    //printf("New_Err2inv_eps: ");
    //printf("%.6f\n", PostSamp.err2inv_eps[0]);
    
}





