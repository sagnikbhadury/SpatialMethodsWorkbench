

void phi2inv_samp(const int K, const int nobs, const int L,
                  struct Sampling PostSamp, const struct Inputdata data, const struct BasisFunc BF,
                  gsl_rng *r){
    gsl_set_error_handler_off();
    int k, i, j, l, Ml, parcel_len;
    double post_a, post_b, tempsum, tempxb, tempxb2;
    
    /*printf("\n#########\nOld_Phi\n");
    for(i=0; i<K; i++){
        for(j=0; j<L; j++){
            printf("%.3f ", PostSamp.Phi2Inv[i*L+j]);
        }
    }
    printf("\n");*/

    
    for(l=0; l<L; l++){
        Ml = BF.Ml[l];
        parcel_len = data.parcel_len[l];
        post_a = PostSamp.Phi2Inv_a + 0.5 * (float)(nobs);
        
        for(k=0; k<K; k++){
            post_b = PostSamp.Phi2Inv_b;
            tempsum = 0.0;
            for(i=0; i<nobs; i++){
                tempxb = 0.0;
                for(j=0; j<Ml; j++){
                    tempxb += PostSamp.alphastar[l][k*Ml+j] * PostSamp.rxb[l][i*Ml+j];
                }
                tempxb2 = tempxb + PostSamp.mustar[l][i*K+k] - PostSamp.latentstar[l][i*K+k];
                tempsum += tempxb2 * tempxb2 * PostSamp.err2inv_eps[0] * 0.5;
            }
            post_b = post_b+tempsum;

            // Sampling
            PostSamp.Phi2Inv[k*L+l] = 1.0 / gsl_ran_gamma(r, post_a, 1.0/post_b);
        }
    }
    
    /*printf("\nNew_Phi\n");
    for(i=0; i<K; i++){
        for(j=0; j<L; j++){
            printf("%.3f ", PostSamp.Phi2Inv[i*L+j]);
        }
    }
    printf("\n");*/

}


