//
// Post sampling of beta
//

void beta_samp_approx(const int nobs, const int L, const int K,
               struct Sampling PostSamp, const struct Inputdata data, const struct BasisFunc BF, gsl_rng *r){
    gsl_set_error_handler_off();
    int i, j, l, k, Ml, parcel_len;

    /*printf("\n#########\nOld_BetaStar\n");
    for(i=0; i<3; i++){
        for(j=0; j<K; j++){
            printf("%.3f ", PostSamp.betastar[0][i*K+j]);
        }
        printf("\n");
    }*/
    
    // Beta = basis * alpha
    double tempsum;
    for(k=0; k<K; k++){
        for(l=0; l<L; l++){
            Ml = BF.Ml[l];
            parcel_len = data.parcel_len[l];
            for(j=0; j<parcel_len; j++){
                tempsum = 0.0;
                for(i=0; i<Ml; i++){
                    tempsum += BF.basis[l][j*Ml+i] * PostSamp.alphastar[l][k*Ml+i];
                }
                PostSamp.betastar[l][j*K+k]=tempsum;
            }
        }
    }
    
    // Transfer back to beta
    for(k=0; k<K; k++){
        for(l=0; l<L; l++){
            Ml = BF.Ml[l];
            parcel_len = data.parcel_len[l];
            double diag = PostSamp.loadstar[l][k*K + k];
            double sign;
            if(diag >=0){
                sign = 1.0;
            }else{
                sign = -1.0;
            }
            double phi = 1.0/sqrt((double)PostSamp.Phi2Inv[k*L+l]);
            
            for(j=0; j<parcel_len; j++){
                PostSamp.beta[l][j*K + k] = sign * phi * PostSamp.betastar[l][j*K + k];
            }
        }
    }
    
    /*printf("\nNew_BetaStar\n");
    for(i=0; i<3; i++){
        for(j=0; j<K; j++){
            printf("%.3f ", PostSamp.betastar[0][i*K+j]);
        }
        printf("\n");
    }*/

}


