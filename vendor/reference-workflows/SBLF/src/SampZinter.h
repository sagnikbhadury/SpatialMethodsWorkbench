//
// Sampling the intercept in submodel 1
//

void err2inv_u_samp(const int L, struct Sampling PostSamp, const struct Inputdata data, const gsl_rng *r);
void Zinter_samp(int nobs, int L, int K, struct Sampling PostSamp, struct Inputdata data, const struct BasisFunc BF, gsl_rng *r, int * singular);


//////////////////////////////////////////////////////////////////////////////////////////
void Zinter_samp(int nobs, int L, int K, struct Sampling PostSamp, struct Inputdata data, const struct BasisFunc BF, gsl_rng *r, int * singular){
    gsl_set_error_handler_off();
    int i, j, l, k, rr, Ml, parcel_len;
    double tmp, tempsum;
    
    /*printf("########\nOld_Intercept:\n");
    for(i=0; i<3; i++){
        for(j=0; j<2; j++){
            printf("%.3f ", PostSamp.Zinter[0][i* data.parcel_len[0]+j]);
        }
        printf("\n");
    }
    printf("\n");*/
    
    
    // Gibbs Sampling
    for(l=0; l<L; l++){
        parcel_len = data.parcel_len[l];
        Ml = BF.Ml[l];
        double var = 1.0 / ((double)nobs * PostSamp.err2inv_e[0] + PostSamp.err2inv_u[l]);
        double var2 = var * (double)PostSamp.err2inv_e[0];
        // compute precision matrix
        for(i=0; i<parcel_len; i++){
            // mean
            tempsum = 0.0;
            for (k=0; k<nobs; k++) {
                //theta * basis
                tmp = 0.0;
                for (rr=0; rr<Ml; rr++) {
                    tmp += (double)PostSamp.theta[l][k*Ml+rr] * (double)BF.basis[l][i*Ml+rr];
                }
                tempsum += (double)data.Zl[l][k*parcel_len+i] - tmp;
            }
            tempsum *= var2;
            PostSamp.Zinter[l][i] = gsl_ran_gaussian(r, sqrtf((float)var)) + tempsum ;
        }
    }

    // Update zb2= (Z-Zinter)*basis
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
    
    /*printf("New Intercepts:\n");
    for(i=0; i<3; i++){
        for(j=0; j<2; j++){
            printf("%.3f ", PostSamp.Zinter[0][i* data.parcel_len[0]+j]);
        }
        printf("\n");
    }*/
    
}


void err2inv_u_samp(const int L, struct Sampling PostSamp, const struct Inputdata data, const gsl_rng *r){
    gsl_set_error_handler_off();
    int l, v, parcel_len;
    double post_a, post_b;
    post_a = (double)PostSamp.err2inv_u_a + 0.5;
    
    for(l = 0; l<L; l++){
        parcel_len = data.parcel_len[l];
        post_b = (double)PostSamp.err2inv_u_b;
        for(v=0; v<parcel_len; v++){
            post_b += PostSamp.Zinter[l][v] * PostSamp.Zinter[l][v] * 0.5;
        }
        PostSamp.err2inv_u[l] =  gsl_ran_gamma(r, post_a, 1.0/post_b);
    }
}
