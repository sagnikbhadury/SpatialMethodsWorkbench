void err2inv_e_samp(const int nobs, const int L, struct Sampling PostSamp, const struct Inputdata data, const struct BasisFunc BF, const gsl_rng *r){
    int i, j, l,m;
    gsl_set_error_handler_off();
    //Rprintf("\n########\nOld_Err2inv_e: ");
    //Rprintf("%.6f\n", PostSamp.err2inv_e[0]);
    
    double post_b = 0.0;
    double temp=0.0;
    double tempsum2=0.0;
    double tempsum=0.0;
    double tmpllh=0.0;
    double post_a = 0.0;
    
    int Ml, parcel_len, image_len;
    image_len = 0;
    
    for(l=0; l<L; l++){
        Ml = BF.Ml[l];
        parcel_len = data.parcel_len[l];
        image_len += parcel_len;
        
        // sampling
        post_a = (double)PostSamp.err2inv_e_a + (double) (parcel_len * nobs) * 0.5;
        post_b = PostSamp.err2inv_e_b;
        
        for(i=0; i<nobs; i++){
            //basis*theta
            tempsum2 = 0.0;
            for(j=0; j<parcel_len; j++){
                tempsum = 0.0;
                for(m=0; m<Ml; m++){
                    tempsum += BF.basis[l][j*Ml+m] * PostSamp.theta[l][i*Ml+m];
                }
                temp = data.Zl[l][i*parcel_len + j] - PostSamp.Zinter[l][j] - tempsum;
                tempsum2 += temp * temp;
            }
            post_b += 0.5 * tempsum2;
            tmpllh += 0.5 * tempsum2 * PostSamp.err2inv_e[0];
        }
        
        PostSamp.err2inv_e[0] =  gsl_ran_gamma(r, post_a, 1.0/post_b);
        //Rprintf("%.6e\n", post_a);
        //Rprintf("%.6e\n", post_b);
    }

    //Rprintf("New_Err2inv_e: ");
    //Rprintf("%.6f\n", PostSamp.err2inv_e[0]);
}



