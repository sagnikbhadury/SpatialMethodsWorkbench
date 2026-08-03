void weight_samp(struct Sampling PostSamp, const int P, const int L, gsl_rng *r){
    gsl_set_error_handler_off();
    int p, l;
    double gm;
    double pp = (double)P;
    /*printf("old Weight\n");
    for(l=0; l<L; l++){
        printf("%.3f ", PostSamp.w[l]);
    }
    printf("\n");*/

    for(l=0; l<L; l++){
        gm = 0.0f;
        for(p=0; p<P; p++){
            gm += PostSamp.gamma[p*L+l];
        }  
        PostSamp.w[l] = gsl_ran_beta(r, gm+PostSamp.w_a, PostSamp.w_b + pp - gm);
       //PostSamp.w[p*L+l] = (gm+PostSamp.w_a) / (PostSamp.w_a+ PostSamp.w_b + pp);
    }
    
    /*printf("New Weight\n");
    for(l=0; l<L; l++){
        printf("%.3f ", PostSamp.w[l]);
    }
    printf("\n");*/
    
}


