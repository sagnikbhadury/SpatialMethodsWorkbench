struct Inputdata{
    int * sizes;
    float * Z;
    float ** Zl;
    float * X;
    float ** Xl;
    int * parcel_len;
    int * parcel_len_sum;
    int * axes;
    float * Z_test;
    float ** Zl_test;
    float * X_test;
    float ** Xl_test;
    float * wt_test_avg;
    float * wt_test;
    float ** wt_voxel;
    int dim_image;
};


struct Inputdata input(const char * inpathx, const char * inpath, const int sim){
    struct Inputdata data;
    int i, j, l, parcel_len;
    
    ///// Input sizes of images and datasets
    data.sizes = (int *)calloc(15, sizeof(int));
    char *sizes = (char *)calloc(500,sizeof(char));
    strcat(sizes, inpathx);
    strcat(sizes, "ROI_sizes.txt");
    FILE *rsizes = fopen(sizes, "r");
    
    if(rsizes == NULL){
        Rprintf("Error");
        Rprintf("Cannot open file for ROI_sizes.txt\n");
        exit(1);
    }
    for(i=0; i<5; i++) {
        if(fscanf(rsizes, "%d", &data.sizes[i]) !=1 ){
            break;
        }
        if(fscanf(rsizes, "\n") != 0){
            break;
        }
    }
    fclose(rsizes);
    free(sizes);

    int L = data.sizes[0]; // number of parcels
    int nobs = data.sizes[1]; // number of subjects
    int image_len = data.sizes[2]; // total number of voxels for the whole parcels
    int P = data.sizes[3]; // number of predictors
    int nts = data.sizes[4];  // number of testing subjects
    
    data.parcel_len = (int *)calloc(L, sizeof(int));
    data.parcel_len_sum = (int *)calloc(L, sizeof(int));
    if (L == 1){
        data.parcel_len[0] = image_len;
        data.parcel_len_sum[0] = 0;
    }else{
        Rprintf("Multiple parcels are not workable. Currently the code works for a single parcel only!\n");
        exit(0);
    }
    
    Rprintf("Data Info: \n");
    Rprintf("    sample size: %d (training), %d (test)\n", nobs, nts);
    Rprintf("    voxels: %d\n", image_len);
    Rprintf("    parcels: %d\n", L);
    Rprintf("    predictors: %d\n", P);
    
    ///// Input axes of images
    int dim_image = 3; // images from real data are 3-Dimensional
    dim_image -= sim; // simulated images are 2-dimensional
    data.dim_image = dim_image;
    
    data.axes = (int *)calloc(dim_image*image_len, sizeof(int));
    char *sa = (char *)calloc(500,sizeof(char));
    strcat(sa, inpathx);
    strcat(sa, "ROI_axes.txt");
    
    FILE *ra = fopen(sa, "r");
    if(ra == NULL){
        Rprintf("Cannot open file for ROI_axes.txt\n");
    }
    for(i=0; i<image_len; i++){
        for(j=0; j<dim_image; j++){
            if(fscanf(ra, "%d", &data.axes[i*data.dim_image+j]) !=1 ){
                break;
            }
        }
        if (fscanf(ra, "\n") != 0){
            break;
        }
    }
    fclose(ra);
    free(sa);

    
    ///// Input task images - outcomes
    Rprintf("Reading imaging outcomes\n");
    
    // training
    data.Z = (float *)calloc(nobs*image_len, sizeof(float));
    char *sz = (char *)calloc(500,sizeof(char));
    strcat(sz, inpath);
    strcat(sz, "ROI_task.txt");
    FILE *rz = fopen(sz, "r");
    if(rz == NULL){
        Rprintf("Cannot open file for ROI_task.txt\n");
    }
    for(i=0; i<nobs; i++){
        for(j=0; j<image_len; j++){
            if(fscanf(rz, "%f", &data.Z[i*image_len+j]) !=1 ){
                break;
            }
        }
        if (fscanf(rz, "\n") != 0){
            break;
        }
    }
    fclose(rz);
    free(sz);

    // sorted by parcels
    data.Zl = (float **)malloc(L * sizeof(float *));
    int parcel_len_sum_l;
    
    for(l=0; l<L; l++){
        parcel_len = data.parcel_len[l];
        parcel_len_sum_l = data.parcel_len_sum[l];
        data.Zl[l] = (float *)malloc(nobs*parcel_len * sizeof(float));
        for(i=0; i<nobs; i++){
            for(j=0; j<parcel_len; j++){
                data.Zl[l][i*parcel_len+j] = data.Z[i*image_len+parcel_len_sum_l+j];
            }
        }
    }
    
    // testing
    if(nts > 0){
        data.Z_test = (float *)calloc(nts*image_len, sizeof(float));
        char *szts = (char *)calloc(500,sizeof(char));
        strcat(szts, inpath);
        strcat(szts, "ROI_task_test.txt");
        FILE *rzts = fopen(szts, "r");
        if(rzts == NULL){
            Rprintf("Cannot open file for ROI_task.txt for testing set\n");
        }
        for(i=0; i<nts; i++){
            for(j=0; j<image_len; j++){
                if(fscanf(rzts, "%f", &data.Z_test[i*image_len+j]) !=1 ){
                    break;
                }
            }
            if (fscanf(rzts, "\n") != 0){
                break;
            }
        }
        fclose(rzts);
        free(szts);
        
        // sorted by parcels
        data.Zl_test = (float **)malloc(L * sizeof(float *));
        for(l=0; l<L; l++){
            parcel_len = data.parcel_len[l];
            parcel_len_sum_l = data.parcel_len_sum[l];
            data.Zl_test[l] = (float *)malloc(nts*parcel_len * sizeof(float));
            for(i=0; i<nts; i++){
                for(j=0; j<parcel_len; j++){
                    data.Zl_test[l][i*parcel_len+j] = data.Z_test[i*image_len+parcel_len_sum_l+j];
                }
            }
        }
    }
    
    
    
    ///// Input Feature images - predictors
    Rprintf("Reading imaging predictors\n");
    
    // training
    data.X = (float *)calloc(nobs*image_len*P, sizeof(float));
    char *sx = (char *)calloc(500,sizeof(char));
    strcat(sx, inpathx);
    strcat(sx, "ROI_dat_mat.txt");
    FILE *rx = fopen(sx, "r");
    if(rx == NULL){
        Rprintf("Cannot open file for ROI_dat_mat.txt\n");
    }
    for(i=0; i<nobs; i++){
        for(j=0; j<image_len*P; j++){
            if(fscanf(rx, "%f", &data.X[i*image_len*P+j]) !=1 ){
                break;
            }
        }
        if (fscanf(rx, "\n") != 0){
            break;
        }
    }
    fclose(rx);
    free(sx);
    
    // sorted by parcels
    data.Xl = (float **)malloc(L * sizeof(float *));
    int p;
    for(l=0; l<L; l++){
        parcel_len = data.parcel_len[l];
        parcel_len_sum_l = data.parcel_len_sum[l];
        data.Xl[l] = (float *)malloc(nobs * parcel_len * P * sizeof(float));
        for(i=0; i<nobs; i++){
            for(p=0; p<P; p++){
                for(j=0; j<parcel_len; j++){
                    data.Xl[l][i*parcel_len*P+p*parcel_len+j] = data.X[i*image_len*P + p*image_len + parcel_len_sum_l + j];
                }
            }
        }
    }

    // testing
    if(nts > 0){
        data.X_test = (float *)calloc(nts*image_len*P, sizeof(float));
        char *sxts = (char *)calloc(500,sizeof(char));
        strcat(sxts, inpathx);
        strcat(sxts, "ROI_dat_mat_test.txt");
        FILE *rxts = fopen(sxts, "r");
        if(rxts == NULL){
            Rprintf("Cannot open file for ROI_dat_mat.txt for testing set\n");
        }
        for(i=0; i<nts; i++){
            for(j=0; j<image_len*P; j++){
                if(fscanf(rxts, "%f", &data.X_test[i*image_len*P+j]) !=1 ){
                    break;
                }
            }
            if (fscanf(rxts, "\n") != 0){
                break;
            }
        }
        fclose(rxts);
        free(sxts);
        
        // sorted by parcels
        data.Xl_test = (float **)malloc(L * sizeof(float *));
        for(l=0; l<L; l++){
            parcel_len = data.parcel_len[l];
            parcel_len_sum_l = data.parcel_len_sum[l];
            data.Xl_test[l] = (float *)malloc(nts * parcel_len * P * sizeof(float));
            for(i=0; i<nts; i++){
                for(p=0; p<P; p++){
                    for(j=0; j<parcel_len; j++){
                        data.Xl_test[l][i*parcel_len*P+p*parcel_len+j] = data.X_test[i*image_len*P + p*image_len + parcel_len_sum_l + j];
                    }
                }
            }
        }
    }
    
    return data;
}
