struct BasisFunc{
    float ** basis;
    int * Ml; // number of basis per parcel
    float ** kernel_loc;
    int M; // total number of basis for all parcels
};

struct BasisFunc genBasis(int L, const char * basispath, struct Inputdata data, const double bandwidth, const int dd){
    
    /*
     generate basis kernels
     */
    
    struct BasisFunc BF;
    
    BF.Ml = (int*)calloc(L, sizeof(int));
    BF.basis = (float **)malloc(L * sizeof(float *));
    BF.kernel_loc = (float **)malloc(L * sizeof(float *));
    int dim_image = data.dim_image;
    int l, i, j, k, m, h;
    int loc, minx, maxx, miny, maxy, minz, maxz, nx, ny, nz;
    int parcel_len_l;
    double dx, dy, dz, dist, tempdist;
    double cut = sqrt((double)(dd*dd)*2.0);
    
    
    // Generate basis kernel locations and values parcel by parcel
    for(l=0; l<L; l++){
        
        parcel_len_l = data.parcel_len[l];
        
        // axis range of each parcel
        int * parcel_loc_x = (int *)calloc(parcel_len_l, sizeof(int));
        int * parcel_loc_y = (int *)calloc(parcel_len_l, sizeof(int));
        int * parcel_loc_z = (int *)calloc(parcel_len_l, sizeof(int));
        loc = data.parcel_len_sum[l];
        minx = 1000; maxx = -1000;
        miny = 1000; maxy = -1000;
        if(dim_image == 3){
            minz = 1000; maxz = -1000;
        }else{
            minz = 0; maxz = 0;
        }
        
        for(i=0; i<parcel_len_l; i++){
            parcel_loc_x[i] = data.axes[(loc+i)*dim_image+0];
            parcel_loc_y[i] = data.axes[(loc+i)*dim_image+1];
            minx = minx < parcel_loc_x[i] ? minx : parcel_loc_x[i];
            miny = miny < parcel_loc_y[i] ? miny : parcel_loc_y[i];
            maxx = maxx > parcel_loc_x[i] ? maxx : parcel_loc_x[i];
            maxy = maxy > parcel_loc_y[i] ? maxy : parcel_loc_y[i];
            
            if(dim_image == 3){
                parcel_loc_z[i] = data.axes[(loc+i)*dim_image+2];
                minz = minz < parcel_loc_z[i] ? minz : parcel_loc_z[i];
                maxz = maxz > parcel_loc_z[i] ? maxz : parcel_loc_z[i];
            }
        }

        // number of basis kernel locations in each parcel
        nx = (int)ceil(((double)(maxx-minx)/(double)dd)) + 1;
        ny = (int)ceil(((double)(maxy-miny)/(double)dd)) + 1;
        if(dim_image == 3){
            nz = (int)ceil(((double)(maxz-minz)/(double)dd)) + 1;
        }else{
            nz = 1;
        }

        int xs[nx];
        int ys[ny];
        int zs[nz];
        for(i=0; i<nx; i++){
            xs[i] = minx + i*dd;
        }
        for(i=0; i<ny; i++){
            ys[i] = miny + i*dd;
        }
        if(dim_image == 3){
            for(i=0; i<nz; i++){
                zs[i] = minz + i*dd;
            }
        }
        
        
        // exclude kernel locations outside the parcel image range
        int * tmp_kernel_locx = (int *)calloc(nx*ny*nz, sizeof(int));
        int * tmp_kernel_locy = (int *)calloc(nx*ny*nz, sizeof(int));
        int * tmp_kernel_locz = (int *)calloc(nx*ny*nz, sizeof(int));
        
        h = 0; // count for the basis
        switch (dim_image){
            case 2:
                h = 0;
                for(i=0; i<nx; i++){
                    for(j=0; j<ny; j++){
                        
                        dist = 1000.0;
                        m = 0;
                        while(dist >= cut && m < parcel_len_l){
                            dx = parcel_loc_x[m]-xs[i];
                            dy = parcel_loc_y[m]-ys[j];
                            tempdist = sqrt(dx*dx + dy*dy);
                            m++;
                            dist = tempdist<dist ? tempdist:dist;
                        }
                        
                        if(dist <= cut){
                            tmp_kernel_locx[h] = xs[i];
                            tmp_kernel_locy[h] = ys[j];
                            h++;
                        }
                    }
                }
                break;
            case 3:
                h = 0;
                for(i=0; i<nx; i++){
                    for(j=0; j<ny; j++){
                        for(k=0; k<nz; k++){
                            
                            dist = 1000.0;
                            m = 0;
                            
                            while(dist >= cut && m < parcel_len_l){
                                dx = parcel_loc_x[m]-xs[i];
                                dy = parcel_loc_y[m]-ys[j];
                                dz = parcel_loc_z[m]-zs[k];
                                tempdist = sqrt(dx*dx + dy*dy + dz*dz);
                                m++;
                                dist = tempdist<dist ? tempdist:dist;
                            }
                            
                            if(dist <= cut){
                                tmp_kernel_locx[h] = xs[i];
                                tmp_kernel_locy[h] = ys[j];
                                tmp_kernel_locz[h] = zs[k];
                                h++;
                            }
                            
                        }
                    }
                }
                break;
        }
        BF.Ml[l] = h;

        
        // calculate basis values
        BF.kernel_loc[l] = (float *)malloc(h * dim_image * sizeof(float));
        BF.basis[l] = (float *)malloc(h*parcel_len_l * sizeof(float));
        if(dim_image == 3){
            for(i=0; i<h; i++){
                BF.kernel_loc[l][i*3+0] = (float)tmp_kernel_locx[i];
                BF.kernel_loc[l][i*3+1] = (float)tmp_kernel_locy[i];
                BF.kernel_loc[l][i*3+2] = (float)tmp_kernel_locz[i];
                for(m=0; m<parcel_len_l; m++){
                    dx = parcel_loc_x[m]-tmp_kernel_locx[i];
                    dy = parcel_loc_y[m]-tmp_kernel_locy[i];
                    dz = parcel_loc_z[m]-tmp_kernel_locz[i];
                    dist = dx*dx+dy*dy+dz*dz;
                    BF.basis[l][m*h+i] = exp(-bandwidth*dist);
                }
            }
        }
        
        if(dim_image == 2){
            for(i=0; i<h; i++){
                BF.kernel_loc[l][i*2+0] = (float)tmp_kernel_locx[i];
                BF.kernel_loc[l][i*2+1] = (float)tmp_kernel_locy[i];
                for(m=0; m<parcel_len_l; m++){
                    dx = parcel_loc_x[m]-tmp_kernel_locx[i];
                    dy = parcel_loc_y[m]-tmp_kernel_locy[i];
                    dist = dx*dx+dy*dy;
                    BF.basis[l][m*h+i] = exp(-bandwidth*dist);
                }
            }
        }
        
        
        free(parcel_loc_x);
        free(parcel_loc_y);
        free(parcel_loc_z);
        free(tmp_kernel_locx);
        free(tmp_kernel_locy);
        free(tmp_kernel_locz);
    }
    
    // total basis kernels of all parcels
    BF.M = 0;
    for(l=0; l<L; l++){
        BF.M += BF.Ml[l];
    }

    
    ///// Output Basis kernel locations and vasis values
    char * outb = (char*)calloc(500, sizeof(char));
    char sl[20];
    strcat(outb, basispath);
    strcat(outb, "Basis_");
    
    for(l=0; l<L; l++){
        char * pb = (char*)calloc(500, sizeof(char));
        strcpy(pb, outb);
        sprintf(sl, "%d", l+1);
        strcat(pb, sl);
        strcat(pb, ".txt");
        FILE * fb = fopen(pb, "wb");
        if(fb == NULL){
            Rprintf("Cannot open file for basis.\n");
        }
        for(i=0; i<data.parcel_len[l]; i++){
            for(j=0; j<BF.Ml[l]; j++){
                fprintf(fb, "%f ", BF.basis[l][i*BF.Ml[l]+j]);
            }
            fprintf(fb, "\n");
        }
        fclose(fb);
        free(pb);
        
        char * pbk = (char*)calloc(500, sizeof(char));
        strcpy(pbk, outb);
        strcat(pbk, "Kernel_");
        strcat(pbk, sl);
        strcat(pbk, ".txt");
        FILE * fbk = fopen(pbk, "wb");
        if(fbk == NULL){
            Rprintf("Cannot open file for basis_kernel.\n");
        }
        for(i=0; i<BF.Ml[l]; i++){
            for(j=0; j<dim_image; j++){
                fprintf(fbk, "%f ", BF.kernel_loc[l][i*dim_image+j]);
            }
            fprintf(fbk, "\n");
        }
        fclose(fbk);
        free(pbk);
    }
    
    return BF;
}
