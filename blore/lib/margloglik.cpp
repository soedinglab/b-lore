#define DLLEXPORT extern "C"
#include <stdio.h>	/* printf */
#include <math.h>	/* log, exp */
#include <stdlib.h>	/* abs */

void compute_grad(int nvar,
                 int nhyp,
                 int nfeat,
                 int zlen,
                 double *pi,
                 double *mu,
                 double *sig2,
                 double *feat_flat,
                 int *zarr,
                 int *znorm,
                 double *zcomps,
                 double *intmu,
                 double *intlam_inv_flat,
                 bool is_covariate,
                 double *grad
                );

void gauss_jordan_inversion(int n, double **matrix, double &logdet, double &signdet) {

    int i, j, k;
    double ratio;
    double a;

    /* hardcoded for fast calculation of single elements */
    if (n == 1) {
        a = matrix[0][0];
        matrix[0][1] = 1 / a;
        logdet = log(abs(a));
        signdet = 1;
        if (a < 0) {
            signdet *= -1;
        }

    } else {
    
        /* creating the identity matrix beside the original one */
        for(i = 0; i < n; ++i){
            for(j = n; j < 2*n; ++j){
                if(i==(j-n))
                    matrix[i][j] = 1.0;
                else
                    matrix[i][j] = 0.0;
            }
        }
    
        /* row operations to reduce original matrix to a diagonal matrix */
        for(i = 0; i < n; ++i){
            for(j = 0; j < n; ++j){
                if(i != j){
                    ratio = matrix[j][i]/matrix[i][i];
                    for(k = 0; k < 2*n; ++k){
                        matrix[j][k] -= ratio * matrix[i][k];
                    }
                }
            }
        }
    
        /* reducing to unit matrix, determinant is calculated from the diagonal matrix before the reduction */
        logdet = 0.0;
        signdet = 1;
        for(i = 0; i < n; ++i){
            a = matrix[i][i];
            if (a < 0) {
                signdet *= -1;
                logdet += log(-a);
            } else {
                logdet += log(a);
            }
            for(j = 0; j < 2*n; ++j){
                matrix[i][j] /= a;
            }
        }
    }
}


DLLEXPORT void margloglik_zcomps(int nvar,			// number of SNPs
                                int nhyp,			// number of hyperparameters (required only for gradient)
                                int nfeat,                      // number of features (required only for gradient)
                                int zlen,			// number of zstates
                                int zmax,			// max norm of the zstates
                                double *pi,
                                double *mu,
                                double *sig2,
                                double *v,
                                double *precflat,		// flattened array of precision matrix
                                double *feat_flat,		// flattened array of features (required only for gradient)
                                int *zarr,			// array containing the causal picks of all zstates
                                int *znorm,			// array with the number of causal picks of all zstates
                                double reg2,
                                double mureg,
                                double *zcomps,			// output
                                double *grad,                   // output
                                bool get_gradient,
                                bool is_covariate
                               ) {

    int i, j;
    int ki, kj;
    int z;
    int dim; // an integer to hold znorm[i]
    int zindx; // count of number of states already considered
    int sqindx;
    double ai, aij;
    double sumlogpi;
    double hscale; // a large number of the order of sig2 for scaling Lamda_z.
    double precdotv[nvar];
    double intlam[nvar][nvar];
    double *intmu;
    double *intlam_inv_flat;
    double logpz, logfz;
    double *amat[zmax];
    double logdet, signdet;
    double vzinner[zmax];
    double fzterm1, fzterm2, fzterm3, fzterm4;
    double fztmp;
    double tmp;

    //printf ("%s %i zstates\n", "In C++", zlen);
    //printf ("In C++, sig2 = %g\n", sig2[0]);

    // How many elements will there be in intmu and intlam_inv?
    ki = 0;
    kj = 0;
    for (i = 0; i < zlen; ++i) {
        ai = znorm[i];
        ki += ai;
        kj += ai * ai;
    }
    intmu = new double[ki];
    intlam_inv_flat = new double[kj];
    ki = 0;
    kj = 0;

    /* Calculate the sum of log (1 - pi) only once
       Calculate Lambda_z only once, but note: this will give speed benefit only when zlen > n^2, which is possible if znorm > 2
       We also precalculate prec.dot.v 
    */
    sumlogpi = 0.0;
    if (!is_covariate) {
        for (i = 0; i < nvar; ++i) {
            tmp = log(1 - pi[i]);
            if (-tmp == INFINITY) {
                printf ("SNP %i: %g %g\n", i + 1, tmp, pi[i]);
            }
            sumlogpi += tmp;
        }
    }
    if ((sumlogpi == INFINITY) or (-sumlogpi == INFINITY)) {
        printf ("Sum of logpi: %g\n", sumlogpi);
        exit(0);
    }

    for (i = 0; i < nvar; ++i) {
        precdotv[i] = 0.0;
        for (j = 0; j < nvar; ++j) {
            aij = precflat[i * nvar + j];
            intlam[i][j] = aij;
            precdotv[i] += aij * v[j];
        }
        intlam[i][i] += (1 / sig2[i]) - (1 / reg2);
    }

    //printf ("Lambdaz[0][0]: %g\n", intlam[0][0]);
    // Start the loop over all zstates for the calculation of marginal likelihood
    zindx = 0;
    sqindx = 0;
    for (z = 0; z < zlen; ++z) {


        dim = znorm[z];
        //printf ("%i %i\n", z, dim);
        /* if znorm is 0, then we need to calculate only zcomps. 
           intmu and intlam are empty, hence not updated.
           Note that zindx and sqindx also does not change.
        */
        if (dim == 0) {
            zcomps[z] = exp(sumlogpi + 1.0);

        } else {

            // Log of probability of z given Xi and Theta => ln p(z | Xi, Theta)
            logpz = sumlogpi;
            if (!is_covariate) {
                for (i = 0; i < dim; ++i) {
                    ki = zarr[zindx + i];
                    //logpz += log (pi[ki] / (1 - pi[ki]));
                    tmp = log (pi[ki] / (1 - pi[ki]));
                    if ((tmp != tmp) or (tmp == INFINITY) or (-tmp == INFINITY)) {
                        printf ("Error in pi value for SNP %i: %g\n", i, pi[ki]);
                        exit(0);
                    }
                    logpz += tmp;
                }
            }
            //printf ("Debug logpz: %g %g\n", logpz, sumlogpi);

            /* Select Lambda_z and invert
               Note that amat should hold the original matrix A,
               plus enough space for an identity matrix I of same dimension.
               A will be irreversibly changed to I.
               We also get the sign and log of the determinant of A.
               and I will change to  inv(A)
               After inversion, we retrieve intlam_inv from the original position of I
            */
            hscale = 1 / sig2[zarr[zindx]];
            for (i = 0; i < dim; ++i) {
                ki = zarr[zindx + i];
                amat[i] = new double[2*dim];
                for (j = 0; j < dim; ++j) {
                    kj = zarr[zindx + j];
                    amat[i][j] = intlam[ki][kj] / hscale;
                }
            }
            //gauss_jordan_inversion(dim, amat, logdet, signdet, is_covariate, sig2);
            gauss_jordan_inversion(dim, amat, logdet, signdet);
            if (signdet < 0) {
                printf ("Error: Precision matrix is not positive semi-definite.");
                exit (EXIT_FAILURE);
            }
            logdet += dim * log(hscale);
            for (i = 0; i < dim; ++i) {
                for (j = 0; j < dim; ++j) {
                    intlam_inv_flat[sqindx + (i * dim + j)] = amat[i][j+dim] / hscale;
                }
            }
            for (i = 0; i < dim; ++i) {
                delete [] amat[i];
            }

            /* V_z | Means of the multivariate Gaussian obtained from integration
               The inner part is a vector prec.dot.v + mu / sig2 - mureg / reg2
               Note that vzinner is an array with zmax elements, so positions other than i can hold spurious values
            */
            for (i = 0; i < dim; ++i) {
                ki = zarr[zindx + i];
                vzinner[i] = precdotv[ki];
                vzinner[i] += mu[ki] / sig2[ki];
                vzinner[i] -= mureg / reg2;
            }
            for (i = 0; i < dim; ++i) {
                ai = 0.0;
                for (j = 0; j < dim; ++j) {
                    ai += intlam_inv_flat[sqindx + (i * dim + j)] * vzinner[j];
                }
                intmu[zindx + i] = ai;
            }

            /* Log of F(z, theta)
               We calculate the prefactor and terms separately and add them up
            */
            fzterm1 = 0.0;
            fzterm2 = logdet;
            fzterm3 = 0.0;
            fzterm4 = 0.0;
            for (i = 0; i < dim; ++i) {
                ki = zarr[zindx + i];
                fzterm1 += log (1 / sig2[ki]);
                fztmp = log (1 / sig2[ki]);
                if (fztmp == INFINITY) {
                    printf ("%i: %g\n", ki, sig2[ki]);
                }
                fzterm3 += mu[ki] * mu[ki] / sig2[ki];
                fzterm4 += intmu[zindx + i] * vzinner[i];
            }
        
            logfz = 0.5 * (fzterm1 - fzterm2 - fzterm3 + fzterm4);
            zcomps[z] = exp(logpz + logfz);
            if (zcomps[z] != zcomps[z]) {
                if (! is_covariate) {
                    printf ("Not a covariate locus\n");
                }
                printf ("Error in zstate %i: %g %g %g\n", z, zcomps[z], logpz, logfz);
                exit(0);
            }

            // update the indices
            zindx += dim;
            sqindx += dim * dim;

        } // end if for znorm = 0
    } // end loop over all zstates

    //printf ("End of z-loop %i\n", zindx);

    if (get_gradient) {
        //printf ("Calculating gradients\n");
        compute_grad(nvar, nhyp, nfeat, zlen, pi, mu, sig2, feat_flat, zarr, znorm, zcomps, intmu, intlam_inv_flat, is_covariate, grad);
    }

    delete [] intmu;
    delete [] intlam_inv_flat;

    return;
}

void compute_grad(int nvar,
                 int nhyp,
                 int nfeat,
                 int zlen,
                 double *pi,
                 double *mu,
                 double *sig2,
                 double *feat_flat,
                 int *zarr,
                 int *znorm,
                 double *zcomps,
                 double *intmu,
                 double *intlam_inv_flat,
                 bool is_covariate,
                 double *grad
                ) {

    int i;
    int k;
    int ki;
    int z;
    int dim;
    int zindx, sqindx;
    double zcompsum;
    double picompsum;
    double grad_pi, innersumpi;
    double grad_mu, innersummu;
    double grad_sig, innersumsig;
    double feat;
    double mudiff, mudiff2, lamdiag;
    double pz[zlen];

    //printf ("Entering gradient calculation\n");

    // Initiate to zero for this locus
    for (i = 0; i < nhyp; ++i) {
        grad[i] = 0.0;
    }

    zcompsum = 0.0;
    for (i = 0; i < zlen; ++i) {
        zcompsum += zcomps[i];
    }

    for (i = 0; i < zlen; ++i) {
        pz[i] = zcomps[i] / zcompsum;
    }

    if (!is_covariate) {
        // Loop over all features for the gradients of pi
        for (k = 0; k < nfeat; ++k) {
    
            // For each feature row, precalculate (- pi .dot. feat[k])
            picompsum = 0.0;
            for (i = 0; i < nvar; ++i) {
                picompsum += - pi[i] * feat_flat[k * nvar + i];
            }
    
            grad_pi = 0.0;
            zindx = 0;
     
            for (z = 0; z < zlen; ++z) {
                innersumpi  = picompsum;
                dim = znorm[z];
                for (i = 0; i < dim; ++i) {
                    ki = zarr[zindx + i];
                    feat = feat_flat[k * nvar + ki]; //feat[k][ki];
                    innersumpi  += feat;
                }
                grad_pi  += pz[z] * innersumpi;
                zindx += dim;
            }
            grad[k] = grad_pi;
        }
    }

    // Gradient of mu and sigma

    grad_mu  = 0.0;
    grad_sig = 0.0;
    zindx = 0;
    sqindx = 0;

    for (z = 0; z < zlen; ++z) {

        innersummu  = 0.0;
        innersumsig = 0.0;

        dim = znorm[z];
        for (i = 0; i < dim; ++i) {
            ki = zarr[zindx + i];
            mudiff = intmu[zindx + i] - mu[ki];
            mudiff2 = mudiff * mudiff;
            lamdiag = intlam_inv_flat[sqindx + (i * dim + i)];

            innersummu  += mudiff / sig2[ki];
            innersumsig += (mudiff2 - sig2[ki] + lamdiag) / sig2[ki];
            //innersumsig += 2.0 * (mudiff2 - sig2[ki] + lamdiag) / (sig2[ki] * sqrt(sig2[ki]));
        }

        grad_mu  += pz[z] * innersummu;
        grad_sig += pz[z] * innersumsig;

        zindx += dim;
        sqindx += dim * dim;
    }

    if (!is_covariate) {
        grad[nhyp - 3] = grad_mu;
        grad[nhyp - 2] = 0.5 * grad_sig;
    } else {
        grad[nhyp - 1] = 0.5 * grad_sig;
    }
    //printf ("grad = %g\n", grad_mu);

    return;
}


DLLEXPORT void freeme(int *ptr) {
    printf("freeing address: %p\n", ptr);
    free(ptr);
}


// for pointer of pointers, see: http://stackoverflow.com/questions/22425921/pass-a-2d-numpy-array-to-c-using-ctypes
