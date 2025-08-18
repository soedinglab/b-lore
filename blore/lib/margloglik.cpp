/*
 * =====================================================================================
 *
 *       Filename:  margloglik.c
 *
 *    Description:  Calculate the log marginal likelihood
 *
 *        Version:  1.0
 *        Created:  25.07.2017 10:35:03
 *       Revision:  15.08.2025 | replace gauss_jordan_inversion with Eigen
 *                  18.08.2025 | replace exit() calls with error codes
 *       Compiler:  gcc
 *
 *         Author:  Saikat Banerjee (banskt), bnrj.saikat@gmail.com
 *   Organization:  Max Planck Institute for Biophysical Chemistry
 *
 * =====================================================================================
 */

#define DLLEXPORT extern "C"
#define ERR_SUCCESS              0
#define ERR_INVALID_PI           1
#define ERR_MATRIX_NOT_PSD       2
#define ERR_DECOMPOSITION_FAILED 3

#include <stdio.h>      /* printf */
#include <math.h>       /* log, exp */
#include <stdlib.h>     /* abs */
#include <iostream>
#include "Eigen/Dense"

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  print_matrix
 *  Description:  print a matrix on cout for debugging
 *                
 * =====================================================================================
 */

void print_matrix (double **matrix, int rows, int columns) {
    for (int i=0; i<rows; ++i){
        for (int j=0; j<columns; ++j) {
            cout << matrix[i][j] << ' ';
        }
        cout << '\n';
    }
    cout << '\n';
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  gauss_jordan_inversion
 *  Description:  calculate matrix inverse using Gauss-Jordan elimination method
 *                fast for small matrices
 *                
 * =====================================================================================
 */
    void
gauss_jordan_inversion ( int n, double **matrix, double &logdet, double &signdet )
{
    int i, j, k;
    double ratio;
    double a;

    const double eps = 1e-12;  // tolerance for zero pivot
    int pivot_row;
    double max_val, tmp, factor, diag;

//  hardcoded for fast calculation of single elements
    if (n == 1) {
        a = matrix[0][0];
        matrix[0][1] = 1 / a;
        signdet = 1;
        if (a < 0) {
            logdet = log(-a);
            signdet = -1;
        } else {
            logdet = log(a);
            signdet = 1;
        }

    } else {

        // Create the augmented matrix: identity matrix beside the original one
        for(i = 0; i < n; ++i){
            for(j = n; j < 2*n; ++j){
                if(i==(j-n)) {
                    matrix[i][j] = 1.0;
                } else {
                    matrix[i][j] = 0.0;
                }
            }
        }

        for (i = 0; i < n; ++i) {
            // Pivoting: find row with largest absolute value in column i
            pivot_row = i;
            max_val = fabs(matrix[i][i]);
            for (j = i + 1; j < n; ++j) {
                if (fabs(matrix[j][i]) > max_val) {
                    max_val = fabs(matrix[j][i]);
                    pivot_row = j;
                }
            }

            // Check for numerical stability
            if (max_val < eps) {
                signdet = 0;
                logdet = -INFINITY;
                printf("Error: Matrix is singular or ill-conditioned (pivot < %.1e).\n", eps);
                return;
            }

            // Swap rows if necessary
            if (pivot_row != i) {
                signdet *= -1;
                for (k = 0; k < 2 * n; ++k) {
                    tmp = matrix[i][k];
                    matrix[i][k] = matrix[pivot_row][k];
                    matrix[pivot_row][k] = tmp;
                }
            }

            // Elimination: make all rows except pivot zero in column i
            diag = matrix[i][i];
            for (j = 0; j < n; ++j) {
                if (j != i) {
                    factor = matrix[j][i] / diag;
                    for (k = 0; k < 2 * n; ++k) {
                            matrix[j][k] -= factor * matrix[i][k];
                    }
                }
            }

            // Normalize the pivot row
            for (k = 0; k < 2 * n; ++k) {
                matrix[i][k] /= diag;
            }

            // Determinant contribution
            logdet += log(fabs(diag));
            if (diag < 0) {
                signdet *= -1;
            }
        }
    }
}		/* -----  end of function gauss_jordan_inversion  ----- */



/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  compute_zcomps
 *  Description:  compute (logPz + logNz) for every z-state, and store them in zcomps
 *                also compute mu*_z and inverse(lam*_z) for every zstate
 * =====================================================================================
 */
    int
compute_zcomps ( int nvar, int zlen, int zmax,
                 int *zarr, int *znorm, 
                 double *pi, double *mu, double *sig2,
                 double *v, double *precflat,
                 double *intmu, double *intlam_inv_flat, double *zcomps,
                 double reg2, double mureg,
                 bool is_covariate
               )
{
    int i, j, z;
    int zindx, sqindx;
    int dim;
    int ki, kj;
    double sumlogpi;
    double logdet, signdet;
    double logpz, logfz;

    
    double precdotv[nvar];
    double intlam[nvar][nvar];
//    double *amat[zmax];
    double vzinner[zmax];

    double hscale;
    double ai, aij;
    double fzterm1, fzterm2, fzterm3, fzterm4;

//  Calculate the sum of log (1 - pi) only once
//  Calculate Lambda_z only once, but note: this will give speed benefit only when zlen > n^2, which is possible if znorm > 2
//  We also precalculate prec.dot.v 
    sumlogpi = 0.0;
    if (!is_covariate) {
        for (i = 0; i < nvar; ++i) {
            if (pi[i] < 1.0) {
                sumlogpi += log(1 - pi[i]);
            } else {
                std::cerr << "Pi value reached 1. Aborting...\n";
                return ERR_INVALID_PI;
            }
        }
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

//  Start the loop over all zstates for the calculation of marginal likelihood
    zindx = 0;
    sqindx = 0;
    for (z = 0; z < zlen; ++z) {
//      if znorm is 0, then we need to calculate only zcomps. 
//      intmu and intlam are empty, hence not updated.
//      Note that zindx and sqindx also does not change.
        dim = znorm[z];
        if (dim == 0) {
            zcomps[z] = sumlogpi + 1.0;

        } else {

//          Log of probability of z given Xi and Theta => ln p(z | Xi, Theta)
            logpz = sumlogpi;
            if (!is_covariate) {
                for (i = 0; i < dim; ++i) {
                    ki = zarr[zindx + i];
                    logpz += log (pi[ki] / (1 - pi[ki]));
                }
            }

//          Select Lambda_z and invert
            hscale = 1 / sig2[zarr[zindx]];
            MatrixXd amat(dim, dim);
            for (i = 0; i < dim; ++i) {
                ki = zarr[zindx + i];
                for (j = 0; j < dim; ++j) {
                    kj = zarr[zindx + j];
                    amat(i,j) = intlam[ki][kj] / hscale;
                }
            }

//          Attempt Cholesky decomposition
            Eigen::LLT<MatrixXd> llt(amat);
            MatrixXd amat_inv;
            logdet = 0.0;
            signdet = 1.0;

            if (llt.info() == Eigen::Success) {
                amat_inv = llt.solve(MatrixXd::Identity(dim, dim));
                double det = amat.determinant();
                logdet = std::log(std::abs(det));
                if (det < 0) signdet = -1.0;
            } else {

//              Fallback to LDLT 
                Eigen::LDLT<MatrixXd> ldlt(amat);
                if (ldlt.info() != Eigen::Success) {
                    if (!ldlt.isPositive()) {
                        std::cerr << "Error: Precision matrix is not positive semi-definite.\n";
                        return ERR_MATRIX_NOT_PSD;
                    } else {
                        std::cerr << "Error: LDLT decomposition failed (numerical instability).\n";
                        return ERR_DECOMPOSITION_FAILED;
                    }
                } else {
                    amat_inv = ldlt.solve(MatrixXd::Identity(dim, dim));
                    double det = amat.determinant();
                    logdet = std::log(std::abs(det));
                    if (det < 0) signdet = -1.0;
                }
            }

/*
            //gauss_jordan_inversion(dim, amat, logdet, signdet, is_covariate, sig2);
            gauss_jordan_inversion(dim, amat, logdet, signdet);
            if (signdet < 0) {
                printf ("Error: Precision matrix is not positive semi-definite.\n");
                printf ("Info: Value of signdet is %f\n", signdet);
                printf ("Info: Value of logdet is %f\n", logdet);
                printf ("Info: Value of dim is %d\n", dim);
                printf ("Info: Value of hscale is %f\n", hscale);
                printf ("Info: Value of sig2 is %f\n", 1 / hscale);
                printf ("Info: Value of reg2 is %f\n", reg2);
                print_matrix (amat, dim, 2*dim);
                //Debug | get input matrix
                for (i = 0; i < dim; ++i) {
                    ki = zarr[zindx + i];
                    for (j = 0; j<dim; ++j) {
                        kj = zarr[zindx + j];
                        amat[i][j] = intlam[ki][kj] / hscale;
                    }
                }
                print_matrix (amat, dim, 2*dim);
                exit (EXIT_FAILURE);
            }
*/
            logdet += dim * log(hscale);
            for (i = 0; i < dim; ++i) {
                for (j = 0; j < dim; ++j) {
                    intlam_inv_flat[sqindx + (i * dim + j)] = amat_inv(i, j) / hscale;
                }
            }
            /*
            for (i = 0; i < dim; ++i) {
                delete [] amat[i];
            }
            */

//          V_z | Means of the multivariate Gaussian obtained from integration
//          The inner part is a vector prec.dot.v + mu / sig2 - mureg / reg2
//          Note that vzinner is an array with zmax elements, so positions other than i can hold spurious values
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

//          Log of F(z, theta)
//          We calculate the prefactor and terms separately and add them up
            fzterm1 = 0.0;
            fzterm2 = logdet;
            fzterm3 = 0.0;
            fzterm4 = 0.0;
            for (i = 0; i < dim; ++i) {
                ki = zarr[zindx + i];
                fzterm1 += log (1 / sig2[ki]);
                fzterm3 += mu[ki] * mu[ki] / sig2[ki];
                fzterm4 += intmu[zindx + i] * vzinner[i];
            }

            logfz = 0.5 * (fzterm1 - fzterm2 - fzterm3 + fzterm4);
            zcomps[z] = logpz + logfz;

//          update the indices
            zindx += dim;
            sqindx += dim * dim;

        } // end if (dim != 0)
    } // end loop over all zstates
    return ERR_SUCCESS;

}		/* -----  end of function compute_zcomps  ----- */



/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  compute_grad
 *  Description:  
 * =====================================================================================
 */
    void
compute_grad ( int nvar, int nhyp, int nfeat, int zlen,
               double *pi, double *mu, double *sig2, double *feat_flat,
               int *zarr, int *znorm, double *pz,
               double *intmu, double *intlam_inv_flat,
               bool is_covariate, double *grad
             )
{

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

//  Initiate to zero for this locus
    for (i = 0; i < nhyp; ++i) {
        grad[i] = 0.0;
    }

    if (!is_covariate) {
//      Loop over all features for the gradients of pi
        for (k = 0; k < nfeat; ++k) {

//          For each feature row, precalculate (- pi .dot. feat[k])
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
                    feat = feat_flat[k * nvar + ki];
                    innersumpi  += feat;
                }
                grad_pi  += pz[z] * innersumpi;
                zindx += dim;
            }
            grad[k] = grad_pi;
        }
    }

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

    //printf ("%f \n", 0.5 * grad_sig);

    //for (i = 0; i < nhyp; ++i) {
    //    printf ("%f ", grad[i]);
    //}

    //printf ("\n");

}		/* -----  end of function compute_grad  ----- */


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  margloglik_zcomps
 *  Description:  Wrapper function for log marginal likelihood and its gradient
 *                This is the interface with python
 * =====================================================================================
 */
DLLEXPORT double
margloglik_zcomps ( int nvar,                   // number of SNPs
                    int nhyp,                   // number of hyperparameters (required only for gradient)
                    int nfeat,                  // number of features (required only for gradient)
                    int zlen,                   // number of zstates
                    int zmax,                   // max norm of the zstates
                    double *pi,
                    double *mu,
                    double *sig2,
                    double *v,
                    double *precflat,           // flattened array of precision matrix
                    double *feat_flat,          // flattened array of features (required only for gradient)
                    int *zarr,                  // array containing the causal picks of all zstates
                    int *znorm,                 // array with the number of causal picks of all zstates
                    double reg2,
                    double mureg,
                    double *zcomps,             // output
                    double *grad,               // output
                    bool get_gradient,
                    bool is_covariate,
                    double *logmL               // output
                  )
{
    int i, j, k, ki, kj;
    int ai;
    int retcode;
    double logk, zcompsum;

    double *intmu;
    double *intlam_inv_flat;

//  How many elements will there be in intmu and intlam_inv?
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

    retcode = compute_zcomps(nvar, zlen, zmax, zarr, znorm, pi, mu, sig2, v, precflat, intmu, intlam_inv_flat, zcomps, reg2, mureg, is_covariate);

    if (retcode != ERR_SUCCESS) {
        *logmL = NAN;
        return retcode;
    }

    logk = zcomps[0];
    for ( i = 1; i < zlen; i++ ) {
        if ( zcomps[i] > logk ) {
            logk = zcomps[i];
        }
    }

    zcompsum = 0.0;
    for (i = 0; i < zlen; i++) {
        zcompsum += exp(zcomps[i] - logk);
    }
    *logmL = log(zcompsum) + logk ;

    for ( i = 0; i < zlen; i++ ) {
        zcomps[i] = exp(zcomps[i] - *logmL);
    }

    if ( get_gradient ) {
        compute_grad(nvar, nhyp, nfeat, zlen, pi, mu, sig2, feat_flat, zarr, znorm, zcomps, intmu, intlam_inv_flat, is_covariate, grad);
    }

    delete [] intmu;
    delete [] intlam_inv_flat;

    return retcode;
}		/* -----  end of function margloglik_zcomps  ----- */
