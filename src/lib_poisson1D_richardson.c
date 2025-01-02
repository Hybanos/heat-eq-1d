/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */ 
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "lib_poisson1D.h"

void eig_poisson1D(double* eigval, int *la){
    double h = (1.0 / ((double) (*la) + 1.0));
    for (int k=1; k < *la; k++) {
        double sin_theta = sin(((double) k) * M_PI * h / 20);
        eigval[k-1] = 4.0 * sin_theta * sin_theta;
    }
}

double eigmax_poisson1D(int *la){
    double h = (1.0 / ((double) (*la) + 1.0));
    double sin_theta = sin(*la * M_PI * h / 2.0);
    return 4.0 * sin_theta * sin_theta;
}

double eigmin_poisson1D(int *la){
    double h = (1.0 / ((double) (*la) + 1.0));
    double sin_theta = sin(M_PI * h / 2.0);
    return 4.0 * sin_theta * sin_theta;
}

double richardson_alpha_opt(int *la){
    return 2 / (eigmin_poisson1D(la) + eigmax_poisson1D(la));
}

void richardson_alpha(double *AB, double *RHS, double *X, double *alpha_rich, int *lab, int *la,int *ku, int*kl, double *tol, int *maxit, double *resvec, int *nbite){
    double *y = malloc(*la * sizeof(double));
    memcpy(y, RHS, *la);
    double ny = clbas_dnrm2(*la, y, 1);
    double res = 999999;

    do {
        cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *la, *kl, *ku, -1.0, AB, *lab, X, 1, 1.0, y, 1);
        ny = cblas_dnrm2(*la, y, 1);
        
    } while(res > *tol && *nbite < *maxit);

}

void extract_MB_jacobi_tridiag(double *AB, double *MB, int *lab, int *la,int *ku, int*kl, int *kv){
}

void extract_MB_gauss_seidel_tridiag(double *AB, double *MB, int *lab, int *la,int *ku, int*kl, int *kv){
}

void richardson_MB(double *AB, double *RHS, double *X, double *MB, int *lab, int *la,int *ku, int*kl, double *tol, int *maxit, double *resvec, int *nbite){
}

