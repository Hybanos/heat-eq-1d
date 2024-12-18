/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */ 
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "lib_poisson1D.h"

void set_GB_operator_colMajor_poisson1D(double* AB, int *lab, int *la, int *kv) {
    if (*kv) {
        for (size_t i = 0; i < *la; i++) {
            AB[0 + i*4] =  0;
            AB[1 + i*4] = -1;
            AB[2 + i*4] =  2;
            AB[3 + i*4] = -1;
        }
    } else {
        for (size_t i = 0; i < *la; i++) {
            AB[0 + i*3] = -1;
            AB[1 + i*3] =  2;
            AB[2 + i*3] = -1;
        }
    }
    AB[*kv] = 0;
    AB[*la*((*lab-1) + *kv)-1] = 0;
}

void set_GB_operator_colMajor_poisson1D_Id(double* AB, int *lab, int *la, int *kv){
    for (size_t i = 0; i < *la; i++) {
        AB[1 + *kv] = 1;
    }
}

void set_dense_RHS_DBC_1D(double* RHS, int* la, double* BC0, double* BC1){

}

void set_analytical_solution_DBC_1D(double* EX_SOL, double* X, int* la, double* BC0, double* BC1){
}  

void set_grid_points_1D(double* x, int* la){
}

double relative_forward_error(double* x, double* y, int* la){
}

int indexABCol(int i, int j, int *lab){
    // printf("i: %d, j: %d, lab: %d\n", i, j, *lab);
    return (j * (*lab) + i);
}

int dgbtrftridiag(int *la, int*n, int *kl, int *ku, double *AB, int *lab, int *ipiv, int *info){
    
    double p = 0.;
    *info = 0;

    if (!*kl || !*ku) {
        *info = -1;
        return *info;
    }

    for (int j = 0; j < *n - 1; j++) {
        p = AB[indexABCol(0, j+1, lab)] / AB[indexABCol(1, j, lab)];
        AB[indexABCol(0, j+1, lab)] = p;
        AB[indexABCol(1, j+1, lab)] -= p * AB[indexABCol(2, j, lab)];
    }

    for (int j = 0; j < *n; j++) {
        ipiv[j] = j+1;
    }

    return *info;
}
