#include <stdio.h>

#include "problem.h"
#include "usefull.h"
#include "poisson.h"

int main(int argc, char *argv[]){
    char *fileProblem = "/mnt/d/cours/Q8/LMECA2660/Project/Data/problem.txt";
    // Poisson_data poisson_data;
    // PetscInitialize(&argc, &argv, 0, 0);
    // initialize_poisson_solver(&poisson_data);


    double h_hill = 500.0;
    double u_hill = 50.0;
    double y0 = 0.1;
    double t_max = 3.0;
    double dt = 0.01;
    int n_iter = (int)(t_max / dt);
    // double Re = 100000000;


    Problem *theProblem;

    theProblem = initProblem(h_hill, u_hill, y0);
    Mesh *U = theProblem->u;
    Mesh *V = theProblem->v;
    Mesh *U_star = theProblem->u_star;
    Mesh *V_star = theProblem->v_star;
    Mesh *Hx = theProblem->Hx;
    Mesh *Hy = theProblem->Hy;
    Mesh *Hx_old = theProblem->Hx_old;
    Mesh *Hy_old = theProblem->Hy_old;
    Mesh *Divx = theProblem->divx;
    Mesh *Divy = theProblem->divy;
    Mesh *Grad_px = theProblem->grad_px;
    Mesh *Grad_py = theProblem->grad_py;


    double **u = U->grid;
    double **v = V->grid;
    double **u_star = U_star->grid;
    double **v_star = V_star->grid;
    double **hx = Hx->grid;
    double **hy = Hy->grid;
    double **hx_old = Hx_old->grid;
    double **hy_old = Hy_old->grid;
    double **divx = Divx->grid;
    double **divy = Divy->grid;
    double **grad_px = Grad_px->grid;
    double **grad_py = Grad_py->grid;

    saveProblem(theProblem, fileProblem);
    // printf("%f\n", theProblem->u_tau);
    initialCondition(theProblem);
    boundaryConditions(theProblem);

    // First iteration using a EE2 scheme for the advective term
    advective(theProblem);
    diffusive(theProblem);
    gradP(theProblem);


    for(int j = 0; j < Grad_px->Ny; j ++){  // The reference indices are taken for grad_p
        for(int i = 0; i < Grad_px->Nx; i ++){
            u_star[j+1][i+1] = u[j+1][i+1] + dt * (- hx[j+1][i+1] + 2.0 * divx[j][i] - grad_px[j][i]);
            hx_old[j+1][i+1] = hx[j+1][i+1];
        }
    }

    for(int j = 0; j < Grad_py->Ny; j ++){
        for(int i = 0; i < Grad_px->Nx; i ++){
            v_star[j+1][i+1] = v[j+1][i+1] + dt * (- hx[j+1][i+1] + 2.0 * divy[j][i] - grad_py[j][i]);
            hy_old[j+1][i+1] = hy[j+1][i+1];
        }
    }


    vorticity(theProblem);

    printf("Size of u: (%d, %d)\n", theProblem->u->Nx, theProblem->u->Ny);
    printf("Size of v: (%d, %d)\n", theProblem->v->Nx, theProblem->v->Ny);
    printf("Size of nu: (%d, %d)\n", theProblem->nu->Nx, theProblem->nu->Ny);
    printf("Number of iterations %d\n", n_iter);

    // Time loop

    for(int it = 0; it < n_iter; it++){

    }

    // printMat(theProblem->u->Nx, theProblem->u->Ny, theProblem->u->grid);

    saveMat(theProblem->u->Nx, theProblem->u->Ny, theProblem->u->grid, "u");
    saveMat(theProblem->v->Nx, theProblem->v->Ny, theProblem->v->grid, "v");
    saveMat(theProblem->Hx->Nx, theProblem->Hx->Ny, theProblem->Hx->grid, "Hx");
    saveMat(theProblem->Hy->Nx, theProblem->Hy->Ny, theProblem->Hy->grid, "Hy");
    saveMat(theProblem->nu->Nx, theProblem->nu->Ny, theProblem->nu->grid, "w");
    freeProblem(theProblem);

    // PetscFinalize();

}
