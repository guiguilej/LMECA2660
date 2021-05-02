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
    
    Problem *theProblem;

    theProblem = initProblem(h_hill, u_hill, y0);
    saveProblem(theProblem, fileProblem);
    // printf("%f\n", theProblem->u_tau);
    initialCondition(theProblem);
    boundaryConditions(theProblem);
    advective(theProblem);
    diffusive(theProblem);
    gradP(theProblem);

    printf("Size of u: (%d, %d)\n", theProblem->u->Nx, theProblem->u->Ny);
    printf("Size of v: (%d, %d)\n", theProblem->v->Nx, theProblem->v->Ny);
    printf("Size of nu: (%d, %d)\n", theProblem->nu->Nx, theProblem->nu->Ny);

    // printMat(theProblem->u->Nx, theProblem->u->Ny, theProblem->u->grid);

    saveMat(theProblem->u->Nx, theProblem->u->Ny, theProblem->u->grid, "u");
    saveMat(theProblem->v->Nx, theProblem->v->Ny, theProblem->v->grid, "v");
    saveMat(theProblem->Hx->Nx, theProblem->Hx->Ny, theProblem->Hx->grid, "Hx");
    saveMat(theProblem->Hy->Nx, theProblem->Hy->Ny, theProblem->Hy->grid, "Hy");
    freeProblem(theProblem);

    // PetscFinalize();

}
