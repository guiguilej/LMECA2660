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
    printf("%f\n", theProblem->u_tau);
    initialCondition(theProblem);
    boundaryConditions(theProblem);
    saveMat(theProblem->u->Nx, theProblem->u->Ny, theProblem->u->grid, "u");
    saveMat(theProblem->v->Nx, theProblem->v->Ny, theProblem->v->grid, "v");
    freeProblem(theProblem);

    // PetscFinalize();

}
