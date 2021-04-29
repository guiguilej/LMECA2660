#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "problem.h"
#include "usefull.h"


Problem *initProblem(double h_hill, double u_hill, double y0){
    Problem *theProblem = malloc(sizeof(Problem));

    theProblem->h_hill = h_hill;
    theProblem->H = 4.0 * h_hill;
    theProblem->L = 20.0 * h_hill;
    theProblem->d_hill = 8.0 * h_hill;
    theProblem->sigma_hill = 3.0 / 2.0 * h_hill;
    theProblem->u_hill = u_hill / 3.6;                                              // [m/s]

    theProblem->y0 = y0;
    theProblem->kappa = 0.41;
    theProblem->u_tau = theProblem->kappa * u_hill / log((h_hill + y0) / y0) / 3.6; // [m/s]
    theProblem->C = 0.5;
    theProblem->h = 500.0;                                                          // To change later
    theProblem->Nx = (int) theProblem->L / theProblem->h;
    theProblem->Ny = (int) theProblem->H / theProblem->h;

    // For u and v, the first point of the mesh corresponds to the u and v which  get in the first cell
    // The addidtional points (of index 0) are thus the v's which lie on the bottom boundary and the u's
    // which line on the left hand side boundary
    theProblem->u = initMesh(theProblem->Nx+1, theProblem->Ny+2, theProblem->h);
    theProblem->v = initMesh(theProblem->Nx+2, theProblem->Ny+1, theProblem->h);
    theProblem->u_p = initU_p(theProblem);
    theProblem->p = initMesh(theProblem->Nx, theProblem->Ny, theProblem->h);
    theProblem->d = initMesh(theProblem->Nx+1, theProblem->Ny+1, theProblem->h);
    theProblem->d1 = initMesh(theProblem->Nx+1, theProblem->Ny+1, theProblem->h);
    theProblem->d2 = initMesh(theProblem->Nx+1, theProblem->Ny+1, theProblem->h);
    theProblem->d3 = initMesh(theProblem->Nx+1, theProblem->Ny+1, theProblem->h);
    theProblem->d4 = initMesh(theProblem->Nx+1, theProblem->Ny+1, theProblem->h);
    theProblem->nu = initMesh(theProblem->Nx+1, theProblem->Ny+1, theProblem->h);
    theProblem->w = initMesh(theProblem->Nx+1, theProblem->Ny+1, theProblem->h);

    return theProblem;

}

double *initU_p(Problem *theProblem){
    double *u_p = calloc(theProblem->u->Ny, sizeof(double));
    double h = theProblem->h;
    printf("%d\n", theProblem->u->Ny);
    printf("%d\n", theProblem->u->Nx);
    for(int j = 1; j < theProblem->u->Ny - 1; j++){
            double y = h / 2.0 + (j - 1) * h;
            //printf("%f\n", y);
            u_p[j] = theProblem->u_tau / theProblem->kappa * log((y + theProblem->y0) / theProblem->y0);
        }
    return u_p;
}

Mesh *initMesh(int Nx, int Ny, double h){
    Mesh *mesh = malloc(sizeof(Mesh));

    mesh->Nx = Nx;
    mesh->Ny = Ny;
    mesh->h = h;
    mesh->grid = calloc(Ny, sizeof(double*));
    for(int j = 0; j < Ny; j++){
        mesh->grid[j] = calloc(Nx, sizeof(double*));
    }
    return mesh;
}

void boundaryConditions(Problem *theProblem){
    Mesh *U = theProblem->u;
    Mesh *V = theProblem->v;

    double **u = U->grid;
    double **v = V->grid;

    // Left conditions:
    // u is already computed in the function initialCondition
    for(int j = 1; j < V->Ny; j++){
        v[j][0] = - 1.0 / 5.0 * (v[j][3] - 5.0 * v[j][2] + 15.0 * v[j][1]);
    }

    // Top conditions:
    for(int i = 0; i < V->Nx; i++){
        v[V->Ny-1][i] = 0.0;
    }
    for(int i = 0; i < U->Nx; i++){
        u[U->Ny-1][i] = u[U->Ny-2][i];
    }

    // Bottom conditions:
    for(int i = 0; i < U->Nx; i++){
        u[0][i] = - 1.0 / 5.0 * (u[3][i] - 5.0 * u[2][i] + 15.0 * u[1][i]);
    }
    for(int i = 0; i < V->Nx; i++){
        v[0][i] = 0.0;
    }

    // Right conditions:
    // u conditions still need to be added;
    for(int j = 0; j < V->Ny; j++){
        v[j][V->Nx-1] = v[j][V->Nx-2];
    }
}

void diffusiveTerm(Problem *theProblem){
    Mesh *U = theProblem->u;
    Mesh *V = theProblem->v;
    Mesh *D = theProblem->d;
    Mesh *D1 = theProblem->d1;
    Mesh *D2 = theProblem->d2;
    Mesh *D3 = theProblem->d3;
    Mesh *D4 = theProblem->d4;
    Mesh *Nu = theProblem->nu;

    double **u = U->grid;
    double **v = V->grid;
    double **d = D->grid;
    double **d1 = D1->grid;
    double **d2 = D2->grid;
    double **d3 = D3->grid;
    double **d4 = D4->grid;
    double **nu = Nu->grid;

    double h = theProblem->h;

    for(int j = 1; j < D->Ny; j++){
        for(int i = 1; i < D->Nx; i++){
            
        }
    }
}

void freeMesh(Mesh *mesh){
    for(int j = 0; j < mesh->Ny; j++){
        free(mesh->grid[j]);
    }
    free(mesh->grid);
    free(mesh);
}

void freeProblem(Problem *theProblem){
    free(theProblem->u_p);
    freeMesh(theProblem->u);
    freeMesh(theProblem->v);
    freeMesh(theProblem->p);
    freeMesh(theProblem->d);
    freeMesh(theProblem->d1);
    freeMesh(theProblem->d2);
    freeMesh(theProblem->d3);
    freeMesh(theProblem->d4);
    freeMesh(theProblem->nu);
    freeMesh(theProblem->w);
    free(theProblem);
}


void initialCondition(Problem *theProblem){
    for(int j = 1; j < theProblem->u->Ny - 1; j++){
        for(int i = 0; i < theProblem->u->Nx; i++){
            theProblem->u->grid[j][i] = theProblem->u_p[j];
        }
    }
}
