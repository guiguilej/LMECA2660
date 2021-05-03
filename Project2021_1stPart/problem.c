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
    theProblem->d_hill = 4.0 * h_hill;
    theProblem->sigma_hill = 3.0 / 2.0 * h_hill;
    theProblem->u_hill = u_hill / 3.6;                                              // [m/s]

    theProblem->y0 = y0;
    theProblem->kappa = 0.41;
    theProblem->u_tau = theProblem->kappa * u_hill / log((h_hill + y0) / y0) / 3.6; // [m/s]
    theProblem->C = 0.5;

    double h = 500;
    // double h = h_hill / 40.0;
    theProblem->h = h;                                                          // To change later
    int Nx = (int) theProblem->L / theProblem->h;
    int Ny = (int) theProblem->H / theProblem->h;
    theProblem->Nx = Nx;
    theProblem->Ny = Ny;

    // For u and v, the first point of the mesh corresponds to the u and v which  get in the first cell
    // The addidtional points (of index 0) are thus the v's which lie on the bottom boundary and the u's
    // which lie on the left hand side boundary
    theProblem->u = initMesh(Nx+1, Ny+2, h);
    theProblem->v = initMesh(Nx+2, Ny+1, h);
    theProblem->u_star = initMesh(Nx+1, Ny+2, h);
    theProblem->v_star = initMesh(Nx+2, Ny+1, h);
    theProblem->u_p = initU_p(theProblem);
    theProblem->p = initMesh(Nx, Ny, h);
    theProblem->phi = initMesh(Nx, Ny, h);
    theProblem->Hx = initMesh(Nx+1, Ny+2, h);
    theProblem->Hy = initMesh(Nx+2, Ny+1, h);
    theProblem->Hx_old = initMesh(Nx+1, Ny+2, h);
    theProblem->Hy_old = initMesh(Nx+2, Ny+1, h);
    theProblem->d1 = initMesh(Nx+1, Ny+1, h);
    theProblem->d2 = initMesh(Nx+1, Ny+1, h);
    theProblem->d3 = initMesh(Nx+1, Ny+1, h);
    theProblem->nu = initMesh(Nx+1, Ny+1, h);
    theProblem->grad_px = initMesh(Nx-1, Ny, h);
    theProblem->grad_py = initMesh(Nx, Ny-1, h);
    theProblem->grad_phix = initMesh(Nx-1, Ny, h);
    theProblem->grad_phiy = initMesh(Nx, Ny-1, h);
    theProblem->divx = initMesh(Nx+1, Ny+1, h);
    theProblem->divy = initMesh(Nx+1, Ny+1, h);
    theProblem->w = initMesh(Nx+1, Ny+1, h);

    return theProblem;

}

double *initU_p(Problem *theProblem){
    double *u_p = calloc(theProblem->u->Ny, sizeof(double));
    double h = theProblem->h;
    for(int j = 1; j < theProblem->u->Ny - 1; j++){
        double y = h / 2.0 + (j - 1) * h;
        u_p[j] = theProblem->u_tau / theProblem->kappa * log((y + theProblem->y0) / theProblem->y0);
        theProblem->massFlow += u_p[j] * h;
    }
    return u_p;
}

Mesh *initMesh(int Nx, int Ny, double h){
    Mesh *mesh = malloc(sizeof(Mesh));

    mesh->Nx = Nx;
    mesh->Ny = Ny;
    mesh->h = h;
    mesh->grid = (double **) malloc(Ny * sizeof(double*));
    for(int j = 0; j < Ny; j++){
        mesh->grid[j] = (double *) calloc(Nx, sizeof(double));
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
    freeMesh(theProblem->u_star);
    freeMesh(theProblem->v_star);
    freeMesh(theProblem->p);
    freeMesh(theProblem->phi);
    freeMesh(theProblem->Hx);
    freeMesh(theProblem->Hy);
    freeMesh(theProblem->Hx_old);
    freeMesh(theProblem->Hy_old);
    freeMesh(theProblem->d1);
    freeMesh(theProblem->d2);
    freeMesh(theProblem->d3);
    freeMesh(theProblem->nu);
    freeMesh(theProblem->grad_px);
    freeMesh(theProblem->grad_py);
    freeMesh(theProblem->grad_phix);
    freeMesh(theProblem->grad_phiy);
    freeMesh(theProblem->divx);
    freeMesh(theProblem->divy);
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

void advective(Problem *theProblem){
    Mesh *U = theProblem->u;
    Mesh *V = theProblem->v;
    Mesh *Hx = theProblem->Hx;
    Mesh *Hy = theProblem->Hy;

    double **u = U->grid;
    double **v = V->grid;
    double **hx = Hx->grid;
    double **hy = Hy->grid;

    double h = theProblem->h;
    double udu;
    double vdu;
    for(int j = 1; j < Hx->Ny-1; j++){
        for(int i = 1; i < Hx->Nx-1; i++){
            udu = 0.5 * (0.5 * (u[j][i+1] + u[j][i]) * (u[j][i+1] - u[j][i]) / h + 0.5 * (u[j][i] + u[j][i-1]) * (u[j][i] - u[j][i-1]) / h);
            vdu = 0.5 * (0.5 * (v[j][i+1] + v[j][i]) * (u[j+1][i] - u[j][i]) / h + 0.5 * (v[j-1][i+1] + v[j-1][i]) * (u[j][i] - u[j-1][i]) / h);
            hx[j][i] = udu + vdu;
        }
    }

    // To change, this is only a copy-paste from Hx -> Should be correct now
    double udv;
    double vdv;
    for(int j = 1; j < Hy->Ny-1; j++){
        for(int i = 1; i < Hy->Nx-1; i++){
            udv = 0.5 * (0.5 * (u[j][i-1] + u[j+1][i-1]) * (v[j][i] - v[j][i-1]) / h + 0.5 * (u[j][i] + u[j+1][i]) * (v[j][i+1] - v[j][i]) / h);
            vdv = 0.5 * (0.5 * (v[j][i] + v[j-1][i]) * (v[j][i] - v[j-1][i]) / h + 0.5 * (v[j][i] + v[j+1][i]) * (v[j+1][i] - v[j][i]) / h);
            hy[j][i] = udv + vdv;
        }
    }
}

void gradP(Problem *theProblem){

    Mesh *P = theProblem->p;
    Mesh *Grad_px = theProblem->grad_px;
    Mesh *Grad_py = theProblem->grad_py;

    double **p = P->grid;
    double **grad_px = Grad_px->grid;
    double **grad_py = Grad_py->grid;

    double h = theProblem->h;

    // Computing grad_px
    for(int j = 0; j < Grad_px->Ny; j++){
        for(int i = 0; i < Grad_px->Nx; i++){
            grad_px[j][i] = (p[j][i+1] - p[j][i]) / h;
        }
    }

    // Computing grad_py
    for(int j = 0; j < Grad_py->Ny; j++){
        for(int i = 0; i < Grad_py->Nx; i++){
            grad_py[j][i] = (p[j+1][i] - p[j][i]) / h;
        }
    }
}

void gradPhi(Problem *theProblem){

    Mesh *Phi = theProblem->phi;
    Mesh *Grad_phix = theProblem->grad_phix;
    Mesh *Grad_phiy = theProblem->grad_phiy;

    double **phi = Phi->grid;
    double **grad_phix = Grad_phix->grid;
    double **grad_phiy = Grad_phiy->grid;

    double h = theProblem->h;

    // Computing grad_phix
    for(int j = 0; j < Grad_phix->Ny; j++){
        for(int i = 0; i < Grad_phix->Nx; i++){
            grad_phix[j][i] = (phi[j][i+1] - phi[j][i]) / h;
        }
    }

    // Computing grad_phiy
    for(int j = 0; j < Grad_phiy->Ny; j++){
        for(int i = 0; i < Grad_phiy->Nx; i++){
            grad_phiy[j][i] = (phi[j+1][i] - phi[j][i]) / h;
        }
    }
}

void diffusive(Problem *theProblem){
    Mesh *U = theProblem->u;
    Mesh *V = theProblem->v;
    Mesh *D1 = theProblem->d1;
    Mesh *D2 = theProblem->d2;
    Mesh *D3 = theProblem->d3;
    Mesh *Nu = theProblem->nu;
    Mesh *Divx = theProblem->divx;
    Mesh *Divy = theProblem->divy;

    double **u = U->grid;
    double **v = V->grid;
    double **d1 = D1->grid;
    double **d2 = D2->grid;
    double **d3 = D3->grid;
    double **nu = Nu->grid;
    double **divx = Divx->grid;
    double **divy = Divy->grid;

    double h = theProblem->h;

    double dnudx;
    double dnudy;
    double dd1dx;
    double dd2dx;
    double dd2dy;
    double dd3dy;

    // Computes nu_SGS in the center of the domain.
    for(int j = 1; j < Nu->Ny-1; j++){
        for(int i = 1; i < Nu->Nx-1; i++){
            d1[j][i] = 0.5 * ((u[j][i+1] - u[j][i-1]) / (2.0 * h) + (u[j+1][i+1] - u[j+1][i-1]) / (2.0 * h));
            d2[j][i] = 0.5 * ((u[j+1][i] - u[j][i]) / h + (v[j][i+1] - v[j][i]) / h);
            d3[j][i] = 0.5 * ((v[j+1][i+1] - v[j-1][i+1]) / (2.0 * h) + (v[j+1][i] - v[j-1][i]) / (2.0 * h));
            nu[j][i] = theProblem->C * pow(h, 2.0) * sqrt(2.0 * (pow(d1[j][i], 2.0) + 2.0 * pow(d2[j][i], 2.0) + pow(d3[j][i], 2.0)));
        }
    }

    // Computes nu_SGS on the left boudary (exluding top and bottom corners of the domain)
    // using a second order decentered scheme for d1. d2 and d3 are computed the same way as in the center of the domain
    for(int j = 1; j < Nu->Ny-1; j++){
        int i = 0;
        d1[j][i] = 0.5 * (1.0 / (2.0 * h) * (-3.0 * u[j][i] + 4.0 * u[j][i+1] - u[j][i+2]) + 1.0 / (2.0 * h) * (-3.0 * u[j+1][i] + 4.0 * u[j+1][i+1] - u[j+1][i+2]));
        d2[j][i] = 0.5 * ((u[j+1][i] - u[j][i]) / h + (v[j][i+1] - v[j][i]) / h);
        nu[j][i] = theProblem->C * pow(h, 2.0) * sqrt(2.0 * (pow(d1[j][i], 2.0) + 2.0 * pow(d2[j][i], 2.0) + pow(d3[j][i], 2.0)));
    }

    // Computes nu_SGS on the right boudary (exluding top and bottom corners of the domain)
    // using a second order decentered scheme for d1. d2 and d3 are computed the same way as in the center of the domain
    for(int j = 1; j < Nu->Ny-1; j++){
        int i = Nu->Nx-1;
        d1[j][i] = 0.5 * (1.0 / (2.0 * h) * (3.0 * u[j][i] - 4.0 * u[j][i-1] + u[j][i-2]) + 1.0 / (2.0 * h) * (3.0 * u[j+1][i] - 4.0 * u[j+1][i-1] + u[j+1][i-2]));
        d2[j][i] = 0.5 * ((u[j+1][i] - u[j][i]) / h + (v[j][i+1] - v[j][i]) / h);
        nu[j][i] = theProblem->C * pow(h, 2.0) * sqrt(2.0 * (pow(d1[j][i], 2.0) + 2.0 * pow(d2[j][i], 2.0) + pow(d3[j][i], 2.0)));
    }

    // Computes nu_SGS on the bottom boudary (exluding left and right corners of the domain)
    // using a second order decentered scheme for d1. d1 and d2 are computed the same way as in the center of the domain
    for(int i = 1; i < Nu->Nx-1; i++){
        int j = 0;
        d1[j][i] = 0.5 * (1.0 / (2.0 * h) * (-3.0 * v[j][i] + 4.0 * v[j+1][i] - v[j+2][i]) + 1.0 / (2.0 * h) * (-3.0 * v[j][i+1] + 4.0 * v[j+1][i+1] - v[j+2][i+1]));
        d2[j][i] = 0.5 * ((u[j+1][i] - u[j][i]) / h + (v[j][i+1] - v[j][i]) / h);
        nu[j][i] = theProblem->C * pow(h, 2.0) * sqrt(2.0 * (pow(d1[j][i], 2.0) + 2.0 * pow(d2[j][i], 2.0) + pow(d3[j][i], 2.0)));
    }

    // Computes nu_SGS on the top boudary (exluding left and right corners of the domain)
    // using a second order decentered scheme for d1. d1 and d2 are computed the same way as in the center of the domain
    for(int i = 1; i < Nu->Nx-1; i++){
        int j = Nu->Ny-1;
        d1[j][i] = 0.5 * (1.0 / (2.0 * h) * (3.0 * v[j][i] - 4.0 * v[j-1][i] + v[j-2][i]) + 1.0 / (2.0 * h) * (3.0 * v[j][i+1] - 4.0 * v[j-1][i+1] + v[j-2][i+1]));
        d2[j][i] = 0.5 * ((u[j+1][i] - u[j][i]) / h + (v[j][i+1] - v[j][i]) / h);
        nu[j][i] = theProblem->C * pow(h, 2.0) * sqrt(2.0 * (pow(d1[j][i], 2.0) + 2.0 * pow(d2[j][i], 2.0) + pow(d3[j][i], 2.0)));
    }


    for(int j = 1; j < Nu->Ny-1; j++){
        for(int i = 1; i < Nu->Nx-1; i++){
            dnudx = (nu[j][i+1] - nu[j][i-1]) / (2.0 * h);
            dnudy = (nu[j+1][i] - nu[j-1][i]) / (2.0 * h);

            dd1dx = (d1[j][i+1] - d1[j][i-1]) / (2.0 * h);

            dd2dx = (d2[j][i+1] - d2[j][i-1]) / (2.0 * h);
            dd2dy = (d2[j+1][i] - d2[j-1][i]) / (2.0 * h);

            dd3dy = (d3[j+1][i] - d3[j-1][i]) / (2.0 * h);

            divx[j][i] = nu[j][i] * (dd1dx + dd2dx) + (d1[j][i] + d2[j][i]) * dnudx;
            divy[j][i] = nu[j][i] * (dd2dy + dd3dy) + (d2[j][i] + d3[j][i]) * dnudy;
        }
    }
}

void vorticity(Problem *theProblem){
    Mesh *U = theProblem->u;
    Mesh *V = theProblem->v;
    Mesh *W = theProblem->w;

    double **u = U->grid;
    double **v = V->grid;
    double **w = W->grid;

    double dvdx;
    double dudy;

    for(int j = 0; j < theProblem->nu->Ny; j++){
        for(int i = 0; i < theProblem->nu->Nx; i++){
            dvdx = (v[j][i+1] - v[j][i]) / theProblem->h;
            dudy = (u[j+1][i] - u[j][i]) / theProblem->h;
            w[j][i] = dvdx - dudy;
        }
    }
}

void outFlow(Problem *theProblem){
    Mesh *U = theProblem->u_star;
    Mesh *U_star = theProblem->u_star;

    double **u = U->grid;
    double **u_star = U_star->grid;
    double *u_p = theProblem->u_p;

    int Ny = U->Ny;
    int Nx = U->Nx;
    // Applies the right border boudary condition (equation 6))
    for(int j = 0; j < Ny; j++){
        u_star[j][Nx-1] = u[j][Nx-1] - theProblem->dt / theProblem->h * u_p[j] * (u[j][Nx-1] - u[j][Nx-2]);
    }

    // Computes the massflow at the outlet and corrects and spreads the excess over all uelocities
    double massFlowIn = theProblem->massFlow;
    double massFlowOut = 0.0;
    for(int j = 0; j < Ny; j++){
        massFlowOut += u_star[j][Nx-1] * theProblem->h;
    }

    for(int j = 0; j < Ny; j++){
        u_star[j][Nx-1] -= (massFlowOut - massFlowIn) / theProblem->H;
    }
}
