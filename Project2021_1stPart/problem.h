#ifndef _PROBLEM_H_
#define _PROBLEM_H_


typedef struct{
    int Nx;
    int Ny;
    double h;
    double **grid;
} Mesh;

typedef struct{
    double h_hill;
    double H;
    double L;
    double d_hill;
    double sigma_hill;
    double u_hill;
    double y0;
    double kappa;
    double u_tau;
    double C;
    int Nx;
    int Ny;
    double h;
    double *u_p;

    double dt;
    double dtau;
    double Re;

    Mesh *Hx;
    Mesh *Hy;

    Mesh *divx;
    Mesh *divy;

    Mesh *u;
    Mesh *v;
    Mesh *p;
    Mesh *d1;
    Mesh *d2;
    Mesh *d3;
    Mesh *nu;
    Mesh *w;
    Mesh *grad_px;
    Mesh *grad_py;
} Problem;


Problem *initProblem(double h_hill, double u_hill, double y0);
Mesh *initMesh(int Nx, int Ny, double h);
void initialCondition(Problem *theProblem);
double *initU_p(Problem *theProblem);
void boundaryConditions(Problem *theProblem);
void advective(Problem *theProblem);
void gradP(Problem *theProblem);
void diffusive(Problem *theProblem);
void freeMesh(Mesh *mesh);
void freeProblem(Problem *theProblem);

#endif
