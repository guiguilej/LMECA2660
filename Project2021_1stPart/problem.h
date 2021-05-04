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
    int it;
    double h;
    double *u_p;
    double *rhs;
    double *x_poisson;

    double dt;
    double dtau;
    double Re;
    double t_max;

    double massFlow;

    Mesh *Hx;
    Mesh *Hy;

    Mesh *Hx_old;
    Mesh *Hy_old;

    Mesh *divx;
    Mesh *divy;

    Mesh *u;
    Mesh *v;
    Mesh *u_star;
    Mesh *v_star;
    Mesh *p;
    Mesh *phi;
    Mesh *d1;
    Mesh *d2;
    Mesh *d3;
    Mesh *nu;
    Mesh *w;
    Mesh *grad_px;
    Mesh *grad_py;
    Mesh *grad_phix;
    Mesh *grad_phiy;
} Problem;


Problem *initProblem(double h_hill, double u_hill, double y0);
Mesh *initMesh(int Nx, int Ny, double h);
void initialCondition(Problem *theProblem);
double *initU_p(Problem *theProblem);
void boundaryConditions(Problem *theProblem);
void outFlow(Problem *theProblem);
void advective(Problem *theProblem);
void gradP(Problem *theProblem);
void gradPhi(Problem *theProblem);
void diffusive(Problem *theProblem);
void vorticity(Problem *theProblem);
void freeMesh(Mesh *mesh);
void freeProblem(Problem *theProblem);

#endif
