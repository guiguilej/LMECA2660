
// #include <mpi.h>
#include "poisson.h"
#include "usefull.h"

/*Called by poisson_solver at each time step*/
/*More than probably, you should need to add arguments to the prototype ... */
/*Modification to do :*/
/*    -Impose zero mass flow here by changing value of U_star*/
/*    -Fill vector rhs*/
void computeRHS(double *rhs, PetscInt rowStart, PetscInt rowEnd, Problem* theProblem){

  int Nx = theProblem -> Nx;
  int Ny = theProblem -> Ny;
  Mesh *U = theProblem -> u_star;
  Mesh *V = theProblem -> v_star;
  double** u_star = U -> grid;
  double** v_star = Mesh -> grid;
  double h = theProblem -> h;
  double dt = theProblem -> dt;

  double IntegralUout;
  double IntegralUin;
  double cor;
  for(int j=0;j<Ny;j++){
       IntegralUin += u_star[j][1];
       IntegralUout += u_star[j][Nx];
  }

  cor = (IntegralUin - IntegralUout)/(double)Ny;

  for(int j=0;j<Ny;j++){
    u_star[j][Nx] += cor;
  }

  for(int j=0;j<Ny;j++){
    for(int i=0;i<Nx;i++){
        rhs[j*Nx + i] = (h/dt) * (u_star[j][i+1] - u_star[j][i] + v_star[j+1][i] - v_star[j][i]); 
       }
  }
  rhs[0]=0;
}

/*To call at each time step after computation of U_star. This function solves the poisson equation*/
/*and copies the solution of the equation into your vector Phi*/
/*More than probably, you should need to add arguments to the prototype ... */
/*Modification to do :*/
/*    - Change the call to computeRHS as you have to modify its prototype too*/
/*    - Copy solution of the equation into your vector PHI*/
void poisson_solver(Poisson_data *data, Problem* theProblem){

    /* Solve the linear system Ax = b for a 2-D poisson equation on a structured grid */
    int its;
    PetscInt rowStart, rowEnd;
    PetscScalar *rhs, *sol;

    KSP sles = data->sles;
    Vec b = data->b;
    Vec x = data->x;
    double *x_poisson = theProblem -> x_poisson;

    /* Fill the right-hand-side vector : b */
    VecGetOwnershipRange(b, &rowStart, &rowEnd);
    VecGetArray(b, &rhs);
    computeRHS(rhs, rowStart, rowEnd, Problem* theProblem); /*MODIFY THE PROTOTYPE HERE*/
    VecRestoreArray(b, &rhs);


    /*Solve the linear system of equations */
    KSPSolve(sles, b, x);
    KSPGetIterationNumber(sles, &its);
    PetscPrintf(PETSC_COMM_WORLD, "Solution to Poisson eqn in %d iterations \n", its);

    VecGetArray(x, &sol);

    int r;
    for(r=rowStart; r<rowEnd; r++){
        x_poisson[r] = sol[r];
    }

    VecRestoreArray(x, &sol);
    free(sol);
    free(rhs);
}

/*This function is called only once during the simulation, i.e. in initialize_poisson_solver.*/
/*In its current state, it inserts unity on the main diagonal.*/
/*More than probably, you should need to add arguments to the prototype ... .*/
/*Modification to do in this function : */
/*   -Insert the correct factor in matrix A*/
void computeLaplacianMatrix(Mat A, int N_row, int N_col){

    int Nx = theProblem -> Nx;
    int Ny = theProblem -> Ny;
    MatSetValue(A,0,0,1,INSERT_VALUES);
    MatSetValue(A,Nx,0,1,INSERT_VALUES);

    for (int r = 1; r < Nx*Ny; r++){

        if (r<Nx-1){
          MatSetValue(A,r,r,-3,INSERT_VALUES);
          MatSetValue(A,r,r+1,1,INSERT_VALUES);
          MatSetValue(A,r,r-1,1,INSERT_VALUES);
        }
        else if(r==Nx-1){
          MatSetValue(A,r,r,-2,INSERT_VALUES);
          MatSetValue(A,r,r-1,1,INSERT_VALUES);
        }

        else if(r%Nx == 0){
          if (r == (Nx*Ny)-Nx){
            MatSetValue(A,r,r,-2,INSERT_VALUES);
            MatSetValue(A,r,r+1,1,INSERT_VALUES);
          }
          else{
            MatSetValue(A,r,r,-3,INSERT_VALUES);
            MatSetValue(A,r,r+1,1,INSERT_VALUES);
          }
        }

        else if ( (r+1)%Nx == 0){
          if(r==(Nx*Ny)-1){
            MatSetValue(A,r,r,-2,INSERT_VALUES);
            MatSetValue(A,r,r-1,1,INSERT_VALUES);
          }
          else{
            MatSetValue(A,r,r,-3,INSERT_VALUES);
            MatSetValue(A,r,r-1,1,INSERT_VALUES);
          }
        }

        if(r<(Nx*Ny - Nx)){
          MatSetValue(A,r,r+Nx,1,INSERT_VALUES);
          MatSetValue(A,r+Nx,r,1,INSERT_VALUES);
        }
      }


        /*USING MATSETVALUE FUNCTION, INSERT THE GOOD FACTOR AT THE GOOD PLACE*/
        /*Be careful; the solution from the system solved is defined within a constant.
        One point from Phi must then be set to an abritrary constant.*/
    // printMat(N_row*N_row,N_col* N_col, matjk);
}

/*To call during the initialization of your solver, before the begin of the time loop*/
/*Maybe you should need to add an argument to specify the number of unknows*/
/*Modification to do in this function :*/
/*   -Specify the number of unknows*/
/*   -Specify the number of non-zero diagonals in the sparse matrix*/
PetscErrorCode initialize_poisson_solver(Poisson_data* data){
    PetscInt rowStart; /*rowStart = 0*/
    PetscInt rowEnd; /*rowEnd = the number of unknows*/
    PetscErrorCode ierr;

	  int nphi = Nx*Ny; /*WRITE HERE THE NUMBER OF UNKNOWS*/

    /* Create the right-hand-side vector : b */
    VecCreate(PETSC_COMM_WORLD, &(data->b));
    VecSetSizes(data->b, PETSC_DECIDE, nphi);
    VecSetType(data->b,VECSTANDARD);

    /* Create the solution vector : x */
    VecCreate(PETSC_COMM_WORLD, &(data->x));
    VecSetSizes(data->x, PETSC_DECIDE, nphi);
    VecSetType(data->x,VECSTANDARD);

    /* Create and assemble the Laplacian matrix : A  */
    MatCreate(PETSC_COMM_WORLD, &(data->A));
    MatSetSizes(data->A, PETSC_DECIDE, PETSC_DECIDE, nphi , nphi);
    MatSetType(data->A, MATAIJ);
    MatSeqAIJSetPreallocation(data->A,5, NULL); // /*SET HERE THE NUMBER OF NON-ZERO DIAGONALS*/
    MatGetOwnershipRange(data->A, &rowStart, &rowEnd);

    computeLaplacianMatrix(data->A, 5, 5);
    ierr = MatAssemblyBegin(data->A, MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr);
    ierr = MatAssemblyEnd(data->A, MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr);

    /* Create the Krylov context */
    KSPCreate(PETSC_COMM_WORLD, &(data->sles));
    KSPSetOperators(data->sles, data->A, data->A);
    KSPSetType(data->sles,KSPGMRES); //KSPGMRES seems the best, it will not be used if PC LU.
    PC prec;
    KSPGetPC(data->sles, &prec);
    PCSetType(prec,PCLU);
    //KSPSetFromOptions(data->sles); // to uncomment if we want to specify the solver to use in command line. Ex: mpirun -ksp_type gmres -pc_type gamg
    KSPSetTolerances(data->sles, 1.e-12, 1e-12, PETSC_DEFAULT, PETSC_DEFAULT);
    KSPSetReusePreconditioner(data->sles,PETSC_TRUE);
    KSPSetUseFischerGuess(data->sles,1,4);
    KSPGMRESSetPreAllocateVectors(data->sles);

    PetscPrintf(PETSC_COMM_WORLD, "Assembly of Mattrix and Vectors is done \n");

    return ierr;
}

/*To call after the simulation to free the vectors needed for Poisson equation*/
/*Modification to do : nothing */
void free_poisson_solver(Poisson_data* data){
    MatDestroy(&(data->A));
    VecDestroy(&(data->b));
    VecDestroy(&(data->x));
    KSPDestroy(&(data->sles));
}
