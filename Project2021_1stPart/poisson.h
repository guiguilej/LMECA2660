#ifndef _POISSON_H_
#define _POISSON_H_

/*To include in the file in which you will call initialize_poisson_solver and poisson_solver*/

#include <petsc.h>
#include <petscsys.h>

//Structure storing petsc vectors

typedef struct {

	Vec b;
	Vec x;
	Mat A;
	KSP sles;

} Poisson_data;

PetscErrorCode initialize_poisson_solver(Poisson_data* data, Problem *theProblem);
void poisson_solver(Poisson_data *data, Problem *theProblem);
void free_poisson_solver(Poisson_data* data);
void computeRHS(double *rhs, Problem* theProblem);
void computeLaplacianMatrix(Mat A, Problem *theProblem);

#endif
