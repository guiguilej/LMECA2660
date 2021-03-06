#ifndef _USEFULL_H_
#define _USEFULL_H_

#include "problem.h"

double** matrix(int Nx, int Ny);
void printVec(int N, double *vec);
void printMat(int N_row, int N_col, double **mat);
void saveProblem(Problem *theProblem, char *filename);
void saveMat(int Nx, int Ny, double **mat, char *name, int it);
void saveMatTest(int Nx, int Ny, double **mat);
void vecToMat(Mesh *mesh, double *vec);
#endif
