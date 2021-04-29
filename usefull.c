#include <stdio.h>
#include <stdlib.h>

#include "problem.h"
#include "usefull.h"

double** matrix(int Nx, int Ny)
{
    double **Mat = calloc(Ny, sizeof(double*));
    for(int i = 0; i < Ny; i++)
    {
        Mat[i] = calloc(Nx, sizeof(double));
    }
    return Mat;
}

void printVec(int N, double *vec){
    printf("%s", "{");
    for (int i = 0; i < N; i++) {
        printf("%f ", vec[i]);
    }
    printf("%s\n", "}");
}


void printMat(int Nx, int Ny, double **mat){
    printf("%s\n", "{");
    for(int i = 0; i < Nx; i ++){
        for(int j = 0; j < Ny; j ++){
            printf("%3.2f  ", mat[i][j]);
        }
        printf("%s", "\n");
    }
    printf("%s\n", "}");
}

void saveProblem(Problem *theProblem, char *filename){
    printf("Writing in file: %s\n", filename);
    FILE *file = fopen(filename, "w");
    fprintf(file, "%.10f\n", theProblem->h);
    fprintf(file, "%.10f\n", theProblem->H);
    fprintf(file, "%.10f\n", theProblem->L);
    fprintf(file, "%d\n", theProblem->Nx);
    fprintf(file, "%d\n", theProblem->Ny);
    fprintf(file, "%.10f\n", theProblem->h_hill);
    fprintf(file, "%.10f\n", theProblem->d_hill);
    fprintf(file, "%.10f\n", theProblem->sigma_hill);
    fprintf(file, "%.10f\n", theProblem->u_hill);
    fprintf(file, "%.10f\n", theProblem->u_tau);
    fprintf(file, "%.10f\n", theProblem->C);
    fclose(file);
}

void saveMat(int Nx, int Ny, double **mat, char *name){
    char *basename = "/mnt/d/cours/Q8/LMECA2660/Project/Data/Results_%s.txt";
    char filename[50];
    sprintf(filename, basename, name);
    printf("Writing in file: %s\n", filename);
    FILE *file = fopen(filename, "w");
    for(int i = 0; i < Ny; i++){
        for(int j = 0; j <Nx; j++){
            fprintf(file, "%.8f ", mat[i][j]);
        }
        fprintf(file, "\n");
    }
    fclose(file);
}
