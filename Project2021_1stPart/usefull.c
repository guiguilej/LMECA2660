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
    for(int j = 0; j < Ny; j ++){
        for(int i = 0; i < Nx; i ++){
            printf("%3.2f  ", mat[j][i]);
        }
        printf("%s", "\n");
    }
    printf("%s\n", "}");
}

void saveProblem(Problem *theProblem, char *filename){

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

void saveMat(int Nx, int Ny, double **mat, char *name, int it){
    char *basename = "/mnt/d/cours/Q8/LMECA2660/Project/Data/Results2/%s_%d.txt";
    char filename[100];
    sprintf(filename, basename, name, it);
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

void saveMatTest(int Nx, int Ny, double **mat){
    char *filename = "/mnt/d/cours/Q8/LMECA2660/Project/Data/mat.txt";
    printf("Writing in file: %s\n", filename);
    FILE *file = fopen(filename, "w");
    for(int i = 0; i < Ny; i++){
        for(int j = 0; j <Nx; j++){
            double val = mat[i][j];
            if(val == 0){
                fprintf(file, "%s", "  ");
            }
            else{
                fprintf(file, "%.1f ", val);
            }
        }
        fprintf(file, "  %d \n", i+1);
    }
    fclose(file);
}

void vecToMat(Mesh *mesh, double *vec){
    double **grid = mesh->grid;
    int Ny = mesh->Ny;
    int Nx = mesh->Nx;
    int index = 0;

    for(int j = 0; j < Ny; j++){
        // printf("%d\n", j);
        for(int i = 0; i < Nx; i++){
            index = j * Nx + i;
            grid[j][i] = vec[index];
        }
    }
}
