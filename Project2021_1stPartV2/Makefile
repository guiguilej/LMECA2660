
# INC_DIR := -I/usr/local/opt/petsc/include
# LIB_DIR := -L/usr/local/opt/petsc/lib

INC_DIR := -I/home/william/Project2660/lib_petsc/include
LIB_DIR := -L/home/william/Project2660/lib_petsc/lib -Wl,-rpath=/home/william/Project2660/lib_petsc/lib

#LIB_DIR := -L/home/alain/Bureau/Projet_2660/lib_petsc/lib -Wl,-rpath=/home/alain/Bureau/Projet_2660/lib_petsc/lib
LIB := -lpetsc

CXX_FLAGS := -O0 -Wall -Werror #-g

#Compilation
all :
	gcc -o project -ggdb3 project.c poisson.c usefull.c problem.c -lm $(CXX_FLAGS) $(LIB_DIR) $(LIB) $(INC_DIR)
	# gcc -o project project.c poisson.c usefull.c problem.c -lm $(CXX_FLAGS) $(LIB_DIR) $(LIB) $(INC_DIR)


# Valgrind line to debug: valgrind --leak-check=yes --show-leak-kinds=all ./project

#Delete of executable files
clean :
	rm projet

#Delete of results
clean_txt :
	rm -vf results/*.txt results/P-* results/U-* results/V-* results/Reh-* results/Vtx-* results/Rehw-* results/Div-*
