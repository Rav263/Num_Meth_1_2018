#include <stdio.h>
#include <stdlib.h>
#include <math.h>


// Выделение памяти для матрицы
double **init_matrix(int n) {
    double **matrix = calloc(n, sizeof(*matrix));

    if (matrix == NULL) {
        fprintf(stderr, "MEMMORY ERROR\n");
        exit(1);
    }
    
    for (int i = 0; i < n; i++) {
        matrix[i] = calloc(n, sizeof(**matrix));
        if (matrix[i] == NULL) {
            fprintf(stderr, "MEMMORY ERROR\n");
            exit(1);
        }
    }

    return matrix;
}

// печать матрицы
void print_matrix(double **matrix, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            printf("%lf ", matrix[i][j]);
        }
        printf("\n");
    }
}


// Чтение матрицы
void scan_matrix(double **matrix, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            scanf("%lf", &matrix[i][j]);
        }
    }
}

// Умножение строки на константу
void string_diff(double *string_1, double *string_2, int n, double coof) {
    for (int i = 0; i < n; i++) {
        string_1[i] -= coof * string_2[i];
    }
}

// Прямой ход метода Гауса
void gauss_direct_way(double **matrix, double *sol, int n) {
    for (int i = 0; i < n - 1; i++) {
        if (matrix[i][i] == 0) {
            int index = i + 1;

            for (; index < n; index++) {
                if (matrix[index][i] != 0) break;
            }

            double *tmp = matrix[index];
            matrix[index] = matrix[i];
            matrix[i] = tmp;
        }
    
        print_matrix(matrix, n);

        for (int j = i + 1; j < n; j++) {
            double coof = matrix[j][i] / matrix[i][i];        
            printf("Coof: %lf\n J: %d\n", coof, j); 
            string_diff(matrix[j], matrix[i], n, coof);
        }
    }
    print_matrix(matrix, n);
}

// Обратный ход метзода Гауса
void gauss_return_way(double **matrix, double *sol, int n) {
    
}

// Подсчёт определителя
double det(double **matrix, int n) {
    gauss_direct_way(matrix, NULL ,n);

    double result = 1;

    for (int i = 0; i < n; i++){
        result *= matrix[i][i];
    }

    return result;
}

double **copy_matrix(double **matrix, int n) {
    double **new_matrix = init_matrix(n);

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            new_matrix[i][j] = matrix[i][j];
        }
    }

    return new_matrix;
}


int main(int argc, char *argv[]) {
    if (argc < 3) {
        fprintf(stderr, "Too few arguments!\nUsage ./prog type mtrix_size\n");
        exit(1);
    }

    int type = atoi(argv[1]); // Чтение со стандартного потока или из файла
    int n    = atoi(argv[2]); // Размер матрицы

    double **matrix = init_matrix(n);
    double **rev_matrix = init_matrix(n);

    if (type == 1) scan_matrix(matrix, n);

    double determ = det(matrix, n);

    printf("%lf\n", determ);

}
