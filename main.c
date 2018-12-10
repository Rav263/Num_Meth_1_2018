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
void scan_matrix(double **matrix, double *right, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            scanf("%lf", &matrix[i][j]);
        }

        scanf("%lf", &right[i]);
    }
}

// Умножение строки на константу
void string_diff(double *string_1, double *string_2, int n, double coof) {
    for (int i = 0; i < n; i++) {
        string_1[i] -= coof * string_2[i];
    }
}


// умножение строки на коэффицент
void string_mult(double *string, int n, double coof) {
    for (int i = 0; i < n; i++) {
        string[i] *= coof;
    }
}


// Прямой ход метода Гауса
void gauss_direct_way(double **matrix, double **rev_matrix, double *right, int n) {
    for (int i = 0; i < n - 1; i++) {
        if (matrix[i][i] == 0) {
            int index = i + 1;

            for (; index < n; index++) {
                if (matrix[index][i] != 0) break;
            }
            if (index == n) continue;

            double *tmp = matrix[index];
            matrix[index] = matrix[i];
            matrix[i] = tmp;

            // Для случая поиска обратной матрицы
            if (rev_matrix != NULL) {
                tmp = rev_matrix[index];
                rev_matrix[index] = rev_matrix[i];
                rev_matrix[i] = tmp;
            }

            // Для случая решения СЛАУ
            if (right != NULL) {
                double num_tmp = right[index];
                right[index] = right[i];
                right[i] = num_tmp;
            }
        }

        //if (rev_matrix != NULL) string_mult(rev_matrix[i], n, 1 / matrix[i][i]);

        for (int j = i + 1; j < n; j++) {
            double coof = matrix[j][i] / matrix[i][i];        
            //printf("Coof: %lf\n J: %d\n", coof, j); 
            if (rev_matrix != NULL) string_diff(rev_matrix[j], rev_matrix[i], n, coof);
            string_diff(matrix[j], matrix[i], n, coof);
            if (right != NULL) right[j] -= right[i] * coof;
        }

        if (rev_matrix != NULL) {
            print_matrix(matrix, n);
            printf("-------------\n");
            print_matrix(rev_matrix, n);
            printf("\n\n\n");
        }
    }
}



// Обратный ход метзода Гаусса
double *gauss_return_way(double **matrix, double **rev_matrix, double *right, int n) {
    double *answer = calloc(n, sizeof(*answer));
    printf("Return Way\n");

    for (int i = n - 1; i>= 0; i--) {
        if (right != NULL) answer[i] = right[i];
        if (rev_matrix != NULL) string_mult(rev_matrix[i], n, 1 / matrix[i][i]);
        

        for (int j = i + 1; j < n; j++) {
            if (rev_matrix != NULL) {
                string_diff(rev_matrix[i], rev_matrix[j], n, matrix[i][j] / matrix[i][i]);
            }
            answer[i] -= matrix[i][j] * answer[j];
        }

        answer[i] /= matrix[i][i];
        if(rev_matrix != NULL) print_matrix(rev_matrix, n);
        printf("-------\n");
        print_matrix(matrix, n);
        printf("\n\n");
    }

    return answer;
}

// Подсчёт определителя
double det(double **matrix, int n) {
    gauss_direct_way(matrix, NULL, NULL, n);

    double result = 1;

    for (int i = 0; i < n; i++){
        result *= matrix[i][i];
    }

    return result;
}


// Метод Гаусса
double *gauss(double **matrix, double *right, int n) {
    gauss_direct_way(matrix, NULL, right, n);

    return gauss_return_way(matrix, NULL, right, n);
}

// Копирование матрицы
double **copy_matrix(double **matrix, int n) {
    double **new_matrix = init_matrix(n);

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            new_matrix[i][j] = matrix[i][j];
        }
    }

    return new_matrix;
}


// Поиск обратной матрицы
void matrix_reverce(double **matrix, double **rev_matrix, int n) {
    gauss_direct_way(matrix, rev_matrix, NULL, n);
    print_matrix(rev_matrix, n);
    printf("\n");
    print_matrix(matrix, n);

    gauss_return_way(matrix, rev_matrix, NULL, n);
}

// Инициализация единичной матрицы
void init_unit_matrix(double **matrix, int n) {
    for (int i = 0; i < n; i++) {
        matrix[i][i] = 1;
    }
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
    init_unit_matrix(rev_matrix, n);
    double *right = calloc(n, sizeof(*right));

    if (type == 1) scan_matrix(matrix, right, n);

    double **tmp_matrix = copy_matrix(matrix, n);
    double determ = det(tmp_matrix, n);
    
    tmp_matrix = copy_matrix(matrix, n);
    double *answer = gauss(tmp_matrix, right, n);

    tmp_matrix = copy_matrix(matrix, n);
    matrix_reverce(tmp_matrix, rev_matrix, n);

    printf("Determ: %lf\n", determ);
    print_matrix(rev_matrix, n);

    for(int i = 0; i < n; i++) {
        printf("%lf ", answer[i]);
    }
    printf("\n");
}
