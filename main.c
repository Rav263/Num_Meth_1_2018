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
void print_matrix(double **matrix, int n, int accur) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            printf("%*.*lf ",accur + 3, accur, matrix[i][j]);
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

void scan_file_matrix(char *file_name, double **matrix, double *right, int n) {
    FILE *file = fopen(file_name, "r");

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            fscanf(file, "%lf", &matrix[i][j]);
        }
        
        fscanf(file, "%lf", &right[i]);
    }
    fclose(file);
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
            if (index == n) { 
                fprintf(stderr, "Wrong Input\n");
                exit(1);
            }

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

        for (int j = i + 1; j < n; j++) {
            double coof = matrix[j][i] / matrix[i][i];        
            if (rev_matrix != NULL) string_diff(rev_matrix[j], rev_matrix[i], n, coof);
            string_diff(matrix[j], matrix[i], n, coof);
            if (right != NULL) right[j] -= right[i] * coof;
        }
    }
}



// Обратный ход метзода Гаусса
double *gauss_return_way(double **matrix, double **rev_matrix, double *right, int n) {
    double *answer = calloc(n, sizeof(*answer));
    
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
    }

    return answer;
}

double d_abs(double a) {
    return a < 0 ? -a : a;
}


// Поиск максимума по модулю в строке 
int max_abs(double *string, int start_index, int n) {
    double max = string[start_index];
    int max_index = start_index;

    for (int i = start_index; i < n; i++) {
        max_index = d_abs(string[i]) > d_abs(max) ? i : max_index;
        max = d_abs(string[i]) > d_abs(max) ? string[i] : max;
    }

    return max_index;
}

// Смена столбцов местами
void swap_col(double **matrix, int index_1, int index_2, int n) {
    for (int i = 0; i < n; i++) {
        double tmp = matrix[i][index_1];
        matrix[i][index_1] = matrix[i][index_2];
        matrix[i][index_2] = tmp;
    }
}


// Метод Гаусса с выбором главного элемента прямой ход
int *mainel_gauss_direct_way(double **matrix, double *right, int n) {
    int *rev_num = calloc(n, sizeof(*rev_num));
    
    for (int i = 0; i < n; i++)
        rev_num[i] = i;
    
    for (int i = 0; i < n; i++) {
        int max_index = max_abs(matrix[i], i, n);

        if (matrix[i][max_index] == 0) {
            fprintf(stderr, "Wrong Input\n");
            exit(1);
        }

        swap_col(matrix, i, max_index, n);

        int tmp = rev_num[i];
        rev_num[i] = rev_num[max_index];
        rev_num[max_index] = tmp;

        for (int j = i + 1; j < n; j++) {
            double coof = matrix[j][i] / matrix[i][i];

            string_diff(matrix[j], matrix[i], n, coof);
            right[j] -= right[i] * coof;
        }
    }

    return rev_num;
}

// Обратный ход метода Гаусса с выбором главного элемента
double *mainel_gauss_return_way(double **matrix, double *right, int *rev_num, int n) {
    double *answer = calloc(n, sizeof(*answer));

    for (int i = n - 1; i >= 0; i--) {
        answer[i] = right[i];

        for (int j = i + 1; j < n; j++) 
            answer[i] -= matrix[i][j] * answer[j];

        answer[i] /= matrix[i][i];
    }

    double *tmp_answer = calloc(n, sizeof(*tmp_answer));

    for (int i = 0; i < n; i++) 
        tmp_answer[rev_num[i]] = answer[i];

    free(answer);
    return tmp_answer;
}

// Метод Гаусса с выбором главного элемента
double *mainel_gauss(double **matrix, double *right, int n) {
    int *rev_num = mainel_gauss_direct_way(matrix, right, n);

    double *answer = mainel_gauss_return_way(matrix, right, rev_num, n);
    free(rev_num);

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
    
    gauss_return_way(matrix, rev_matrix, NULL, n);
}

// Инициализация единичной матрицы
void init_unit_matrix(double **matrix, int n) {
    for (int i = 0; i < n; i++) {
        matrix[i][i] = 1;
    }
}

double *copy_string(double *string, int n) {
    double *new_string = calloc(n, sizeof(*new_string));

    for (int i = 0; i < n; i++) 
        new_string[i] = string[i];

    return new_string;
}

void printf_string(double *string, int n, int accur) {
    for (int i = 0; i < n; i++) 
        printf("%*.*lf ", accur + 3, accur, string[i]);
    printf("\n");
}

void free_matrix(double **matrix, int n) {
    for (int i = 0; i < n; i++) 
        free(matrix[i]);
    free(matrix);
}


void count_rev_matrix(double **matrix, int n, int accur) {
    double **rev_matrix = init_matrix(n);
    init_unit_matrix(rev_matrix, n);
    
    double **tmp_matrix = copy_matrix(matrix, n);
    matrix_reverce(tmp_matrix, rev_matrix, n);

    printf("Reverce matrix:\n");
    print_matrix(rev_matrix, n, accur);
    printf("---------------\n");

    free_matrix(rev_matrix, n);
    free_matrix(tmp_matrix, n);
}


void count_det(double **matrix, int n, int accur) {
    double **tmp_matrix = copy_matrix(matrix, n);

    double deter = det(tmp_matrix, n);

    printf("Determinant: %*.*lf\n\n", accur + 3, accur, deter);

    free_matrix(tmp_matrix, n);
}

void count_gauss(double **matrix, double *right, int n, int accur) {
    double **tmp_matrix = copy_matrix(matrix, n);
    double *tmp_right = copy_string(right, n);

    double *answer = gauss(tmp_matrix, right, n);

    printf("Gauss answer: ");

    printf_string(answer, n, accur);

    free(tmp_right);
    free(answer);
    free_matrix(tmp_matrix, n);
}

void count_mainel_gauss(double **matrix, double *right, int n, int accur) {
    double **tmp_matrix = copy_matrix(matrix, n);
    double *tmp_right = copy_string(right, n);

    double *answer = mainel_gauss(tmp_matrix, right, n);

    printf("Main element Gauss: ");
    printf_string(answer, n, accur);

    free(tmp_right);
    free(answer);
    free_matrix(tmp_matrix, n);
}


int generate_matrix(double **matrix, double *right, double x, int n) {
    int M = 3;
    double qm = 1.001 - 2 * M * 0.001;


    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= n; j++) {
            if (i == j) matrix[i - 1][j - 1] = pow(qm - 1, i + j);
            else matrix[i - 1][j - 1] = pow(qm, i + j) + 0.1 * (i - j);
        }
    }

    for (int i = 1; i <= n; i++) {
        right[i - 1] = x * exp(x / i) * cos(x / i);
    }

    return n;
}

double *next_iteration(double **matrix, double *x, double *right, double w, int n) {
    double *new_x = calloc(n, sizeof(*new_x));

    for (int i = 0; i < n; i++) {
        double fir = 0;
        
        for (int j = 0; j < i; j++) {
            fir += matrix[j][i] * new_x[j];
        }

        double sec = 0;

        for (int j = i; j < n; j++) {
            sec += matrix[j][i] * x[j];
        }

        new_x[i] = x[i] + w / matrix[i][i] * (right[i] - fir - sec);
    }

    free(x);
    return new_x;
}


double norma(double **matrix, double *right, double *x, int n) {
    double *new_f = copy_string(right, n);

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            new_f[i] -= matrix[i][j] * x[j];
        }
    }

    double norm = 0;

    for (int i = 0; i < n; i++) {
        norm += new_f[i] * new_f[i];
    }

    norm = sqrt(norm);

    free(new_f);

    return norm;
}


double *iteration_method(double **matrix, double *right, double w, double eps, int max_iter, int n) {
    double *x = calloc(n, sizeof(*x));

    for (int i = 0; i < max_iter; i++) {
        x = next_iteration(matrix, x, right, w, n);

        double norm = norma(matrix, right, x, n);

        if (norm < eps) {
            fprintf(stderr, "Number of iterations: %d\n", i);
            break;
        }
    }

    return x;
}

void count_iter(double **matrix, double *right, int n, int accur) {
    int max_iter = 10000;
    double eps = 0.1;
    printf("Please enter max iterations: ");
    scanf("%d", &max_iter);
    printf("Please enter eps: ");
    scanf("%lf", &eps);

    for (double w = 0.1; w < 2; w += 0.1) {
        double *x = iteration_method(matrix, right, w, eps, max_iter, n);

        printf("W: %*.*lf, answ: ", accur + 3, accur, w);
        printf_string(x, n, accur);
        free(x);
    }
}

double l_norm(double **matrix, int n) {
    double *sums = calloc(n, sizeof (*sums));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            sums[i] += d_abs(matrix[i][j]);
        }
    }

    double max = sums[max_abs(sums, 0, n)];
    
    free(sums);

    return max;
}


void count_l_norm(double **matrix, int n, int accur) {
    double **tmp_matrix = copy_matrix(matrix, n);
    double **rev_matrix = init_matrix(n);

    init_unit_matrix(rev_matrix, n);
    
    matrix_reverce(tmp_matrix, rev_matrix, n);

    free(tmp_matrix);

    double norm_1 = l_norm(matrix, n);
    double norm_2 = l_norm(rev_matrix, n);

    double norm = norm_1 * norm_2;

    printf("Condition number of matrix: %*.*lf\n", accur + 3, accur, norm);

    free(rev_matrix);
}



int main(int argc, char *argv[]) {
    if (argc < 3) {
        fprintf(stderr, "Too few arguments!\nUsage ./prog type mtrix_size\n");
        fprintf(stderr, "Ussage: \n1 -- reverce matrix\n2 -- determinant\n3 -- normal gauss\n4 -- main element gauss\n");
        exit(1);
    }
    int type = atoi(argv[1]);
    int n = 0;
    if (type != 3) n= atoi(argv[type + 1]);
    else  n = 10;

    double **matrix = init_matrix(n);
    double *right   = calloc(n, sizeof(*right));

    int accur;
    printf("Please enter accurate: ");
    scanf("%d", &accur);

    if (type == 1) scan_matrix(matrix, right, n);
    else if (type == 2) scan_file_matrix(argv[2], matrix, right, n);
    else if (type == 3) {
        double x = atof(argv[2]);
        

        n = generate_matrix(matrix, right, x, n);

        print_matrix(matrix, n, accur);
        printf("\n\n");
        printf_string(right, n, accur);

        type -= 2;
    }


    for (int i = type + 2; i < argc; i++) {
        int now = atoi(argv[i]);

        if (now == 1) count_rev_matrix(matrix, n, accur);
        else if (now == 2) count_det(matrix, n, accur);
        else if (now == 3) count_gauss(matrix, right, n, accur);
        else if (now == 4) count_mainel_gauss(matrix, right, n, accur);
        else if (now == 5) count_iter(matrix, right, n, accur);
        else if (now == 6) count_l_norm(matrix, n, accur);
    }


    free(right);
    free_matrix(matrix, n);
}
