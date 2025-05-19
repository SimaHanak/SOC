#include "functions.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>


double h = 1e-4;
double acc = 1e-5;
double max_h = 0.1;

double* initialize_velocity(double* state_vector, double E, double L_z) {
    double** g_val = make_g(state_vector[2], 0);
    double** g_inv_val = make_g_inv(state_vector[2], 0);
    state_vector[4] = - g_inv_val[0][0]*E + g_inv_val[1][0]*L_z;
    state_vector[5] = g_inv_val[1][1]*L_z - g_inv_val[1][0]*E;
    state_vector[7] = sqrt((- 1
                            - g_val[0][0] * state_vector[4] * state_vector[4]
                            - g_val[1][1] * state_vector[5] * state_vector[5] 
                            - 2*g_val[1][0] * state_vector[4] * state_vector[5] 
                            - g_val[2][2] * state_vector[6] * state_vector[6])/g_val[3][3]);

    free_g(g_val);
    free_g_inv(g_inv_val);

    return state_vector;
}

double* eq_of_motion(double* state_vector) {
    double*** Christoffel = generate_Christoffel_symbols(state_vector[2], state_vector[3]);

    double* dydt = (double*)malloc(8* sizeof(double));
    double vel[4];
    for (int i = 0; i < 4; i++) {
        dydt[i] = state_vector[i+4];
        vel[i] = state_vector[i+4];
    }
    for (int i = 4; i < 8; i++) {
        dydt[i] = 0;
    }

    for (int coords = 0; coords < 4; coords++) {
        dydt[coords + 4] = 0;
        for (int kappa = 0; kappa < 4; kappa++) {
            for (int lambda = 0; lambda < 4; lambda++) {
                dydt[coords + 4] -= Christoffel[coords][kappa][lambda] * vel[kappa] * vel[lambda];
            }
        }
    }
    free_Christoffel(Christoffel);

    return dydt;
}

double* rk4(double* state_vector) {
    double input[8];
    double* k1 = eq_of_motion(state_vector);
    for (int i = 0; i < 8; i++) {
        input[i] = state_vector[i] + h * k1[i]/2.0;
    }
    double* k2 = eq_of_motion(input);
    for (int i = 0; i < 8; i++) {
        input[i] = state_vector[i] + h * k2[i]/2.0;
    }
    double* k3 = eq_of_motion(input);
    for (int i = 0; i < 8; i++) {
        input[i] = state_vector[i] + h * k3[i];
    }
    double* k4 = eq_of_motion(input);

    double* result = (double*)malloc(8 * sizeof(double));
    for (int i = 0; i < 8; i++) {
        result[i] = state_vector[i] + h * (k1[i] + 2*k2[i] + 2*k3[i] + k4[i])/6.0;
    }

    free(k1);
    free(k2);
    free(k3);
    free(k4);
    return result;
}

double** make_CashKarp_b() {
    double** b = malloc(6 * sizeof(double*));
    for (int i = 0; i < 6; i++) {
        b[i] = calloc(5, sizeof(double));
    }
    b[1][0] = 1.0/5.0;
    b[2][0] = 3.0/40.0;
    b[3][0] = 3.0/10.0;
    b[4][0] = - 11.0/54.0;
    b[5][0] = 1631.0/55296.0;
    b[2][1] = 9.0/40.0;
    b[3][1] = - 9.0/10.0;
    b[4][1] = 5.0/2.0;
    b[5][1] = 175.0/512.0;
    b[3][2] = 6.0/5.0;
    b[4][2] = - 70.0/27.0;
    b[5][2] = 575.0/13824.0;
    b[4][3] = 35.0/27.0;
    b[5][3] = 44275.0/110592.0;
    b[5][4] = 253.0/4096.0;
    return b;
}

double** make_CashKarp_c() {
    double** c = malloc(6 * sizeof(double*));
    for (int i = 0; i < 6; i++) {
        c[i] = calloc(2, sizeof(double));
    }
    c[0][0] = 37.0/378.0;
    c[2][0] = 250.0/621.0;
    c[3][0] = 125.0/594.0;
    c[5][0] = 512.0/1771.0;
    c[0][1] = 2825.0/27648.0;
    c[2][1] = 18575.0/48384.0;
    c[3][1] = 13525.0/55296.0;
    c[4][1] = 277.0/14336.0;
    c[5][1] = 1.0/4.0;
    return c;
}

void free_CashKarp(double** b, double** c) {
    for (int i = 0; i < 6; i++) {
        free(b[i]);
        free(c[i]);
    }
    free(b);
    free(c);
}

double find_max(double* arr, int n) {
    double max = arr[0];
    for (int i = 1; i < n; i++) {
        if (arr[i] > max) {
            max = arr[i]; 
        }
    }
    return max;
}

double* rk45(double* state_vector, double** b, double** c, double* h) {
    double input[8];
    double* new_state_vector = (double *)calloc(8, sizeof(double));
    double star_state_vector[8];
    double new_step;
    double error[8];
    double** k = malloc(6 * sizeof(double*));
    for (int i = 0; i < 6; i++) {
        k[i] = calloc(8, sizeof(double));
    }

    while (1) {
        memset(new_state_vector, 0, 8 * sizeof(double));
        memset(star_state_vector, 0, 8 * sizeof(double));
        double* tmp = eq_of_motion(state_vector);
        memcpy(k[0], tmp, 8 * sizeof(double));
        free(tmp);

        for (int k_it = 1; k_it < 6; k_it++) {
            for (int y_it = 0; y_it < 8; y_it++) {
                input[y_it] = state_vector[y_it];
                for (int b_it = 0; b_it < k_it; b_it++) {
                    input[y_it] += (*h) * b[k_it][b_it] * k[b_it][y_it];
                }
            }
            double* tmp = eq_of_motion(input);
            memcpy(k[k_it], tmp, 8 * sizeof(double));
            free(tmp);
        }

        for (int y_it = 0; y_it < 8; y_it++) {
            for (int c_it = 0; c_it < 6; c_it++) {
                new_state_vector[y_it] += k[c_it][y_it] * c[c_it][0];
                star_state_vector[y_it] += k[c_it][y_it] * c[c_it][1];
            }
            error[y_it] = new_state_vector[y_it] - star_state_vector[y_it];
        }

        new_step = (*h) * pow(acc / (find_max(error, 8) + 1e-10), 0.2);
        if (new_step >= (*h)) {
            (*h) = fmin(new_step, max_h);
            break;
        }
        (*h) = new_step;
    }

    for (int i = 0; i < 6; i++) {
        free(k[i]);
    }
    free(k);
    return new_state_vector;
}

double std(double* arr, int n) {
    if (n <= 0) return 0.0;
    double avg = 0;
    int N = n+1;
    for (int i = 0; i < N; i++) {
        avg += arr[i];
    }
    avg = avg/N;
    double res = 0;
    for (int i = 0; i < N; i++) {
        res += (arr[i] - avg)*(arr[i] - avg);
    }
    return sqrt(res/N);
}

void print_array(double *arr, int len) {
    for (int i = 0; i < len; i++)
        printf("%.5f, ", arr[i]);
    printf("\n");
}

double** make_trajectory(int N) {
    double** trajectory = malloc(N * sizeof(double*));
    for (int i = 0; i < N; i++) {
        trajectory[i] = malloc(8 * sizeof(double));
    }
    return trajectory;
}

void free_trajectory(double** trajectory, int N) {
    for (int i = 0; i < N; i++) {
        free(trajectory[i]);
    }
    free(trajectory);
}

double norm_vel(double* state_vector) {
    double** g_val = make_g(state_vector[2], state_vector[3]);
    double norm = + g_val[0][0] * state_vector[4] * state_vector[4]
                  + g_val[1][1] * state_vector[5] * state_vector[5] 
                  + 2*g_val[1][0] * state_vector[4] * state_vector[5] 
                  + g_val[2][2] * state_vector[6] * state_vector[6]
                  + g_val[3][3] * state_vector[7] * state_vector[7];
    free_g(g_val);
    return norm; 
}

double calculate_E(double* state_vector) {
    double** g_val = make_g(state_vector[2], state_vector[3]);
    double E = - g_val[0][0]*state_vector[4] - g_val[0][1]*state_vector[5];
    free_g(g_val);
    return E;
}

double calculate_L_z(double* state_vector) {
    double** g_val = make_g(state_vector[2], state_vector[3]);
    double L_z = g_val[1][1]*state_vector[5] + g_val[0][1]*state_vector[4];
    free_g(g_val);
    return L_z;
}

void write2csv(double** arr, int N) {
    FILE *ftpr;
    ftpr = fopen("trajectory.csv", "w");
    for (int i = 0; i < N; i++) {
        double* dpt = arr[i];
        fprintf(ftpr, "%.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f\n", dpt[0], dpt[1], dpt[2], dpt[3], dpt[4], dpt[5], dpt[6], dpt[7]);
    }
    fclose(ftpr);
}

int main() {
    const int N = (int)1e4;
    const int save_interval = 10;
    double** b = make_CashKarp_b();
    double** c = make_CashKarp_c();
    double** trajectory = make_trajectory(N/save_interval);
    double* arr_norm = (double*)malloc(N * sizeof(double));
    double* arr_E = (double*)malloc(N * sizeof(double));
    double* arr_L_z = (double*)malloc(N * sizeof(double));
    double* state_vector = (double*)calloc(8, sizeof(double));
    state_vector[2] = 7;
    state_vector = initialize_velocity(state_vector, 0.945, 3.29);
    print_array(state_vector, 8);

    for (int n = 0; n < N; n++) {
        double* new_state = rk45(state_vector, b, c, &h);
        free(state_vector);
        state_vector = new_state;

        arr_norm[n] = norm_vel(state_vector);
        arr_E[n] = calculate_E(state_vector);
        arr_L_z[n] = calculate_L_z(state_vector);

        if (n%save_interval == 0) {
            for (int j = 0; j < 8; j++) {
                trajectory[n/save_interval][j] = state_vector[j];
            }
            printf("Step %d/%d: normalization %.4f  E %.4f  L_z %.4f  std norm %.4f  std E %.4f  std L_z %.4f\n", 
                n, N, arr_norm[n], arr_E[n], arr_L_z[n], std(arr_norm, n), std(arr_E, n), std(arr_L_z, n));
        }
    }

    write2csv(trajectory, N/save_interval);

    free(state_vector);
    free(arr_norm);
    free(arr_E);
    free(arr_L_z);
    free_trajectory(trajectory, N/save_interval);
    free_CashKarp(b, c);

    return 0;
} 

