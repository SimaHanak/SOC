#include "functions.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

const double L_z = 3;
const double E = 0.95;
const double init_r = 5;
const double init_ur = 0;
const Params p = {.M = 1.0, .J = 0.3, .M2 = -0.1, .S3 = 0.05, .M4 = 0.01};

double h = 1e-4;
const double abs_tol = 1e-10;
const double rel_tol = 1e-8;
const double max_h = 1e-1;
const double min_h = 1e-8;
const double poincare_tol = 5e-5;
#define M_PI 3.14159265358979323846

static const double b[6][5] = {
    { 0.0,       0.0,        0.0,        0.0,       0.0 },
    { 1.0/5.0,   0.0,        0.0,        0.0,       0.0 },
    { 3.0/40.0,  9.0/40.0,   0.0,        0.0,       0.0 },
    { 3.0/10.0, -9.0/10.0,   6.0/5.0,    0.0,       0.0 },
    {-11.0/54.0, 5.0/2.0,   -70.0/27.0,  35.0/27.0, 0.0 },
    {1631.0/55296.0, 175.0/512.0, 575.0/13824.0, 44275.0/110592.0, 253.0/4096.0}
};

static const double c[6][2] = {
    { 37.0/378.0,     2825.0/27648.0 },
    { 0.0,            0.0 },
    { 250.0/621.0,    18575.0/48384.0 },
    { 125.0/594.0,    13525.0/55296.0 },
    { 0.0,            277.0/14336.0 },
    { 512.0/1771.0,   1.0/4.0 }
};

double* initialize_velocity(double* state_vector) {
    double** g_val = make_g(state_vector[2], 0, &p);
    double** g_inv_val = make_g_inv(state_vector[2], 0, &p);
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

void print_array(double *arr, int len, char* text) {
    printf("%s ", text);
    for (int i = 0; i < len; i++)
        printf("%.5f, ", arr[i]);
    printf("\n");
}

double* eq_of_motion(double* state_vector) {
    double*** Christoffel = generate_Christoffel_symbols(state_vector[2], state_vector[3], &p);

    double* dydt = (double*)calloc(8, sizeof(double));
    double vel[4];
    for (int i = 0; i < 4; i++) {
        dydt[i] = state_vector[i+4];
        vel[i] = state_vector[i+4];
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

double calculate_E(double* state_vector) {
    double** g_val = make_g(state_vector[2], state_vector[3], &p);
    double E = - g_val[0][0]*state_vector[4] - g_val[0][1]*state_vector[5];
    free_g(g_val);
    return E;
}

double calculate_L_z(double* state_vector) {
    double** g_val = make_g(state_vector[2], state_vector[3], &p);
    double L_z = g_val[1][1]*state_vector[5] + g_val[0][1]*state_vector[4];
    free_g(g_val);
    return L_z;
}

double find_max(double* arr, int n) {
    double max = fabs(arr[0]);
    for (int i = 1; i < n; i++) {
        double val = fabs(arr[i]);
        if (val > max) {
            max = val; 
        }
    }
    return max;
}

/*
double* rk45(double* state_vector, double* h_ptr) {
    double h = *h_ptr;
    double* new_state = calloc(8, sizeof(double));
    double* star_state = calloc(8, sizeof(double));
    double* error = malloc(8 * sizeof(double));
    double k[6][8];
    double temp[8];

    while (1) {
        // k1 = f(t, y)
        double* k_temp = eq_of_motion(state_vector);
        memcpy(k[0], k_temp, 8 * sizeof(double));
        free(k_tem&p);

        // Compute k2 to k6
        for (int i = 1; i < 6; i++) {
            for (int j = 0; j < 8; j++) {
                temp[j] = state_vector[j];
                for (int m = 0; m < i; m++) {
                    temp[j] += h * b[i][m] * k[m][j];
                }
            }
            k_temp = eq_of_motion(tem&p);
            memcpy(k[i], k_temp, 8 * sizeof(double));
            free(k_tem&p);
        }

        // Compute new state (4th and 5th order solutions)
        for (int i = 0; i < 8; i++) {
            new_state[i] = state_vector[i];
            star_state[i] = state_vector[i];
            for (int j = 0; j < 6; j++) {
                new_state[i] += h * c[j][0] * k[j][i];
                star_state[i] += h * c[j][1] * k[j][i];
            }
            // Error estimate (normalized)
            double scale = abs_tol + rel_tol * fmax(fabs(state_vector[i]), fabs(new_state[i]));
            error[i] = fabs(new_state[i] - star_state[i]) / scale;
        }

        double max_err = fmax(find_max(error, 8), 1e-8);

        // Compute adaptive step size
        double factor = 0.9 * pow(1.0 / (max_err + 1e-10), 0.2);
        factor = fmin(2.0, fmax(0.5, factor)); // Clamp factor

        double new_h = h * factor;

        double dL_z = fabs(calculate_L_z(new_state) - calculate_L_z(state_vector));
        double dE = fabs(calculate_E(new_state) - calculate_E(state_vector));
        if (dL_z > 1e-9 || dE > 1e-9) {
            h = h*0.1;
        } else if (max_err <= 1.0) {
            *h_ptr = fmin(fmax(new_h, min_h), max_h);
            break; // Accept step
        } else {
            if (new_h < min_h) {
                //fprintf(stderr, "Warning: step size too small, forcing min_h\n");
                *h_ptr = fmin(fmax(new_h, min_h), max_h);
                break;
            }
            h = new_h; // Retry with smaller h
        }
    }

    free(star_state);
    free(error);
    return new_state;
}
*/

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

double norm_vel(double* state_vector) {
    double** g_val = make_g(state_vector[2], state_vector[3], &p);
    double norm = + g_val[0][0] * state_vector[4] * state_vector[4]
                  + g_val[1][1] * state_vector[5] * state_vector[5] 
                  + 2*g_val[1][0] * state_vector[4] * state_vector[5] 
                  + g_val[2][2] * state_vector[6] * state_vector[6]
                  + g_val[3][3] * state_vector[7] * state_vector[7];
    free_g(g_val);
    return norm; 
}

int main() {
    FILE *ftpr;
    ftpr = fopen("trajectory.csv", "w");
    const unsigned long int N = (int)1e9;
    size_t save_interval = (int)1e3;
    double* arr_norm = (double*)malloc(N/save_interval * sizeof(double));
    double* arr_E = (double*)malloc(N/save_interval * sizeof(double));
    double* arr_L_z = (double*)malloc(N/save_interval * sizeof(double));
    double* state_vector = (double*)calloc(8, sizeof(double));
    state_vector[2] = init_r;
    state_vector[6] = init_ur;
    state_vector = initialize_velocity(state_vector);
    fprintf(ftpr, "E%f, L_z%f, r%f, ur%f, M%f, J%f, M2%f, S3%f, M4%f\n", E, L_z, init_r, init_ur, p.M, p.J, p.M2, p.S3, p.M4);
    print_array(state_vector, 8, "State_vector: ");
    for (int n = 0; n < N; n++) {
        //double* new_state = rk45(state_vector, &h);
        double* new_state = rk4(state_vector);
        free(state_vector);
        state_vector = new_state;
        if (state_vector[2] < 0.5) {
            return 0;
        }
        double modulo_f = fmod(state_vector[1], M_PI);
        if (modulo_f < poincare_tol || M_PI - poincare_tol < modulo_f) {
            fprintf(ftpr, "%.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f\n", 
                state_vector[0], state_vector[1], state_vector[2], state_vector[3], state_vector[4], state_vector[5], state_vector[6], state_vector[7]);
        }
        if (n%save_interval == 0) {
            int index = n/save_interval;
            arr_norm[index] = norm_vel(state_vector);
            arr_E[index] = calculate_E(state_vector);
            arr_L_z[index] = calculate_L_z(state_vector);

            double percentage = (double)n/(double)N * 100;
            printf("Step %d/%ld (done %.2f percent):  normalization %.4f  E %.4f  L_z %.4f  std norm %.4f  std E %.4f  std L_z %.4f  h %1.4e\n", 
                n+1, N, percentage, arr_norm[index], arr_E[index], arr_L_z[index], std(arr_norm, index), std(arr_E, index), std(arr_L_z, index), h);

            fprintf(ftpr, "%.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f\n", 
                state_vector[0], state_vector[1], state_vector[2], state_vector[3], state_vector[4], state_vector[5], state_vector[6], state_vector[7]);
        }
    }

    free(state_vector);
    free(arr_norm);
    free(arr_E);
    free(arr_L_z);
    fclose(ftpr);

    return 0;
} 

/*
int main_commented() {
    double* state_vector = (double*)calloc(8, sizeof(double));
    state_vector[2] = 6;
    state_vector[6] = 0.01;
    state_vector = initialize_velocity(state_vector);
    for (int i = 0; i < 30; i++) {
        printf("Step %d\n", i);
        print_array(state_vector, 8, "State_vector:");
        print_array(eq_of_motion(state_vector), 8, "EoM(State_vector):");
        printf("norm %f  E %f  L_z %f\n", norm_vel(state_vector), calculate_E(state_vector), calculate_L_z(state_vector));
        printf("h = %f\n\n", h);
        state_vector = rk45(state_vector, &h);
    }
}
*/