#include "functions.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
//#include <windows.h>
#include <time.h>

const long double L_z = 3L;
const long double E = 0.95L;
//const long double init_r = 5;
const long double init_ur = 0L;
const long double M = 1.0L;
const long double J = 0.3L;
const long double alpha_const = 1.5L;
const long double beta_const = 1L;
const long double gamma_const = 1L;

long double h = 1e-6L;
#define M_PI 3.14159265358979323846

static const long double b[6][5] = {
    { 0.0,       0.0,        0.0,        0.0,       0.0 },
    { 1.0/5.0,   0.0,        0.0,        0.0,       0.0 },
    { 3.0/40.0,  9.0/40.0,   0.0,        0.0,       0.0 },
    { 3.0/10.0, -9.0/10.0,   6.0/5.0,    0.0,       0.0 },
    {-11.0/54.0, 5.0/2.0,   -70.0/27.0,  35.0/27.0, 0.0 },
    {1631.0/55296.0, 175.0/512.0, 575.0/13824.0, 44275.0/110592.0, 253.0/4096.0}
};

static const long double c[6][2] = {
    { 37.0/378.0,     2825.0/27648.0 },
    { 0.0,            0.0 },
    { 250.0/621.0,    18575.0/48384.0 },
    { 125.0/594.0,    13525.0/55296.0 },
    { 0.0,            277.0/14336.0 },
    { 512.0/1771.0,   1.0/4.0 }
};

long double* initialize_velocity(long double* state_vector, Params* p, long double** g, long double** g_inv, long double*** dg) {
    update_g(state_vector[2], state_vector[3], g, p);
    update_g_inv(state_vector[2], state_vector[3], g_inv, p);
    state_vector[4] = - g_inv[0][0]*E + g_inv[1][0]*L_z;
    state_vector[5] = g_inv[1][1]*L_z - g_inv[1][0]*E;
    state_vector[7] = sqrtl((- 1
                            - g[0][0] * state_vector[4] * state_vector[4]
                            - g[1][1] * state_vector[5] * state_vector[5] 
                            - 2*g[1][0] * state_vector[4] * state_vector[5] 
                            - g[2][2] * state_vector[6] * state_vector[6])/g[3][3]);

    return state_vector;
}

void print_array(long double *arr, int len, char* text) {
    printf("%s ", text);
    for (int i = 0; i < len; i++)
        printf("%Lf, ", arr[i]);
    printf("\n");
}

long double* eq_of_motion(long double* state_vector, Params* p, long double** g, long double** g_inv, long double*** dg) {
    update_g(state_vector[2], state_vector[3], g, p);
    update_dg(state_vector[2], state_vector[3], dg, p);
    update_g_inv(state_vector[2], state_vector[3], g_inv, p);
    long double*** Christoffel = generate_Christoffel_symbols(state_vector[2], state_vector[3], p, g, g_inv, dg);

    long double* dydt = (long double*)calloc(8, sizeof(long double));
    long double vel[4];
    for (int i = 0; i < 4; i++) {
        dydt[i] = state_vector[i+4];
        vel[i] = state_vector[i+4];
    }

    for (int coords = 0; coords < 4; coords++) {
        for (int kappa = 0; kappa < 4; kappa++) {
            for (int lambda = 0; lambda < 4; lambda++) {
                dydt[coords + 4] -= Christoffel[coords][kappa][lambda] * vel[kappa] * vel[lambda];
            }
        }
    }
    free_Christoffel(Christoffel);

    return dydt;
}

long double* rk4(long double* state_vector, Params* p, long double** g, long double** g_inv, long double*** dg) {
    long double input[8];
    long double* k1 = eq_of_motion(state_vector, p, g, g_inv, dg);
    for (int i = 0; i < 8; i++) {
        input[i] = state_vector[i] + h * k1[i]/2.0;
    }
    long double* k2 = eq_of_motion(input, p, g, g_inv, dg);
    for (int i = 0; i < 8; i++) {
        input[i] = state_vector[i] + h * k2[i]/2.0;
    }
    long double* k3 = eq_of_motion(input, p, g, g_inv, dg);
    for (int i = 0; i < 8; i++) {
        input[i] = state_vector[i] + h * k3[i];
    }
    long double* k4 = eq_of_motion(input, p, g, g_inv, dg);

    long double* result = (long double*)malloc(8 * sizeof(long double));
    for (int i = 0; i < 8; i++) {
        result[i] = state_vector[i] + h * (k1[i] + 2*k2[i] + 2*k3[i] + k4[i])/6.0;
    }

    free(k1);
    free(k2);
    free(k3);
    free(k4);
    update_g(state_vector[2], state_vector[3], g, p);
    update_g_inv(state_vector[2], state_vector[3], g_inv, p);
    return result;
}

long double calculate_E(long double* state_vector, Params* p, long double** g) {
    long double E = - g[0][0]*state_vector[4] - g[0][1]*state_vector[5];
    return E;
}

long double calculate_L_z(long double* state_vector, Params* p, long double** g) {
    long double L_z = g[1][1]*state_vector[5] + g[0][1]*state_vector[4];
    return L_z;
}

long double find_max(long double* arr, int n) {
    long double max = fabs(arr[0]);
    for (int i = 1; i < n; i++) {
        long double val = fabs(arr[i]);
        if (val > max) {
            max = val; 
        }
    }
    return max;
}

/*
long double* rk45(long double* state_vector, long double* h_ptr) {
    long double h = *h_ptr;
    long double* new_state = calloc(8, sizeof(long double));
    long double* star_state = calloc(8, sizeof(long double));
    long double* error = malloc(8 * sizeof(long double));
    long double k[6][8];
    long double temp[8];

    while (1) {
        // k1 = f(t, y)
        long double* k_temp = eq_of_motion(state_vector);
        memcpy(k[0], k_temp, 8 * sizeof(long double));
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
            memcpy(k[i], k_temp, 8 * sizeof(long double));
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
            long double scale = abs_tol + rel_tol * fmax(fabs(state_vector[i]), fabs(new_state[i]));
            error[i] = fabs(new_state[i] - star_state[i]) / scale;
        }

        long double max_err = fmax(find_max(error, 8), 1e-8);

        // Compute adaptive step size
        long double factor = 0.9 * pow(1.0 / (max_err + 1e-10), 0.2);
        factor = fmin(2.0, fmax(0.5, factor)); // Clamp factor

        long double new_h = h * factor;

        long double dL_z = fabs(calculate_L_z(new_state) - calculate_L_z(state_vector));
        long double dE = fabs(calculate_E(new_state) - calculate_E(state_vector));
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

long double norm_vel(long double* state_vector, Params* p, long double** g) {
    long double norm = + g[0][0] * state_vector[4] * state_vector[4]
                  + g[1][1] * state_vector[5] * state_vector[5] 
                  + 2*g[1][0] * state_vector[4] * state_vector[5] 
                  + g[2][2] * state_vector[6] * state_vector[6]
                  + g[3][3] * state_vector[7] * state_vector[7];
    return norm; 
}

long double sgn(long double x) {
    if (x > 0) return 1;
    if (x < 0) return -1;
    return 0;
}

Params make_params(long double M, long double J, long double alpha_const, long double beta_const, long double gamma_const) {
    long double j = J/(M*M);
    Params p = {.M = M, .J = J, .M2 = -alpha_const*j*j*pow(M, 3), .S3 = -beta_const*pow(j, 3)*pow(M, 4), .M4 = gamma_const*pow(j, 4)*pow(M, 5)};
    return p;
}

int main() {
    //SetThreadExecutionState(ES_CONTINUOUS | ES_DISPLAY_REQUIRED);
    time_t start_time;
    time_t cur_time;
    time(&start_time);
    const unsigned long N = 1e9;
    size_t save_interval = (int)1e4;

    long double** g = make_g();
    long double** g_inv = make_g_inv();
    long double*** dg = make_dg();

    FILE *ftpr;
    ftpr = fopen("trajectory.csv", "a");

    for (long double init_r = 10.5; init_r > 2.0; init_r -= 0.5) {
        int logged = 0;

        long double* state_vector = (long double*)calloc(8, sizeof(long double));
        state_vector[2] = init_r;
        state_vector[6] = init_ur;

        Params p = make_params(M, J, alpha_const, beta_const, gamma_const);
        state_vector = initialize_velocity(state_vector, &p, g, g_inv, dg);
        long double prev_z = state_vector[3];

        long double norm_dev = -fabs(norm_vel(state_vector, &p, g) + 1)/1;
        long double E_dev = fabs(calculate_E(state_vector, &p, g) - E)/E;
        long double L_z_dev = fabs(calculate_L_z(state_vector, &p, g) - L_z)/L_z;  
        printf("norm: %Lf    E: %Le   L_z: %Lf \n", norm_dev, E_dev, L_z_dev);

        fprintf(ftpr, "E%Lf, L_z%Lf, r%Lf, ur%Lf, M%Lf, J%Lf, M2%Lf, S3%Lf, M4%Lf\n", E, L_z, init_r, init_ur, p.M, p.J, p.M2, p.S3, p.M4);
        print_array(state_vector, 8, "State_vector: ");

        for (int n = 0; logged < 1000; n++) {
            //long double* new_state = rk45(state_vector, &h);
            long double* new_state = rk4(state_vector, &p, g, g_inv, dg);
            free(state_vector);
            state_vector = new_state;
            if (state_vector[2] < 0.5) {
                break;
            }

            if ((sgn(prev_z) != sgn(state_vector[3])) && (sgn(state_vector[7]) == 1)) {
                logged += 1;
                //printf("Logging...\n");
                fprintf(ftpr, "%Lf, %Lf, %Lf, %Lf, %Lf, %Lf, %Lf, %Lf\n", 
                    state_vector[0], state_vector[1], state_vector[2], state_vector[3], state_vector[4], state_vector[5], state_vector[6], state_vector[7]);
            }

            if (n%save_interval == 0) {
                time(&cur_time);

                norm_dev = - fabs(norm_vel(state_vector, &p, g) + 1)/1;
                E_dev = fabs(calculate_E(state_vector, &p, g) - E)/E;
                L_z_dev = fabs(calculate_L_z(state_vector, &p, g) - L_z)/L_z;
                
                printf("Step %Lf | Logged %d | Time %.4f | Relative deviations:     norm: %Lf    E: %Lf   L_z: %Lf \n", (long double)n, logged, difftime(cur_time, start_time), norm_dev, E_dev, L_z_dev);
            }
            prev_z = state_vector[3];
        }
    }

    free_g(g);
    free_g_inv(g_inv);
    free_dg(dg);
    fclose(ftpr);

    return 0;
} 

/*
int main_commented() {
    long double* state_vector = (long double*)calloc(8, sizeof(long double));
    state_vector[2] = 6;
    state_vector[6] = 0.01;
    state_vector = initialize_velocity(state_vector);
    for (int i = 0; i < 30; i++) {
        printf("Step %d\n", i);
        print_array(state_vector, 8, "State_vector:");
        print_array(eq_of_motion(state_vector), 8, "EoM(State_vector):");
        printf("norm %Lf  E %Lf  L_z %Lf\n", norm_vel(state_vector), calculate_E(state_vector), calculate_L_z(state_vector));
        printf("h = %Lf\n\n", h);
        state_vector = rk45(state_vector, &h);
    }
}
*/