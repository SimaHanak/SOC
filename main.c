#include "functions.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>
#include <time.h>

const double L_z = 3.0;
const double E = 0.95;
//const double init_r = 5;
const double init_ur = 0;
const double M = 1.0;
const double J = 0.3;
const double alpha_const = 0.0;
const double beta_const = 0.2;
const double gamma_const = 3.0;

double h = 1e-2;
#define M_PI 3.14159265358979323846

void initialize_velocity(double state_vector[8], Params* p, double g[4][4], double g_inv[4][4]) {
    update_g(state_vector[2], state_vector[3], g, p);
    update_g_inv(state_vector[2], state_vector[3], g_inv, p);
    state_vector[4] = - g_inv[0][0]*E + g_inv[1][0]*L_z;
    state_vector[5] = g_inv[1][1]*L_z - g_inv[1][0]*E;
    state_vector[7] = sqrtl((- 1
                            - g[0][0] * state_vector[4] * state_vector[4]
                            - g[1][1] * state_vector[5] * state_vector[5] 
                            - 2*g[1][0] * state_vector[4] * state_vector[5] 
                            - g[2][2] * state_vector[6] * state_vector[6])/g[3][3]);
}

void print_array(double *arr, int len, char* text) {
    printf("%s ", text);
    for (int i = 0; i < len; i++)
        printf("%f, ", arr[i]);
    printf("\n");
}

void eq_of_motion(double state_vector[8], Params* p, double g[4][4], double g_inv[4][4], double dg[4][4][4], double Christoffel[4][4][4], double output[8]) {
    update_g(state_vector[2], state_vector[3], g, p);
    update_g_inv(state_vector[2], state_vector[3], g_inv, p);
    update_dg(state_vector[2], state_vector[3], dg, p);
    update_Christoffel_symbols(state_vector[2], state_vector[3], p, g, g_inv, dg, Christoffel);

    double vel[4];
    for (int i = 0; i < 4; i++) {
        output[i] = state_vector[i+4];
        vel[i] = state_vector[i+4];
        output[i+4] = 0.0;
    }

    for (int coords = 0; coords < 4; coords++) {
        for (int kappa = 0; kappa < 4; kappa++) {
            for (int lambda = 0; lambda < 4; lambda++) {
                output[coords + 4] -= Christoffel[coords][kappa][lambda] * vel[kappa] * vel[lambda];
            }
        }
    }
}

void rk4(double state_vector[8], Params* p, double g[4][4], double g_inv[4][4], double dg[4][4][4], double Christoffel[4][4][4]) {
    double k1[8], k2[8], k3[8], k4[8], input[8];
    
    eq_of_motion(state_vector, p, g, g_inv, dg, Christoffel, k1);
    for (int i = 0; i < 8; i++) {
        input[i] = state_vector[i] + h * k1[i]/2.0;
    }
    eq_of_motion(input, p, g, g_inv, dg, Christoffel, k2);
    for (int i = 0; i < 8; i++) {
        input[i] = state_vector[i] + h * k2[i]/2.0;
    }
    eq_of_motion(input, p, g, g_inv, dg, Christoffel, k3);
    for (int i = 0; i < 8; i++) {
        input[i] = state_vector[i] + h * k3[i];
    }
    eq_of_motion(input, p, g, g_inv, dg, Christoffel, k4);

    for (int i = 0; i < 8; i++) {
        state_vector[i] += h * (k1[i] + 2*k2[i] + 2*k3[i] + k4[i])/6.0;
    }

    update_g(state_vector[2], state_vector[3], g, p);
    update_g_inv(state_vector[2], state_vector[3], g_inv, p);
}

double calculate_E(double* state_vector, Params* p, double g[4][4]) {
    double E = - g[0][0]*state_vector[4] - g[0][1]*state_vector[5];
    return E;
}

double calculate_L_z(double* state_vector, Params* p, double g[4][4]) {
    double L_z = g[1][1]*state_vector[5] + g[0][1]*state_vector[4];
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

double norm_vel(double* state_vector, Params* p, double g[4][4]) {
    double norm = + g[0][0] * state_vector[4] * state_vector[4]
                  + g[1][1] * state_vector[5] * state_vector[5] 
                  + 2*g[1][0] * state_vector[4] * state_vector[5] 
                  + g[2][2] * state_vector[6] * state_vector[6]
                  + g[3][3] * state_vector[7] * state_vector[7];
    return norm; 
}

double sgn(double x) {
    if (x > 0) return 1;
    if (x < 0) return -1;
    return 0;
}

Params make_params(double M, double J, double alpha_const, double beta_const, double gamma_const) {
    double j = J/(M*M);
    Params p = {.M = M, .J = J, .M2 = -alpha_const*j*j*pow(M, 3), .S3 = -beta_const*pow(j, 3)*pow(M, 4), .M4 = gamma_const*pow(j, 4)*pow(M, 5)};
    return p;
}

void working_dir(char *folder_name, size_t size){
    time_t start_time;
    time(&start_time);
    strftime(folder_name, size, "%Y-%m-%d_%H-%M-%S", localtime(&start_time));
    printf("Creating new folder %s for this iteration...\n", folder_name);
    _mkdir(folder_name);
    strcat(folder_name, "/");
}

int main() {
    printf("Program started\n");

    printf("Initialization...\n");
    printf("\t Creating variables...\n");
    time_t start_time;
    time_t cur_time;
    time(&start_time);
    size_t save_interval = (int)1e5;

    int len_poincare_iteration = 10000;
    double r_start = 2.6883;
    double r_end = 2.6876;
    double r_step = 2e-5;
    double computation_count = (r_start - r_end)/r_step*len_poincare_iteration;
    int logged_all = 0;
    Params p = make_params(M, J, alpha_const, beta_const, gamma_const);
    //Params p = {.M = 1, .J = 0.33, .M2 = 0.28, .S3 = 0.05, .M4 = 0.01};

    printf("\t Opening files...\n");
    char folder_name[512];
    working_dir(folder_name, sizeof(folder_name));
    strftime(folder_name, sizeof(folder_name), "%Y-%m-%d_%H-%M-%S", localtime(&start_time));
    char metafilepath[512];
    snprintf(metafilepath, sizeof(metafilepath), "C:/Users/simon/Documents/01School/02SOC/SOC/%s/metadata.txt", folder_name);
    FILE *ftprmeta = fopen(metafilepath, "a");
    fprintf(ftprmeta, "#E%f,L_z%f,M%f,J%f,M2%f,S3%f,M4%f\n", E, L_z, p.M, p.J, p.M2, p.S3, p.M4);
    fclose(ftprmeta);
    printf("Initialization complete.\n\n");
    
    printf("Starting main loop...\n");
    int n_r = (int)round((r_start - r_end)/r_step);
    //int n_r = 1;
    printf("\t Starting %d processes...\n", n_r);
    #pragma omp parallel for schedule(dynamic)
    for (int idx = 0; idx < n_r; idx++) {
        printf("\t\t Starting initialization of process ID %d...\n", idx);
        double g[4][4] = {0};
        double g_inv[4][4] = {0};
        double dg[4][4][4] = {0};
        double Christoffel[4][4][4] = {0};
        double init_r = r_start - idx * r_step;

        char filepath[512];
        snprintf(filepath, sizeof(filepath), "%s/trajectory-%f.csv", folder_name, init_r);
        FILE *ftprtra = fopen(filepath, "a");
        fprintf(ftprtra, "init_r,r,ur\n");

        int logged_partial = 0;

        double state_vector[8] = {0};
        state_vector[2] = init_r;
        state_vector[6] = init_ur;
        initialize_velocity(state_vector, &p, g, g_inv);
        double prev_r = state_vector[2];
        double prev_z = state_vector[3];
        double prev_ur = state_vector[6];

        double norm_dev = fabs(norm_vel(state_vector, &p, g) + 1)/1;
        double E_dev = fabs(calculate_E(state_vector, &p, g) - E)/E;
        double L_z_dev = fabs(calculate_L_z(state_vector, &p, g) - L_z)/L_z;  
        // printf("norm: %e    E: %e   L_z: %e \n", norm_dev, E_dev, L_z_dev);

        //print_array(state_vector, 8, "State_vector: ");

        printf("\t\t Starting simulation of process ID %d...\n", idx);
        for (int n = 0; logged_partial < len_poincare_iteration; n++) {
        //for (int n = 0; n < 1e5; n++) {
            rk4(state_vector, &p, g, g_inv, dg, Christoffel);
            if (state_vector[2] < 1.0 || isnan(state_vector[2])) {
                break;
            }

            if ((sgn(prev_z) != sgn(state_vector[3])) && (sgn(state_vector[7]) == 1) && (n != 0)) {
                double r0 = (state_vector[3]*prev_r - prev_z*state_vector[2])/(state_vector[3] - prev_z);
                double ur0 = (state_vector[3]*prev_ur - prev_z*state_vector[6])/(state_vector[3] - prev_z);
                fprintf(ftprtra, "%f,%f,%f\n", init_r, r0, ur0);
                
                logged_partial++;
                #pragma omp atomic
                logged_all++;
            }
            //fprintf(ftprtra, "%f,%f,%f\n", state_vector[1], state_vector[2], state_vector[3]);

            if (n%save_interval == 0) {
                time(&cur_time);

                update_g(state_vector[2], state_vector[3], g, &p);
                update_g_inv(state_vector[2], state_vector[3], g_inv, &p);
                //norm_dev = - fabs(norm_vel(state_vector, &p, g) + 1)/1;
                //E_dev = fabs(calculate_E(state_vector, &p, g) - E)/E;
                //L_z_dev = fabs(calculate_L_z(state_vector, &p, g) - L_z)/L_z;
                double diff_time = difftime(cur_time, start_time);
                double time_per_step = diff_time/(double)logged_all;
                double predicted_time = time_per_step*computation_count - diff_time;
                predicted_time /= 3600;
                
                //printf("Step %e | Logged %d | Time %.4f | Ends in %.4f hours | Relative deviations:     norm: %e    E: %e   L_z: %e \n", (double)n, logged_partial, diff_time, predicted_time, norm_dev, E_dev, L_z_dev);
                //print_array(state_vector, 8, "State_vector: ");
                printf("Elapsed time: %.2f \t Predicted remaining time: %.2f h \t Process %d is %.2f done.\n", diff_time/3600, predicted_time, idx, (double)logged_partial/(double)len_poincare_iteration);
            }

            prev_z = state_vector[3];
            prev_r = state_vector[2];
            prev_ur = state_vector[6];
        }
        printf("\tThread ID %d finished successfully.\n", idx);
        fclose(ftprtra);
    }   
    printf("All threads finished successfully.\n");
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
