#include "functions.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>
#include <time.h>

#ifdef _WIN32
#include <direct.h>
#define MKDIR(x) _mkdir(x)
#else
#include <sys/stat.h>
#include <sys/types.h>
#define MKDIR(x) mkdir(x, 0755)
#endif

const double L_z = 3.0;
const double E = 0.95;
const double init_r = 5.0;
const double init_ur = 0;
const double M = 1.0;
const double J = 0.3;
const double alpha_const = 0.0;
const double beta_const = 0.2;
const double gamma_const = 3.0;

double h = 1e-2;

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

void initialize_velocity(double state_vector[16], Params* p, double g[4][4], double g_inv[4][4]) {
    update_g(state_vector[2], state_vector[3], g, p);
    update_g_inv(state_vector[2], state_vector[3], g_inv, p);
    state_vector[4] = - g_inv[0][0]*E + g_inv[1][0]*L_z;
    state_vector[5] = g_inv[1][1]*L_z - g_inv[1][0]*E;
    state_vector[7] = sqrt((- 1
                            - g[0][0] * state_vector[4] * state_vector[4]
                            - g[1][1] * state_vector[5] * state_vector[5] 
                            - 2*g[1][0] * state_vector[4] * state_vector[5] 
                            - g[2][2] * state_vector[6] * state_vector[6])/g[3][3]);
}

void print_arr(double arr[], int size) {
    for (int i = 0; i < size; i++) {
        printf("%f ", arr[i]);
    }
    printf("\n");
}

void print_NDIM(double *arr, int ndim, int dim, int offset) {
    if (dim == ndim - 1) {
        // Base case: last dimension
        for (int i = 0; i < 4; i++) {
            printf("%f ", arr[offset + i]);
        }
        printf("\n");
        return;
    }

    // Size of a block of the remaining dimensions
    int block = 1;
    for (int i = dim + 1; i < ndim; i++)
        block *= 4;

    // Iterate through this dimension
    for (int i = 0; i < 4; i++) {
        print_NDIM(arr, ndim, dim + 1, offset + i * block);
        if (dim == 0) printf("\n");
    }
}

void eq_of_motion(double state_vector[16], Params* p, double g[4][4], double g_inv[4][4], double dg[4][4][4], double dg_inv[4][4][4], double ddg[4][4][4][4], double Christoffel[4][4][4], double DChristoffel[4][4][4][4], double output[16]) {
    update_g(state_vector[2], state_vector[3], g, p);
    update_g_inv(state_vector[2], state_vector[3], g_inv, p);
    update_dg(state_vector[2], state_vector[3], dg, p);
    update_dg_inv(state_vector[2], state_vector[3], dg_inv, p);
    update_ddg(state_vector[2], state_vector[3], ddg, p);
    update_Christoffel_symbols(state_vector[2], state_vector[3], p, g, g_inv, dg, Christoffel);
    update_DChristoffel_symbols(state_vector[2], state_vector[3], p, g, g_inv, dg_inv, dg, ddg, DChristoffel);

    for (int i = 0; i < 4; i++) {
        output[i] = state_vector[i+4];
        output[i+8] = state_vector[i+12];
        for (int kappa = 0; kappa < 4; kappa++) {
            for (int lambda = 0; lambda < 4; lambda++) {
                output[i+4] -= Christoffel[i][kappa][lambda] * state_vector[kappa+4] * state_vector[lambda+4];
                for (int nu = 0; nu < 4; nu++) {
                    output[i+12] -= DChristoffel[i][kappa][lambda][nu]*state_vector[kappa+4]*state_vector[lambda+4]*state_vector[nu+8];
                    output[i+12] -= 2*Christoffel[i][kappa][lambda]*state_vector[kappa+4]*state_vector[lambda+12];
                }
            }
        }
    }
}

void rk4(double state_vector[16], Params* p, double g[4][4], double g_inv[4][4], double dg[4][4][4], double dg_inv[4][4][4], double ddg[4][4][4][4], double Christoffel[4][4][4], double DChristoffel[4][4][4][4]) {
    double k1[16] = {0}, k2[16] = {0}, k3[16] = {0}, k4[16] = {0}, input[16] = {0};
    eq_of_motion(state_vector, p, g, g_inv, dg, dg_inv, ddg, Christoffel, DChristoffel, k1);
    for (int i = 0; i < 16; i++) {
        input[i] = state_vector[i] + h * k1[i]/2.0;
    }
    eq_of_motion(input, p, g, g_inv, dg, dg_inv, ddg, Christoffel, DChristoffel, k2);
    for (int i = 0; i < 16; i++) {
        input[i] = state_vector[i] + h * k2[i]/2.0;
    }
    eq_of_motion(input, p, g, g_inv, dg, dg_inv, ddg, Christoffel, DChristoffel, k3);
    for (int i = 0; i < 16; i++) {
        input[i] = state_vector[i] + h * k3[i];
    }
    eq_of_motion(input, p, g, g_inv, dg, dg_inv, ddg, Christoffel, DChristoffel, k4);

    for (int i = 0; i < 16; i++) {
        state_vector[i] += h * (k1[i] + 2*k2[i] + 2*k3[i] + k4[i])/6.0;
    }

    update_g(state_vector[2], state_vector[3], g, p);
    update_g_inv(state_vector[2], state_vector[3], g_inv, p);
}

double calculate_E(double state_vector[16], Params* p, double g[4][4]) {
    double E = - g[0][0]*state_vector[4] - g[0][1]*state_vector[5];
    return E;
}

double calculate_L_z(double state_vector[16], Params* p, double g[4][4]) {
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

double norm_vel(double state_vector[16], Params* p, double g[4][4]) {
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
    if (MKDIR(folder_name) != 0) {
        perror("mkdir failed");
        exit(EXIT_FAILURE);
    }
    printf("Folder created.");
}

double compute_measure_of_dev(double state_vector[16], double g[4][4], Params *p) {
    double result = 0;
    for (int mu = 0; mu < 4; mu++) {
        for (int nu = 0; nu < 4; nu++) {
            result += g[mu][nu]*state_vector[mu+8]*state_vector[nu+8];
        }
    }
    return sqrt(fabs(result));
}

int main() {
    printf("Program started\n");

    printf("Initialization...\n");
    printf("\t Creating variables...\n");
    time_t start_time;
    time_t cur_time;
    time(&start_time);
    size_t save_interval = (int)1e3;

    double computation_count = 1e7;
    Params p = make_params(M, J, alpha_const, beta_const, gamma_const);
    //Params p = {.M = 1, .J = 0.33, .M2 = 0.28, .S3 = 0.05, .M4 = 0.01};

    printf("\t Opening files...\n");
    char folder_name[512];
    working_dir(folder_name, sizeof(folder_name));
    char metafilepath[1024];
    snprintf(metafilepath, sizeof(metafilepath), "C:/Users/simon/Documents/01School/02SOC/SOC/%s/metadata.txt", folder_name);
    //snprintf(metafilepath, sizeof(metafilepath), "/home/shanak/Documents/[01] Studium/SOČ/%s/metadata.txt", folder_name);
    FILE *ftprmeta = fopen(metafilepath, "a");
    if (!ftprmeta) {
        perror("fopen failed");
        fprintf(stderr, "Path: %s\n", metafilepath);
        exit(EXIT_FAILURE);
    }
    fprintf(ftprmeta, "#E%f,L_z%f,M%f,J%f,M2%f,S3%f,M4%f\n", E, L_z, p.M, p.J, p.M2, p.S3, p.M4);
    fclose(ftprmeta);
    printf("Initialization complete.\n\n");

    double g[4][4] = {0};
    double g_inv[4][4] = {0};
    double dg[4][4][4] = {0};
    double dg_inv[4][4][4] = {0};
    double ddg[4][4][4][4] = {0};
    double Christoffel[4][4][4] = {0};
    double DChristoffel[4][4][4][4] = {0};

    char filepath[1024];
    snprintf(filepath, sizeof(filepath), "%s/trajectory-%f.csv", folder_name, init_r);
    FILE *ftprtra = fopen(filepath, "a");
    fprintf(ftprtra, "t,sum_log_stretch\n");

    double state_vector[16] = {0};
    state_vector[2] = init_r;
    state_vector[6] = init_ur;
    state_vector[10] = 1.0;
    initialize_velocity(state_vector, &p, g, g_inv);
    double prev_r = state_vector[2];
    double prev_z = state_vector[3];
    double prev_ur = state_vector[6];

    double dev_log_sum = 0;

    double norm_dev = fabs(norm_vel(state_vector, &p, g) + 1)/1;
    double E_dev = fabs(calculate_E(state_vector, &p, g) - E)/E;
    double L_z_dev = fabs(calculate_L_z(state_vector, &p, g) - L_z)/L_z;
    
    printf("norm: %e    E: %e   L_z: %e \n", norm_dev, E_dev, L_z_dev);
    print_arr(state_vector, 16);

    //print_array(state_vector, 8, "State_vector: ");

    printf("Starting simulation...\n");
    for (int n = 0; n < computation_count; n++) {
        rk4(state_vector, &p, g, g_inv, dg, dg_inv, ddg, Christoffel, DChristoffel);
        if (state_vector[2] < 1.0 || isnan(state_vector[2]) || isnan(state_vector[10])) {
            break;
        }

        // if ((sgn(prev_z) != sgn(state_vector[3])) && (sgn(state_vector[7]) == 1) && (n != 0)) {
        //     double r0 = (state_vector[3]*prev_r - prev_z*state_vector[2])/(state_vector[3] - prev_z);
        //     double ur0 = (state_vector[3]*prev_ur - prev_z*state_vector[6])/(state_vector[3] - prev_z);
        //     fprintf(ftprtra, "%f,%f,%f\n", init_r, r0, ur0);
        // }
        //fprintf(ftprtra, "%f,%f,%f\n", state_vector[1], state_vector[2], state_vector[3]);

        if (n%200 == 0) {
            double measure_of_dev = compute_measure_of_dev(state_vector, g, &p);
            dev_log_sum += log(measure_of_dev);
            for (int i = 0; i < 8; i++) {
                state_vector[i+8] /= measure_of_dev;
            }
            fprintf(ftprtra, "%f,%f\n", state_vector[0], dev_log_sum);
        }

        if (n%save_interval == 0) {
            time(&cur_time);

            update_g(state_vector[2], state_vector[3], g, &p);
            norm_dev = - fabs(norm_vel(state_vector, &p, g) + 1);
            E_dev = fabs(calculate_E(state_vector, &p, g) - E)/E;
            L_z_dev = fabs(calculate_L_z(state_vector, &p, g) - L_z)/L_z;
            double diff_time = difftime(cur_time, start_time);
            double time_per_step = diff_time/(double)n;
            double predicted_time = time_per_step*computation_count - diff_time;
            predicted_time /= 3600;
            
            printf("Step %e | Time %.4f | Ends in %.4f hours | norm: %e    E: %e   L_z: %e \n", (double)n, diff_time, predicted_time, norm_dev, E_dev, L_z_dev);
            //print_array(state_vector, 8, "State_vector: ");
        }

        prev_z = state_vector[3];
        prev_r = state_vector[2];
        prev_ur = state_vector[6];
    }
    printf("Simulation finished successfully.\n");
    fclose(ftprtra);
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
