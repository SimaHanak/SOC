#include "functions.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

const double L_z = 3.0;
const double E = 0.95;
const double M = 1.0;
const double J = 0.3;
const double alpha_const = 0.0;
const double beta_const = 0.2;
const double gamma_const = 3.0;

Params make_params(double M, double J, double alpha_const, double beta_const, double gamma_const) {
    double j = J/(M*M);
    Params p = {.M = M, .J = J, .M2 = -alpha_const*j*j*pow(M, 3), .S3 = -beta_const*pow(j, 3)*pow(M, 4), .M4 = gamma_const*pow(j, 4)*pow(M, 5)};
    return p;
}

int main(){
    Params p = make_params(M, J, alpha_const, beta_const, gamma_const);
    double effective_potential;

    FILE *ftpr = fopen("V_eff_0.csv", "w");
    fprintf(ftpr, "r,z\n");

    for (double r=1; r < 15; r+=0.01) {
        for (double z=-6; z < 6; z+=0.01) {
            effective_potential = V_eff(r, z, E, L_z, &p);
            if (abs(effective_potential) < 1e-4) {
                fprintf(ftpr, "%f,%f\n", r, z);
            }
        }
    }
    return 0;
}