#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "functions.h"

double pythagorean(double r, double z, Params *p) {
    return r*r + z*z;
}

double A(double r, double z, Params *p) {
    double term1 = 8*r*r*z*z*(24*(*p).J*(*p).J*(*p).M + 17*(*p).M*(*p).M*(*p).M2 + 21*(*p).M4);
    double term2 = powl(r, 4)*(-10*(*p).J*(*p).J*(*p).M + 7*powl((*p).M, 5) + 32*(*p).M2*(*p).M*(*p).M - 21*(*p).M4);
    double term3 = 8*powl(z, 4)*(20*(*p).J*(*p).J*(*p).M - 7*powl((*p).M, 5) - 22*(*p).M2*(*p).M*(*p).M - 7*(*p).M4);
    return term1 + term2 + term3;
}

double dA_r(double r, double z, Params *p) {
    double term1 = 16*r*z*z*(24*(*p).J*(*p).J*(*p).M + 17*(*p).M*(*p).M*(*p).M2 + 21*(*p).M4);
    double term2 = 4*powl(r, 3)*(-10*(*p).J*(*p).J*(*p).M + 7*powl((*p).M, 5) + 32*(*p).M2*(*p).M*(*p).M - 21*(*p).M4);
    double term3 = 0;
    return term1 + term2 + term3;
}

double dA_z(double r, double z, Params *p) {
    double term1 = 16*r*r*z*(24*(*p).J*(*p).J*(*p).M + 17*(*p).M*(*p).M*(*p).M2 + 21*(*p).M4);
    double term2 = 0;
    double term3 = 32*powl(z, 3)*(20*(*p).J*(*p).J*(*p).M - 7*powl((*p).M, 5) - 22*(*p).M2*(*p).M*(*p).M - 7*(*p).M4);
    return term1 + term2 + term3;
}

double B(double r, double z, Params *p) {
    double term1 = powl(r, 4)*(10*(*p).J*(*p).J + 10*(*p).M2*powl((*p).M, 3) + 21*(*p).M4*(*p).M + 7*(*p).M2*(*p).M2);
    double term2 = 4*powl(z, 4)*(-40*(*p).J*(*p).J*(*p).M*(*p).M - 14*(*p).J*(*p).S3 + 7*powl((*p).M, 6) + 30*(*p).M2*powl((*p).M, 3) + 14*(*p).M4*(*p).M + 7*(*p).M2*(*p).M2);
    double term3 = - 4*r*r*z*z*(27*(*p).J*(*p).J*(*p).M*(*p).M - 21*(*p).J*(*p).S3 + 7*powl((*p).M, 6) + 48*(*p).M2*powl((*p).M, 3) + 42*(*p).M4*(*p).M + 7*(*p).M2*(*p).M2);
    return term1 + term2 + term3;
}

double dB_r(double r, double z, Params *p) {
    double term1 = 4*powl(r, 3)*(10*(*p).J*(*p).J + 10*(*p).M2*powl((*p).M, 3) + 21*(*p).M4*(*p).M + 7*(*p).M2*(*p).M2);
    double term2 = 0;
    double term3 = - 8*r*z*z*(27*(*p).J*(*p).J*(*p).M*(*p).M - 21*(*p).J*(*p).S3 + 7*powl((*p).M, 6) + 48*(*p).M2*powl((*p).M, 3) + 42*(*p).M4*(*p).M + 7*(*p).M2*(*p).M2);
    return term1 + term2 + term3;
}

double dB_z(double r, double z, Params *p) {
    double term1 = 0;
    double term2 = 16*powl(z, 3)*(-40*(*p).J*(*p).J*(*p).M*(*p).M - 14*(*p).J*(*p).S3 + 7*powl((*p).M, 6) + 30*(*p).M2*powl((*p).M, 3) + 14*(*p).M4*(*p).M + 7*(*p).M2*(*p).M2);
    double term3 = - 8*r*r*z*(27*(*p).J*(*p).J*(*p).M*(*p).M - 21*(*p).J*(*p).S3 + 7*powl((*p).M, 6) + 48*(*p).M2*powl((*p).M, 3) + 42*(*p).M4*(*p).M + 7*(*p).M2*(*p).M2);
    return term1 + term2 + term3;
}

double H(double r, double z, Params *p) {
    double term1 = 4*r*r*z*z*((*p).J*((*p).M2 - 2*powl((*p).M, 3)) - 3*(*p).M*(*p).S3);
    double term2 = powl(r, 4)*((*p).J*(*p).M2 + 3*(*p).M*(*p).S3);
    return term1 + term2;
}

double dH_r(double r, double z, Params *p) {
    double term1 = 8*r*z*z*((*p).J*((*p).M2 - 2*powl((*p).M, 3)) - 3*(*p).M*(*p).S3);
    double term2 = 4*powl(r, 3)*((*p).J*(*p).M2 + 3*(*p).M*(*p).S3);
    return term1 + term2;
}

double dH_z(double r, double z, Params *p) {
    double term1 = 8*r*r*z*((*p).J*((*p).M2 - 2*powl((*p).M, 3)) - 3*(*p).M*(*p).S3);
    double term2 = 0;
    return term1 + term2;
}

double G(double r, double z, Params *p) {
    double term1 = powl((*p).J, 3)*(- powl(r, 4)
                                - 8*powl(z, 4)
                                + 12*r*r*z*z);

    double term2 = (*p).J*(*p).M*((+ powl((*p).M, 3) + 2*(*p).M2)*powl(r, 4) 
                             - 8*(3*powl((*p).M, 3) + 2*(*p).M2)*powl(z, 4)
                             + 4*(powl((*p).M, 3) + 10*(*p).M2)*r*r*z*z);

    double term3 = (*p).M*(*p).M*(*p).S3*(+ 3*powl(r, 4)
                                 - 40*powl(z, 4)
                                 + 12*r*r*z*z);
    return r*r*(term1 + term2 + term3);
}

double dG_r(double r, double z, Params *p) {
    double term1 = powl((*p).J, 3)*(- powl(r, 4)
                                - 8*powl(z, 4)
                                + 12*r*r*z*z);

    double term2 = (*p).J*(*p).M*(+ (powl((*p).M, 3) + 2*(*p).M2)*powl(r, 4) 
                            - 8*(3*powl((*p).M, 3) + 2*(*p).M2)*powl(z, 4)
                            + 4*(powl((*p).M, 3) + 10*(*p).M2)*r*r*z*z);

    double term3 = (*p).M*(*p).M*(*p).S3*(+ 3*powl(r, 4)
                                 - 40*powl(z, 4)
                                 + 12*r*r*z*z);

    double left = r*r;
    double dright = + powl((*p).J, 3)*(- 4*powl(r, 3)
                               - 0
                               + 24*r*z*z)
                    + (*p).J*(*p).M*(+ (powl((*p).M, 3) + 2*(*p).M2)*4*powl(r, 3)
                               - 0
                               + 4*(powl((*p).M, 3) + 10*(*p).M2)*2*r*z*z)
                    + (*p).M*(*p).M*(*p).S3*(+ 12*powl(r, 3)
                                    - 0
                                    + 24*r*z*z);
    double right = term1 + term2 + term3;
    double dleft = 2*r;
    return left*dright + right*dleft;
}

double dG_z(double r, double z, Params *p) {
    double term1 = powl((*p).J, 3)*(- 0
                                - 32*powl(z, 3)
                                + 24*r*r*z);

    double term2 = (*p).J*(*p).M*(+ 0 
                            - 32*(3*powl((*p).M, 3) + 2*(*p).M2)*powl(z, 3)
                            + 8*(powl((*p).M, 3) + 10*(*p).M2)*r*r*z);

    double term3 = (*p).M*(*p).M*(*p).S3*(+ 0
                                 - 160*powl(z, 3)
                                 + 24*r*r*z);
    
    return r*r*(term1 + term2 + term3);
}

double F(double r, double z, Params *p) {
    return (+ powl(r, 4)*((*p).S3 - (*p).J*(*p).M*(*p).M)
            - 4*r*r*z*z*((*p).J*(*p).M*(*p).M + (*p).S3));
}

double dF_r(double r, double z, Params *p) {
    return (+ 4*powl(r, 3)*((*p).S3 - (*p).J*(*p).M*(*p).M)
            - 8*r*z*z*((*p).J*(*p).M*(*p).M + (*p).S3)); 
}

double dF_z(double r, double z, Params *p) {
    return (+ 0
            - 8*r*r*z*((*p).J*(*p).M*(*p).M + (*p).S3)); 
}

double f(double r, double z, Params *p) {
    double pyth = pythagorean(r, z, p);
    double term1 = 1;
    double term2 = - (2*(*p).M)/(sqrt(pyth));
    double term3 = (2*(*p).M*(*p).M)/(pyth);
    double term4 = (((*p).M2 - powl((*p).M, 3))*r*r - 2*(powl((*p).M, 3) + (*p).M2)*z*z)/(powl(pyth, 5.0/2.0));
    double term5 = (2*z*z*(-(*p).J*(*p).J + powl((*p).M, 4) + 2*(*p).M2*(*p).M) - 2*(*p).M*(*p).M2*r*r)/(powl(pyth, 3));
    double term6 = (A(r, z, p))/(28*powl(pyth, 9.0/2.0));
    double term7 = (B(r, z, p))/(14*powl(pyth, 5));
    return term1 + term2 + term3 + term4 + term5 + term6 + term7;
}

double df_r(double r, double z, Params *p) {
    double pyth = pythagorean(r, z, p);
    double term1 = 0;
    double term2 = (2*r*(*p).M)/(powl(pyth, 3.0/2.0));
    double term3 = - (4*r*(*p).M*(*p).M)/(pyth*pyth);
    double term4 = (2*r*powl(pyth, 5.0/2.0)*((*p).M2 - powl((*p).M, 3)) - 5*r*powl(pyth, 3.0/2.0)*(((*p).M2 - powl((*p).M, 3))*r*r - 2*(powl((*p).M, 3) + (*p).M2)*z*z))/(powl(pyth, 5));
    double term5 = (-4*r*powl(pyth, 3)*(*p).M*(*p).M2 - 6*r*powl(pyth, 2)*(2*z*z*(-(*p).J*(*p).J + powl((*p).M, 4) + 2*(*p).M2*(*p).M) - 2*(*p).M*(*p).M2*r*r))/(powl(pyth, 6));
    double term6 = (powl(pyth, 9.0/2.0)*dA_r(r, z, p) - 9*r*powl(pyth, 7.0/2.0)*A(r, z, p))/(28*powl(pyth, 9));
    double term7 = (powl(pyth, 5)*dB_r(r, z, p) - 10*r*powl(pyth, 4)*B(r, z, p))/(14*powl(pyth, 10));
    return term1 + term2 + term3 + term4 + term5 + term6 + term7;
}

double df_z(double r, double z, Params *p) {
    double pyth = pythagorean(r, z, p);
    double term1 = 0;
    double term2 = (2*z*(*p).M)/(powl(pyth, 3.0/2.0));
    double term3 = - (4*z*(*p).M*(*p).M)/(pyth*pyth);
    double term4 = (-4*z*powl(pyth, 5.0/2.0)*(powl((*p).M, 3) + (*p).M2) - 5*z*powl(pyth, 3.0/2.0)*(((*p).M2 - powl((*p).M, 3))*r*r - 2*(powl((*p).M, 3) + (*p).M2)*z*z))/(powl(pyth, 5));
    double term5 = (4*z*powl(pyth, 3)*(-(*p).J*(*p).J + powl((*p).M, 4) + 2*(*p).M2*(*p).M) - 6*powl(pyth, 2)*z*(2*z*z*(-(*p).J*(*p).J + powl((*p).M, 4) + 2*(*p).M2*(*p).M) - 2*(*p).M*(*p).M2*r*r))/(powl(pyth, 6));
    double term6 = (powl(pyth, 9.0/2.0)*dA_z(r, z, p) - 9*z*powl(pyth, 7.0/2.0)*A(r, z, p))/(28*powl(pyth, 9));
    double term7 = (powl(pyth, 5)*dB_z(r, z, p) - 10*z*powl(pyth, 4)*B(r, z, p))/(14*powl(pyth, 10));
    return term1 + term2 + term3 + term4 + term5 + term6 + term7;
}

double omega(double r, double z, Params *p) {
    double pyth = pythagorean(r, z, p);
    double term1 = - (2*(*p).J*r*r)/(powl(pyth, 3.0/2.0));
    double term2 = - (2*(*p).J*(*p).M*r*r)/(pyth*pyth);
    double term3 = + (F(r, z, p))/(powl(pyth, 7.0/2.0));
    double term4 = + (H(r, z, p))/(2*powl(pyth, 4));
    double term5 = + (G(r, z, p))/(4*powl(pyth, 11.0/2.0));
    return term1 + term2 + term3 + term4 + term5;
}

double domega_r(double r, double z, Params *p) {
    double pyth = pythagorean(r, z, p);
    double term1 = - (4*r*powl(pyth, 3.0/2.0)*(*p).J - 6*powl(r, 3)*powl(pyth, 1.0/2.0)*(*p).J)/(powl(pyth, 3));
    double term2 = - (4*r*pyth*pyth*(*p).J*(*p).M - 8*powl(r, 3)*pyth*(*p).J*(*p).M)/(powl(pyth, 4));
    double term3 = + (dF_r(r, z, p)*powl(pyth, 7.0/2.0) - 7*r*powl(pyth, 5.0/2.0)*F(r, z, p))/(powl(pyth, 7));
    double term4 = + (2*dH_r(r, z, p)*powl(pyth, 4) - 16*r*powl(pyth, 3)*H(r, z, p))/(4*powl(pyth, 8));
    double term5 = + (4*dG_r(r, z, p)*powl(pyth, 11.0/2.0) - 44*r*powl(pyth, 9.0/2.0)*G(r, z, p))/(16*powl(pyth, 11));
    return term1 + term2 + term3 + term4 + term5;
}

double domega_z(double r, double z, Params *p) {
    double pyth = pythagorean(r, z, p);
    double term1 = + (6*r*r*z*powl(pyth, 1.0/2.0)*(*p).J)/(powl(pyth, 3));
    double term2 = + (8*r*r*z*pyth*(*p).J*(*p).M)/(powl(pyth, 4));
    double term3 = + (dF_z(r, z, p)*powl(pyth, 7.0/2.0) - 7*z*powl(pyth, 5.0/2.0)*F(r, z, p))/(powl(pyth, 7));
    double term4 = + (2*dH_z(r, z, p)*powl(pyth, 4) - 16*z*powl(pyth, 3)*H(r, z, p))/(4*powl(pyth, 8));
    double term5 = + (4*dG_z(r, z, p)*powl(pyth, 11.0/2.0) - 44*z*powl(pyth, 9.0/2.0)*G(r, z, p))/(16*powl(pyth, 11));
    return term1 + term2 + term3 + term4 + term5;
}

double my_gamma(double r, double z, Params *p) {
    double pyth = pythagorean(r, z, p);
    double term1 = + (r*r*((*p).J*(*p).J*(r*r - 8*z*z) + (*p).M*(powl((*p).M, 3) + 3*(*p).M2)*(r*r - 4*z*z)))/(4*powl(pyth, 4));
    double term2 = - ((*p).M*(*p).M*r*r)/(2*pyth*pyth);
    return term1 + term2;
}

double dgamma_r(double r, double z, Params *p) {
    double pyth = pythagorean(r, z, p);
    double left = r*r;
    double dright = (*p).J*(*p).J*2*r + (*p).M*(powl((*p).M, 3) + 3*(*p).M2)*2*r;
    double right = (*p).J*(*p).J*(r*r - 8*z*z) + (*p).M*(powl((*p).M, 3) + 3*(*p).M2)*(r*r - 4*z*z);
    double dleft = 2*r;
    double denom = 4*powl(pyth, 4);
    double term1 = + ((left*dright + right*dleft)*denom - left*right*32*r*powl(pyth, 3))/(denom*denom);
    double term2 = - (4*r*pyth*pyth*(*p).M*(*p).M - 8*powl(r, 3)*pyth*(*p).M*(*p).M)/(denom);
    return term1 + term2;
}

double dgamma_z(double r, double z, Params *p) {
    double pyth = pythagorean(r, z, p);
    double term1 = + (4*r*r*powl(pyth, 4)*(-(*p).J*(*p).J*16*z - (*p).M*(powl((*p).M, 3) + 3*(*p).M2)*8*z) - 32*r*r*z*powl(pyth, 3)*((*p).J*(*p).J*(r*r - 8*z*z) + (*p).M*(powl((*p).M, 3) + 3*(*p).M2)*(r*r - 4*z*z)))/(16*powl(pyth, 8));
    double term2 = + (8*r*r*z*pyth)/(4*powl(pyth, 4));
    return term1 + term2;
}

double g_tt(double r, double z, Params *p) { return - f(r, z, p); }
double dg_tt_r(double r, double z, Params *p) { return - df_r(r, z, p); }
double dg_tt_z(double r, double z, Params *p) { return - df_z(r, z, p); }

double g_tf(double r, double z, Params *p) { return omega(r, z, p)*f(r, z, p); }
double dg_tf_r(double r, double z, Params *p) { return omega(r, z, p)*df_r(r, z, p) + domega_r(r, z, p)*f(r, z, p); }
double dg_tf_z(double r, double z, Params *p) { return omega(r, z, p)*df_z(r, z, p) + domega_z(r, z, p)*f(r, z, p); }

double g_ff(double r, double z, Params *p) { return - f(r, z, p)*powl(omega(r, z, p), 2) + (r*r)/f(r, z, p); }
double dg_ff_r(double r, double z, Params *p) {
    double f_val = f(r, z, p);
    double omega_val = omega(r, z, p);
    return - 2*f_val*omega_val*domega_r(r, z, p) - df_r(r, z, p)*omega_val*omega_val + (2*r*f_val - r*r*df_r(r, z, p))/(f_val*f_val); 
}
double dg_ff_z(double r, double z, Params *p) {
    double f_val = f(r, z, p);
    double omega_val = omega(r, z, p);
    return - 2*f_val*omega_val*domega_z(r, z, p) - df_z(r, z, p)*omega_val*omega_val - (r*r*df_z(r, z, p))/(f_val*f_val);
}

double g_rr(double r, double z, Params *p) { return exp(2*my_gamma(r, z, p))/f(r, z, p); }
double dg_rr_r(double r, double z, Params *p) { 
    double f_val = f(r, z, p);
    double e_2gamma = exp(2*my_gamma(r, z, p));
    return (f_val*e_2gamma*2*dgamma_r(r, z, p) - e_2gamma*df_r(r, z, p))/(f_val*f_val); 
}
double dg_rr_z(double r, double z, Params *p) { 
    double f_val = f(r, z, p);
    double e_2gamma = exp(2*my_gamma(r, z, p));
    return (f_val*e_2gamma*2*dgamma_z(r, z, p) - e_2gamma*df_z(r, z, p))/(f_val*f_val); 
}

double Det(double r, double z, Params *p) {return g_tt(r, z, p)*g_ff(r, z, p) - powl(g_tf(r, z, p), 2); }

double** make_g(){
    double** g = malloc(4 * sizeof(double*));
    for (int i = 0; i < 4; i++) {
        g[i] = calloc(4, sizeof(double));
    }
    return g;
}

void update_g(double r, double z, double** g, Params *p) {
    g[0][0] = g_tt(r, z, p);
    g[0][1] = g_tf(r, z, p);

    g[1][0] = g_tf(r, z, p);
    g[1][1] = g_ff(r, z, p);

    g[2][2] = g_rr(r, z, p);

    g[3][3] = g_rr(r, z, p);
}

void free_g(double** g) {
    for (int i = 0; i < 4; i++) {
        free(g[i]);
    }
    free(g);
}

double** make_g_inv() {
    double** g_inv = malloc(4 * sizeof(double*));
    for (int i = 0; i < 4; i++) {
        g_inv[i] = calloc(4, sizeof(double));
    }

    return g_inv;
}

void update_g_inv(double r, double z, double** g_inv, Params *p) {
    double Det_val = Det(r, z, p);

    g_inv[0][0] = g_ff(r, z, p)/Det_val;
    g_inv[0][1] = - g_tf(r, z, p)/Det_val;

    g_inv[1][0] = - g_tf(r, z, p)/Det_val;
    g_inv[1][1] = g_tt(r, z, p)/Det_val;

    g_inv[2][2] = 1/g_rr(r, z, p);

    g_inv[3][3] = 1/g_rr(r, z, p);
}

void free_g_inv(double** g_inv) {
    for (int i = 0; i < 4; i++) {
        free(g_inv[i]);
    }
    free(g_inv);
}

double*** make_dg() {
    double*** dg = malloc(4 * sizeof(double**));
    for (int i = 0; i < 4; ++i) {
        dg[i] = malloc(4 * sizeof(double*));
        for (int j = 0; j < 4; ++j) {
            dg[i][j] = calloc(4, sizeof(double));
        }
    }

    return dg;
}

void update_dg(double r, double z, double*** dg, Params *p) {
    dg[2][0][0] = dg_tt_r(r, z, p);
    dg[2][0][1] = dg_tf_r(r, z, p);
    dg[2][1][0] = dg_tf_r(r, z, p);
    dg[2][1][1] = dg_ff_r(r, z, p);
    dg[2][2][2] = dg_rr_r(r, z, p);
    dg[2][3][3] = dg_rr_r(r, z, p);
    
    dg[3][0][0] = dg_tt_z(r, z, p);
    dg[3][0][1] = dg_tf_z(r, z, p);
    dg[3][1][0] = dg_tf_z(r, z, p);
    dg[3][1][1] = dg_ff_z(r, z, p);
    dg[3][2][2] = dg_rr_z(r, z, p);
    dg[3][3][3] = dg_rr_z(r, z, p);
}

void free_dg(double*** dg) {
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            free(dg[i][j]);
        }
        free(dg[i]);
    }
    free(dg);
}

double*** generate_Christoffel_symbols(double r, double z, Params *p, double** g, double** g_inv, double*** dg) {

    double*** Christoffel = malloc(4 * sizeof(double**));
    for (int i = 0; i < 4; ++i) {
        Christoffel[i] = malloc(4 * sizeof(double*));
        for (int j = 0; j < 4; ++j) {
            Christoffel[i][j] = calloc(4, sizeof(double));
        }
    }

    for (int mu = 0; mu < 4 ; mu++){
        for (int kappa = 0; kappa < 4 ; kappa++){
            for (int lambda = 0; lambda < 4 ; lambda++){
                for (int sigma = 0; sigma < 4 ; sigma++){
                    Christoffel[mu][kappa][lambda] += 0.5*g_inv[mu][sigma]*(dg[lambda][sigma][kappa] + dg[kappa][lambda][sigma] - dg[sigma][kappa][lambda]);
                }
            }
        }
    }
    return Christoffel;
}

void free_Christoffel(double*** Christoffel) {
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            free(Christoffel[i][j]);
        }
        free(Christoffel[i]);
    }
    free(Christoffel);
}

/* 
const double L_z = 3;
const double E = 0.95;
Params p = {.M = 1.0, .J = 0.3, .M2 = -0.1, .S3 = 0.05, .M4 = 0.01};

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

int main() {
    double* state_vector = (double*)calloc(8, sizeof(double));
    state_vector[2] = 7;
    state_vector[3] = 0.2;
    state_vector = initialize_velocity(state_vector);
    double** g = make_g();
    double** g_inv = make_g_inv();
    double*** dg = make_dg();
    update_g(state_vector[2], state_vector[3], g, &p);
    update_g_inv(state_vector[2], state_vector[3], g_inv, &p);
    update_dg(state_vector[2], state_vector[3], dg, &p);

    for (int i=0;i<8;i++){
        printf("%Le ", state_vector[i]);
    }
    printf("\n");
    double*** Christoffel = generate_Christoffel_symbols(7, 0.2, &p, g, g_inv, dg);
    for (int i=0;i<4;i++) {
        for (int j=0;j<4;j++) {
            for (int k=0;k<4;k++) {
                printf("%.16Lf ", Christoffel[i][j][k]);
            }
            printf("\n");
        }
        printf("\n");   
    }
    // for (int i=0;i<4;i++) {
    //     for (int j=0;j<4;j++) {
    //         for (int k=0;k<4;k++) {
    //             printf("%.16f, ", dg[i][j][k]);
    //         }
    //         printf("\n");
    //     }
    //     printf("\n");
    // }
    // printf("%.16f\n", f(state_vector[2], state_vector[3], &p));
}
 */