#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct {
    double M, J, M2, S3, M4;
} Params;

Params p = {1.0, 0.37, -0.1, 0.09, -0.17};

double pythagorean(double r, double z) {
    return r*r + z*z;
}

double A(double r, double z) {
    double term1 = 8*r*r*z*z*(24*p.J*p.J*p.M + 17*p.M*p.M*p.M2 + 21*p.M4);
    double term2 = pow(r, 4)*(-10*p.J*p.J*p.M + 7*pow(p.M, 5) + 32*p.M2*p.M*p.M - 21*p.M4);
    double term3 = 8*pow(z, 4)*(20*p.J*p.J*p.M - 7*pow(p.M, 5) - 22*p.M2*p.M*p.M - 7*p.M4);
    return term1 + term2 + term3;
}

double dA_r(double r, double z) {
    double term1 = 16*r*z*z*(24*p.J*p.J*p.M + 17*p.M*p.M*p.M2 + 21*p.M4);
    double term2 = 4*pow(r, 3)*(-10*p.J*p.J*p.M + 7*pow(p.M, 5) + 32*p.M2*p.M*p.M - 21*p.M4);
    double term3 = 0;
    return term1 + term2 + term3;
}

double dA_z(double r, double z) {
    double term1 = 16*r*r*z*(24*p.J*p.J*p.M + 17*p.M*p.M*p.M2 + 21*p.M4);
    double term2 = 0;
    double term3 = 32*pow(z, 3)*(20*p.J*p.J*p.M - 7*pow(p.M, 5) - 22*p.M2*p.M*p.M - 7*p.M4);
    return term1 + term2 + term3;
}

double B(double r, double z) {
    double term1 = pow(r, 4)*(10*p.J*p.J + 10*p.M2*pow(p.M, 3) + 21*p.M4*p.M + 7*p.M2*p.M2);
    double term2 = 4*pow(z, 4)*(-40*p.J*p.J*p.M*p.M - 14*p.J*p.S3 + 7*pow(p.M, 6) + 30*p.M2*pow(p.M, 3) + 14*p.M4*p.M + 7*p.M2*p.M2);
    double term3 = - 4*r*r*z*z*(27*p.J*p.J*p.M*p.M - 21*p.J*p.S3 + 7*pow(p.M, 6) + 48*p.M2*pow(p.M, 3) + 42*p.M4*p.M + 7*p.M2*p.M2);
    return term1 + term2 + term3;
}

double dB_r(double r, double z) {
    double term1 = 4*pow(r, 3)*(10*p.J*p.J + 10*p.M2*pow(p.M, 3) + 21*p.M4*p.M + 7*p.M2*p.M2);
    double term2 = 0;
    double term3 = - 8*r*z*z*(27*p.J*p.J*p.M*p.M - 21*p.J*p.S3 + 7*pow(p.M, 6) + 48*p.M2*pow(p.M, 3) + 42*p.M4*p.M + 7*p.M2*p.M2);
    return term1 + term2 + term3;
}

double dB_z(double r, double z) {
    double term1 = 0;
    double term2 = 16*pow(z, 3)*(-40*p.J*p.J*p.M*p.M - 14*p.J*p.S3 + 7*pow(p.M, 6) + 30*p.M2*pow(p.M, 3) + 14*p.M4*p.M + 7*p.M2*p.M2);
    double term3 = - 8*r*r*z*(27*p.J*p.J*p.M*p.M - 21*p.J*p.S3 + 7*pow(p.M, 6) + 48*p.M2*pow(p.M, 3) + 42*p.M4*p.M + 7*p.M2*p.M2);
    return term1 + term2 + term3;
}

double H(double r, double z) {
    double term1 = 4*r*r*z*z*(p.J*(p.M2 - 2*pow(p.M, 3)) - 3*p.M*p.S3);
    double term2 = pow(r, 4)*(p.J*p.M2 + 3*p.M*p.S3);
    return term1 + term2;
}

double dH_r(double r, double z) {
    double term1 = 8*r*z*z*(p.J*(p.M2 - 2*pow(p.M, 3)) - 3*p.M*p.S3);
    double term2 = 4*pow(r, 3)*(p.J*p.M2 + 3*p.M*p.S3);
    return term1 + term2;
}

double dH_z(double r, double z) {
    double term1 = 8*r*r*z*(p.J*(p.M2 - 2*pow(p.M, 3)) - 3*p.M*p.S3);
    double term2 = 0;
    return term1 + term2;
}

double G(double r, double z) {
    double term1 = pow(p.J, 3)*(- pow(r, 4)
                                - 8*pow(z, 4)
                                + 12*r*r*z*z);

    double term2 = p.J*p.M*((+ pow(p.M, 3) + 2*p.M2)*pow(r, 4) 
                             - 8*(2*pow(p.M, 3) + 2*p.M2)*pow(z, 4)
                             + 4*(pow(p.M, 3) + 10*p.M2)*r*r*z*z);

    double term3 = p.M*p.M*p.S3*(+ 3*pow(r, 4)
                                 - 40*pow(z, 4)
                                 + 12*r*r*z*z);
    return r*r*(term1 + term2 + term3);
}

double dG_r(double r, double z) {
    double term1 = pow(p.J, 3)*(- pow(r, 4)
                                - 8*pow(z, 4)
                                + 12*r*r*z*z);

    double term2 = p.J*p.M*(+ (pow(p.M, 3) + 2*p.M2)*pow(r, 4) 
                            - 8*(2*pow(p.M, 3) + 2*p.M2)*pow(z, 4)
                            + 4*(pow(p.M, 3) + 10*p.M2)*r*r*z*z);

    double term3 = p.M*p.M*p.S3*(+ 3*pow(r, 4)
                                 - 40*pow(z, 4)
                                 + 12*r*r*z*z);

    double left = r*r;
    double dright = + pow(p.J, 3)*(- 4*pow(r, 3)
                               - 0
                               + 24*r*z*z)
                    + p.J*p.M*(+ (pow(p.M, 3) + 2*p.M2)*4*pow(r, 3)
                               - 0
                               + 4*(pow(p.M, 3) + 10*p.M2)*2*r*z*z)
                    + p.M*p.M*p.S3*(+ 12*pow(r, 3)
                                    - 0
                                    + 24*r*z*z);
    double right = term1 + term2 + term3;
    double dleft = 2*r;
    return left*dright + right*dleft;
}

double dG_z(double r, double z) {
    double term1 = pow(p.J, 3)*(- 0
                                - 32*pow(z, 3)
                                + 24*r*r*z);

    double term2 = p.J*p.M*(+ 0 
                            - 32*(2*pow(p.M, 3) + 2*p.M2)*pow(z, 3)
                            + 8*(pow(p.M, 3) + 10*p.M2)*r*r*z);

    double term3 = p.M*p.M*p.S3*(+ 0
                                 - 160*pow(z, 3)
                                 + 24*r*r*z);
    
    return r*r*(term1 + term2 + term3);
}

double F(double r, double z) {
    return (+ pow(r, 4)*(p.S3 - p.J*p.M*p.M)
            - 4*r*r*z*z*(p.J*p.M*p.M + p.S3));
}

double dF_r(double r, double z) {
    return (+ 4*pow(r, 3)*(p.S3 - p.J*p.M*p.M)
            - 8*r*z*z*(p.J*p.M*p.M + p.S3)); 
}

double dF_z(double r, double z) {
    return (+ 0
            - 8*r*r*z*(p.J*p.M*p.M + p.S3)); 
}

double f(double r, double z) {
    double pyth = pythagorean(r, z);
    double term1 = 1;
    double term2 = - (2*p.M)/(sqrt(pyth));
    double term3 = (2*p.M*p.M)/(pyth);
    double term4 = ((p.M2 - pow(p.M, 3))*r*r - 2*(pow(p.M, 3) + p.M2)*z*z)/(pow(pyth, 5.0/2.0));
    double term5 = (2*z*z*(-p.J*p.J + pow(p.M, 4) + 2*p.M2*p.M) - 2*p.M*p.M2*r*r)/(pow(pyth, 3));
    double term6 = (A(r, z))/(28*pow(pyth, 9.0/2.0));
    double term7 = (B(r, z))/(14*pow(pyth, 5));
    return term1 + term2 + term3 + term4 + term5 + term6 + term7;
}

double df_r(double r, double z) {
    double pyth = pythagorean(r, z);
    double term1 = 0;
    double term2 = (2*r*p.M)/(pow(pyth, 3.0/2.0));
    double term3 = - (4*r*p.M*p.M)/(pyth*pyth);
    double term4 = (2*r*pow(pyth, 5.0/2.0)*(p.M2 - pow(p.M, 3)) - 5*r*pow(pyth, 3.0/2.0)*((p.M2 - pow(p.M, 3))*r*r - 2*(pow(p.M, 3) + p.M2)*z*z))/(pow(pyth, 5));
    double term5 = (-4*r*pow(pyth, 3)*p.M*p.M2 - 6*r*pow(pyth, 2)*(2*z*z*(-p.J*p.J + pow(p.M, 4) + 2*p.M2*p.M) - 2*p.M*p.M2*r*r))/(pow(pyth, 6));
    double term6 = (pow(pyth, 9.0/2.0)*dA_r(r, z) - 9*r*pow(pyth, 7.0/2.0)*A(r, z))/(28*pow(pyth, 9));
    double term7 = (pow(pyth, 5)*dB_r(r, z) - 10*r*pow(pyth, 4)*B(r, z))/(14*pow(pyth, 10));
    return term1 + term2 + term3 + term4 + term5 + term6 + term7;
}

double df_z(double r, double z) {
    double pyth = pythagorean(r, z);
    double term1 = 0;
    double term2 = (2*z*p.M)/(pow(pyth, 3.0/2.0));
    double term3 = - (4*z*p.M*p.M)/(pyth*pyth);
    double term4 = (-4*z*pow(pyth, 5.0/2.0)*(pow(p.M, 3) + p.M2) - 5*z*pow(pyth, 3.0/2.0)*((p.M2 - pow(p.M, 3))*r*r - 2*(pow(p.M, 3) + p.M2)*z*z))/(pow(pyth, 5));
    double term5 = (4*z*pow(pyth, 3)*(-p.J*p.J + pow(p.M, 4) + 2*p.M2*p.M) - 6*pow(pyth, 2)*z*(2*z*z*(-p.J*p.J + pow(p.M, 4) + 2*p.M2*p.M) - 2*p.M*p.M2*r*r))/(pow(pyth, 6));
    double term6 = (pow(pyth, 9.0/2.0)*dA_z(r, z) - 9*z*pow(pyth, 7.0/2.0)*A(r, z))/(28*pow(pyth, 9));
    double term7 = (pow(pyth, 5)*dB_z(r, z) - 10*z*pow(pyth, 4)*B(r, z))/(14*pow(pyth, 10));
    return term1 + term2 + term3 + term4 + term5 + term6 + term7;
}

double omega(double r, double z) {
    double pyth = pythagorean(r, z);
    double term1 = - (2*p.J*r*r)/(pow(pyth, 3.0/2.0));
    double term2 = - (2*p.J*p.M*r*r)/(pyth*pyth);
    double term3 = + (F(r, z))/(pow(pyth, 7.0/2.0));
    double term4 = + (H(r, z))/(2*pow(pyth, 4));
    double term5 = + (G(r, z))/(4*pow(pyth, 11.0/2.0));
    return term1 + term2 + term3 + term4 + term5;
}

double domega_r(double r, double z) {
    double pyth = pythagorean(r, z);
    double term1 = - (4*r*pow(pyth, 3.0/2.0)*p.J - 6*pow(r, 3)*pow(pyth, 1.0/2.0)*p.J)/(pow(pyth, 3));
    double term2 = - (4*r*pyth*pyth*p.J*p.M - 8*pow(r, 3)*pyth*p.J*p.M)/(pow(pyth, 4));
    double term3 = + (dF_r(r, z)*pow(pyth, 7.0/2.0) - 7*r*pow(pyth, 5.0/2.0)*F(r, z))/(pow(pyth, 7));
    double term4 = + (2*dH_r(r, z)*pow(pyth, 4) - 16*r*pow(pyth, 3)*H(r, z))/(4*pow(pyth, 8));
    double term5 = + (4*dG_r(r, z)*pow(pyth, 11.0/2.0) - 44*r*pow(pyth, 9.0/2.0)*G(r, z))/(16*pow(pyth, 11));
    return term1 + term2 + term3 + term4 + term5;
}

double domega_z(double r, double z) {
    double pyth = pythagorean(r, z);
    double term1 = + (6*r*r*z*pow(pyth, 1.0/2.0)*p.J)/(pow(pyth, 3));
    double term2 = + (8*r*r*z*pyth*p.J*p.M)/(pow(pyth, 4));
    double term3 = + (dF_z(r, z)*pow(pyth, 7.0/2.0) - 7*z*pow(pyth, 5.0/2.0)*F(r, z))/(pow(pyth, 7));
    double term4 = + (2*dH_z(r, z)*pow(pyth, 4) - 16*z*pow(pyth, 3)*H(r, z))/(4*pow(pyth, 8));
    double term5 = + (4*dG_z(r, z)*pow(pyth, 11.0/2.0) - 44*z*pow(pyth, 9.0/2.0)*G(r, z))/(16*pow(pyth, 11));
    return term1 + term2 + term3 + term4 + term5;
}

double my_gamma(double r, double z) {
    double pyth = pythagorean(r, z);
    double term1 = + (r*r*(p.J*p.J*(r*r - 8*z*z) + p.M*(pow(p.M, 3) + 3*p.M2)*(r*r - 4*z*z)))/(4*pow(pyth, 4));
    double term2 = - (p.M*p.M*r*r)/(2*pyth*pyth);
    return term1 + term2;
}

double dgamma_r(double r, double z) {
    double pyth = pythagorean(r, z);
    double left = r*r;
    double dright = p.J*p.J*2*r + p.M*(pow(p.M, 3) + 3*p.M2)*2*r;
    double right = p.J*p.J*(r*r - 8*z*z) + p.M*(pow(p.M, 3) + 3*p.M2)*(r*r - 4*z*z);
    double dleft = 2*r;
    double denom = 4*pow(pyth, 4);
    double term1 = + ((left*dright + right*dleft)*denom - left*right*32*r*pow(pyth, 3))/(denom*denom);
    double term2 = - (4*r*pyth*pyth*p.M*p.M - 8*pow(r, 3)*pyth*p.M*p.M)/(denom);
    return term1 + term2;
}

double dgamma_z(double r, double z) {
    double pyth = pythagorean(r, z);
    double term1 = + (4*r*r*pow(pyth, 4)*(-p.J*p.J*16*z - p.M*(pow(p.M, 3) + 3*p.M2)*8*z) - 32*r*r*z*pow(pyth, 3)*(p.J*p.J*(r*r - 8*z*z) + p.M*(pow(p.M, 3) + 3*p.M2)*(r*r - 4*z*z)))/(16*pow(pyth, 8));
    double term2 = + (8*r*r*z*pyth)/(4*pow(pyth, 4));
    return term1 + term2;
}

double g_tt(double r, double z) { return - f(r, z); }
double dg_tt_r(double r, double z) { return - df_r(r, z); }
double dg_tt_z(double r, double z) { return - df_z(r, z); }

double g_tf(double r, double z) { return omega(r, z)*f(r, z); }
double dg_tf_r(double r, double z) { return omega(r, z)*df_r(r, z) + domega_r(r, z)*f(r, z); }
double dg_tf_z(double r, double z) { return omega(r, z)*df_z(r, z) + domega_z(r, z)*f(r, z); }

double g_ff(double r, double z) { return - f(r, z)*omega(r, z)*omega(r, z) + (r*r)/f(r, z); }
double dg_ff_r(double r, double z) {
    double f_val = f(r,  z);
    double omega_val = omega(r,  z);
    return - 2*f_val*omega_val*domega_r(r, z) - df_r(r, z)*omega_val*omega_val + (2*r*f_val - r*r*df_r(r, z))/(f_val*f_val); 
}
double dg_ff_z(double r, double z) {
    double f_val = f(r,  z);
    double omega_val = omega(r,  z);
    return - 2*f_val*omega_val*domega_z(r, z) - df_z(r, z)*omega_val*omega_val - (r*r*df_z(r, z))/(f_val*f_val);
}

double g_rr(double r, double z) { return exp(2*my_gamma(r, z))/f(r, z); }
double dg_rr_r(double r, double z) { 
    double f_val = f(r, z);
    double e_2gamma = exp(2*my_gamma(r, z));
    return (f_val*e_2gamma*2*dgamma_r(r, z) - e_2gamma*df_r(r, z))/(f_val*f_val); 
}
double dg_rr_z(double r, double z) { 
    double f_val = f(r, z);
    double e_2gamma = exp(2*my_gamma(r, z));
    return (f_val*e_2gamma*2*dgamma_z(r, z) - e_2gamma*df_z(r, z))/(f_val*f_val); 
}

double Det(double r, double z) {return g_tt(r, z)*g_ff(r, z) - g_tf(r, z)*g_tf(r, z); }

double** make_g(double r, double z){
    double** g = malloc(4 * sizeof(double*));
    for (int i = 0; i < 4; i++) {
        g[i] = malloc(4* sizeof(double));
    }

    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            g[i][j] = 0;
        }
    }

    g[0][0] = g_tt(r, z);
    g[0][1] = g_tf(r, z);

    g[1][0] = g_tf(r, z);
    g[1][1] = g_ff(r, z);

    g[2][2] = g_rr(r, z);

    g[3][3] = g_rr(r, z);

    return g;
}

void free_g(double** g) {
    for (int i = 0; i < 4; i++) {
        free(g[i]);
    }
    free(g);
}

double** make_g_inv(double r, double z) {
    double** g_inv = malloc(4 * sizeof(double*));
    for (int i = 0; i < 4; i++) {
        g_inv[i] = malloc(4* sizeof(double));
    }

    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            g_inv[i][j] = 0;
        }
    }

    double Det_val = Det(r, z);

    g_inv[0][0] = g_ff(r, z)/Det_val;
    g_inv[0][1] = - g_tf(r, z)/Det_val;

    g_inv[1][0] = - g_tf(r, z)/Det_val;
    g_inv[1][1] = g_tt(r, z)/Det_val;

    g_inv[2][2] = 1/g_rr(r, z);

    g_inv[3][3] = 1/g_rr(r, z);

    return g_inv;
}

void free_g_inv(double** g_inv) {
    for (int i = 0; i < 4; i++) {
        free(g_inv[i]);
    }
    free(g_inv);
}

double*** make_dg(double r, double z) {
    double*** dg = malloc(4 * sizeof(double**));
    for (int i = 0; i < 4; ++i) {
        dg[i] = malloc(4 * sizeof(double*));
        for (int j = 0; j < 4; ++j) {
            dg[i][j] = malloc(4 * sizeof(double));
        }
    }

    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            for (int k = 0; k < 4; k++) {
                dg[i][j][k] = 0;
            }
        }
    }

    dg[2][0][0] = dg_tt_r(r, z);
    dg[2][0][1] = dg_tf_r(r, z);
    dg[2][1][0] = dg_tf_r(r, z);
    dg[2][2][2] = dg_rr_r(r, z);
    dg[2][3][3] = dg_rr_r(r, z);
    
    dg[3][0][0] = dg_tt_z(r, z);
    dg[3][0][1] = dg_tf_z(r, z);
    dg[3][1][0] = dg_tf_z(r, z);
    dg[3][2][2] = dg_rr_z(r, z);
    dg[3][3][3] = dg_rr_z(r, z);

    return dg;
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

double*** generate_Christoffel_symbols(double r, double z) {
    double** g_val = make_g(r, z);
    double*** dg_val = make_dg(r, z);
    double** g_inv_val = make_g_inv(r, z);
    double*** Christoffel = malloc(4 * sizeof(double**));
    for (int i = 0; i < 4; ++i) {
        Christoffel[i] = malloc(4 * sizeof(double*));
        for (int j = 0; j < 4; ++j) {
            Christoffel[i][j] = malloc(4 * sizeof(double));
        }
    }

    for (int mu = 0; mu < 4 ; mu++){
        for (int kappa = 0; kappa < 4 ; kappa++){
            for (int lambda = 0; lambda < 4 ; lambda++){
                Christoffel[mu][kappa][lambda] = 0;
                for (int sigma = 0; sigma < 4 ; sigma++){
                    Christoffel[mu][kappa][lambda] += 0.5*g_inv_val[mu][sigma]*(dg_val[lambda][sigma][kappa] + dg_val[kappa][lambda][sigma] - dg_val[sigma][kappa][lambda]);
                }
            }
        }
    }

    free_g(g_val);
    free_g_inv(g_inv_val);
    free_dg(dg_val);

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
