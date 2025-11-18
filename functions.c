#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
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
    return term1 + term2;
}

double dA_z(double r, double z, Params *p) {
    double term1 = 16*r*r*z*(24*(*p).J*(*p).J*(*p).M + 17*(*p).M*(*p).M*(*p).M2 + 21*(*p).M4);
    double term2 = 32*powl(z, 3)*(20*(*p).J*(*p).J*(*p).M - 7*powl((*p).M, 5) - 22*(*p).M2*(*p).M*(*p).M - 7*(*p).M4);
    return term1 + term2;
}

double ddA_r_r(double r, double z, Params *p) {
    double term1 = 16*z*z*(24*(*p).J*(*p).J*(*p).M + 17*(*p).M*(*p).M*(*p).M2 + 21*(*p).M4);
    double term2 = 12*r*r*(-10*(*p).J*(*p).J*(*p).M + 7*powl((*p).M, 5) + 32*(*p).M2*(*p).M*(*p).M - 21*(*p).M4);
    return term1 + term2;
}

double ddA_r_z(double r, double z, Params *p) {
    return 32*r*z*(24*(*p).J*(*p).J*(*p).M + 17*(*p).M*(*p).M*(*p).M2 + 21*(*p).M4);
}

double ddA_z_z(double r, double z, Params *p) {
    double term1 = 16*r*r*(24*(*p).J*(*p).J*(*p).M + 17*(*p).M*(*p).M*(*p).M2 + 21*(*p).M4);
    double term2 = 96*z*z*(20*(*p).J*(*p).J*(*p).M - 7*powl((*p).M, 5) - 22*(*p).M2*(*p).M*(*p).M - 7*(*p).M4);
    return term1 + term2;
}


double B(double r, double z, Params *p) {
    double term1 = powl(r, 4)*(10*(*p).J*(*p).J + 10*(*p).M2*powl((*p).M, 3) + 21*(*p).M4*(*p).M + 7*(*p).M2*(*p).M2);
    double term2 = 4*powl(z, 4)*(-40*(*p).J*(*p).J*(*p).M*(*p).M - 14*(*p).J*(*p).S3 + 7*powl((*p).M, 6) + 30*(*p).M2*powl((*p).M, 3) + 14*(*p).M4*(*p).M + 7*(*p).M2*(*p).M2);
    double term3 = - 4*r*r*z*z*(27*(*p).J*(*p).J*(*p).M*(*p).M - 21*(*p).J*(*p).S3 + 7*powl((*p).M, 6) + 48*(*p).M2*powl((*p).M, 3) + 42*(*p).M4*(*p).M + 7*(*p).M2*(*p).M2);
    return term1 + term2 + term3;
}

double dB_r(double r, double z, Params *p) {
    double term1 = 4*powl(r, 3)*(10*(*p).J*(*p).J*(*p).M*(*p).M + 10*(*p).M2*powl((*p).M, 3) + 21*(*p).M4*(*p).M + 7*(*p).M2*(*p).M2);
    double term2 = - 8*r*z*z*(27*(*p).J*(*p).J*(*p).M*(*p).M - 21*(*p).J*(*p).S3 + 7*powl((*p).M, 6) + 48*(*p).M2*powl((*p).M, 3) + 42*(*p).M4*(*p).M + 7*(*p).M2*(*p).M2);
    return term1 + term2;
}

double dB_z(double r, double z, Params *p) {
    double term1 = 16*powl(z, 3)*(-40*(*p).J*(*p).J*(*p).M*(*p).M - 14*(*p).J*(*p).S3 + 7*powl((*p).M, 6) + 30*(*p).M2*powl((*p).M, 3) + 14*(*p).M4*(*p).M + 7*(*p).M2*(*p).M2);
    double term2 = - 8*r*r*z*(27*(*p).J*(*p).J*(*p).M*(*p).M - 21*(*p).J*(*p).S3 + 7*powl((*p).M, 6) + 48*(*p).M2*powl((*p).M, 3) + 42*(*p).M4*(*p).M + 7*(*p).M2*(*p).M2);
    return term1 + term2;
}

double ddB_r_r(double r, double z, Params *p) {
    double term1 = 12*r*r*(10*(*p).J*(*p).J*(*p).M*(*p).M + 10*(*p).M2*powl((*p).M, 3) + 21*(*p).M4*(*p).M + 7*(*p).M2*(*p).M2);
    double term2 = - 8*z*z*(27*(*p).J*(*p).J*(*p).M*(*p).M - 21*(*p).J*(*p).S3 + 7*powl((*p).M, 6) + 48*(*p).M2*powl((*p).M, 3) + 42*(*p).M4*(*p).M + 7*(*p).M2*(*p).M2);
    return term1 + term2;
}

double ddB_r_z(double r, double z, Params *p) {
    return - 16*r*z*(27*(*p).J*(*p).J*(*p).M*(*p).M - 21*(*p).J*(*p).S3 + 7*powl((*p).M, 6) + 48*(*p).M2*powl((*p).M, 3) + 42*(*p).M4*(*p).M + 7*(*p).M2*(*p).M2);
}

double ddB_z_z(double r, double z, Params *p) {
    double term1 = 48*z*z*(-40*(*p).J*(*p).J*(*p).M*(*p).M - 14*(*p).J*(*p).S3 + 7*powl((*p).M, 6) + 30*(*p).M2*powl((*p).M, 3) + 14*(*p).M4*(*p).M + 7*(*p).M2*(*p).M2);
    double term2 = - 8*r*r*(27*(*p).J*(*p).J*(*p).M*(*p).M - 21*(*p).J*(*p).S3 + 7*powl((*p).M, 6) + 48*(*p).M2*powl((*p).M, 3) + 42*(*p).M4*(*p).M + 7*(*p).M2*(*p).M2);
    return term1 + term2;
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
    return 8*r*r*z*((*p).J*((*p).M2 - 2*powl((*p).M, 3)) - 3*(*p).M*(*p).S3);
}

double ddH_r_r(double r, double z, Params *p) {
    double term1 = 8*z*z*((*p).J*((*p).M2 - 2*powl((*p).M, 3)) - 3*(*p).M*(*p).S3);
    double term2 = 12*r*r*((*p).J*(*p).M2 + 3*(*p).M*(*p).S3);
    return term1 + term2;
}

double ddH_r_z(double r, double z, Params *p) {
    return 16*r*z*((*p).J*((*p).M2 - 2*powl((*p).M, 3)) - 3*(*p).M*(*p).S3);
}

double ddH_z_z(double r, double z, Params *p) {
    return 8*r*r*((*p).J*((*p).M2 - 2*powl((*p).M, 3)) - 3*(*p).M*(*p).S3);
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
    double term1 = powl((*p).J, 3)*(- 4*z*z + 3*r*r);
    double term2 = (*p).J*(*p).M*(- 4*z*z*(3*powl((*p).M, 3) + 2*(*p).M2) + (powl((*p).M, 3) + 10*(*p).M2)*r*r);
    double term3 = (*p).M*(*p).M*(*p).S3*(- 20*z*z + 3*r*r);
    
    return 8*z*r*r*(term1 + term2 + term3);
}

double ddG_r_r(double r, double z, Params *p) {
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
    double dright = 4*r*(powl((*p).J, 3)*(- r*r + 8*z*z)
                        + (*p).J*(*p).M*(r*r*(powl((*p).M, 3) + 2*(*p).M2) + 2*z*z*(powl((*p).M, 3) + 10*(*p).M2))
                        + 3*(*p).M*(*p).M*(*p).S3*(r*r + 2*z*z));
    double ddright = 4*(*p).J*powl((*p).M, 4)*(3*r*r + 2*z*z) + 8*(*p).J*(*p).M*(*p).M2*(3*r*r + 10*z*z) + 4*powl((*p).J, 3)*(-3*r*r + 8*z*z) + 6*r*(*p).M*(*p).M*(*p).S3;
    double right = term1 + term2 + term3;
    double dleft = 2*r;
    double ddleft = 2;
    return left*ddright + 2*dleft*dright + right*ddleft;
}

double ddG_r_z(double r, double z, Params *p) {
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
    double dright_z = 8*z*(*p).J*(*p).J*(-4*z*z + 3*r*r) + 8*z*(*p).J*(*p).M*(- 4*z*z*(3*powl((*p).M, 3) + 2*(*p).M2) + r*r*(powl((*p).M, 3) + 10*(*p).M2)) + 8*z*(*p).M*(*p).M*(-20*z*z + 3*r*r);
    double ddright = 4*z*(16*r*powl((*p).J, 3) + 4*r*(*p).J*(*p).M*(powl((*p).M, 3) + 10*(*p).M2) + 3*(*p).M*(*p).M*(*p).S3);
    double right = term1 + term2 + term3;
    double dleft_r = 2*r;
    return left*ddright + dleft_r*dright_z;
}

double ddG_z_z(double r, double z, Params *p) {
    double term1 = - 8*z*powl((*p).J, 3);
    double term2 = - 8*z*(*p).J*(*p).M*(3*powl((*p).M, 3) + 2*(*p).M2);
    double term3 = - 40*z*(*p).M*(*p).M*(*p).S3;

    double dright = term1 + term2 + term3;
    
    term1 = powl((*p).J, 3)*(- 4*z*z + 3*r*r);
    term2 = (*p).J*(*p).M*(- 4*z*z*(3*powl((*p).M, 3) + 2*(*p).M2) + (powl((*p).M, 3) + 10*(*p).M2)*r*r);
    term3 = (*p).M*(*p).M*(*p).S3*(- 20*z*z + 3*r*r);

    double right = term1 + term2 + term3;

    return 8*r*r*(z*dright + right);
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
    return (- 8*r*r*z*((*p).J*(*p).M*(*p).M + (*p).S3)); 
}

double ddF_r_r(double r, double z, Params *p) {
    return (+ 12*r*r*((*p).S3 - (*p).J*(*p).M*(*p).M)
            - 8*z*z*((*p).J*(*p).M*(*p).M + (*p).S3)); 
}

double ddF_r_z(double r, double z, Params *p) {
    return (- 16*r*z*((*p).J*(*p).M*(*p).M + (*p).S3));
}

double ddF_z_z(double r, double z, Params *p) {
    return - 8*r*r*((*p).J*(*p).M*(*p).M + (*p).S3); 
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

double ddf_r_r(double r, double z, Params *p) {
    double pyth = pythagorean(r, z, p);
    double term1 = 0;
    double term2 = 2*(*p).M*(pyth - 3*r*r)/powl(pyth, 5.0/2.0);
    double term3 = - 4*(*p).M*(*p).M*(pyth - 4*r*r)/powl(pyth, 3);
    double term4 = (((*p).M2 - powl((*p).M, 3))*(2*pyth*pyth - 17*r*r*pyth + 15*powl(r, 4)) + 2*z*z*((*p).M2 - powl((*p).M, 3))*(pyth - 3*r*r))/powl(pyth, 5.0/2.0);
    double term5 = 4*((*p).M*(*p).M2*(-pyth*pyth + 19*r*r*pyth-3*powl(r, 4)) + 3*z*z*((*p).J*(*p).J + powl((*p).M, 4) + 2*(*p).M*(*p).M2)*(-pyth + 12*r*r))/powl(pyth, 7);
    double term6 = (ddA_r_r(r, z, p)*pyth*pyth - 18*dA_r(r, z, p)*r*pyth - 9*A(r, z, p)*(pyth - 11*r*r))/powl(pyth, 13.0/2.0);
    double term7 = (ddB_r_r(r, z, p)*pyth*pyth - 20*dB_(r, z, p)*r*pyth - 10*B(r, z, p)*(pyth - 12*r*r))/powl(pyth, 7);
    return term1 + term2 + term3 + term4 + term5 + term6 + term7;
}

double ddf_r_z(double r, double z, Params *p) {
    double pyth = pythagorean(r, z, p);
    double term1 = 0;
    double term2 = - 6*(*p).M*r*z/powl(pyth, 5.0/2.0);
    double term3 = 16*(*p).M*(*p).M*r*z/powl(pyth, 3);
    double term4 = r*(2*z*pyth*(9*(*p).M2 + 11*powl((*p).M, 3)) + 15*z*((*p).M2 - powl((*p).M, 3))*(r*r - 2*z*z))/powl(pyth, 5.0/2.0);
    double term5 = 24*r*(((*p).M*(*p).M2*z)*(pyth - 8 + 4*r*r) + z*(-(*p).J*(*p).J + powl((*p).M, 4))*(-pyth + 4*z*z))/powl(pyth, 5);
    double term6 = (ddA_r_z(r, z, p)*pyth*pyth - 9*dA_r(r, z, p)*z*pyth - 9*dA_z(r, z, p)*r*pyth + 99*A(r, z, p)*r*z)/(28*powl(pyth, 13.0/2.0));
    double term7 = (ddB_r_z(r, z, p)*pyth*pyth - 10*dB_r(r, z, p)*z*pyth - 10*dB_z(r, z, p)*r*pyth + 120*A(r, z, p)*r*z)/(14*powl(pyth, 7));
    return term1 + term2 + term3 + term4 + term5 + term6 + term7;
}

double ddf_z_z(double r, double z, Params *p) {
    double pyth = pythagorean(r, z, p);
    double term1 = 0;
    double term2 = 2*(*p).M*(pyth - 3*z*z)/powl(pyth, 5.0/2.0);
    double term3 = - 4*(*p).M*(*p).M*(pyth - 4*z*z)/powl(pyth, 3);
    double term4 = (2*(powl((*p).M, 3) + (*p).M2)*(-2*pyth*pyth + 18*z*z*pyth + 15*powl(z, 4)) + 5*r*r*((*p).M2 - powl((*p).M, 3))*(- pyth + 3*z*z))/powl(pyth, 5.0/2.0);
    double term5 = 4*((-15*z*z*pyth + pyth*pyth + 24*powl(z, 4))*(-(*p).J*(*p).J + powl((*p).M, 4) + 2*(*p).M2*(*p).M) + 3*r*r*(*p).M*(*p).M2*(pyth - 8*z*z))/powl(pyth, 5);
    double term6 = (ddA_z_z(r, z, p)*pyth*pyth - 18*dA_z(r, z, p)*z*pyth - 9*A(r, z, p)*(pyth - 11*z*z))/(28*powl(pyth, 13.0/2.0));
    double term7 = (ddB_z_z(r, z, p)*pyth*pyth - 20*dB_z(r, z, p)*z*pyth - 10*B(r, z, p)*(pyth - 12*z*z))/(14*powl(pyth, 7));
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

double ddomega_r_r(double r, double z, Params *p) {
    double pyth = pythagorean(r, z, p);
    double term1 = 2*(*p).J*(15*r*r*pyth - 15*powl(r, 4) - 2*pyth*pyth)/powl(pyth, 7.0/2.0);
    double term2 = 2*(*p).J*(*p).M*(10*r*r*pyth - pyth*pyth - 12*powl(pyth, 4))/powl(pyth, 4);
    double term3 = (ddF_r_r(r, z, p)*pyth*pyth - 14*dF_r(r, z, p)*r*pyth - 7*F(r, z, p)*(pyth - 9*r*r))/powl(pyth, 11.0/2.0);
    double term4 = (ddH_r_r(r, z, p)*pyth*pyth - 16*dH_r(r, z, p)*r*pyth - 8*H(r, z, p)*(pyth - 10*r*r))/(2*powl(pyth, 6));
    double term5 = (ddG_r_r(r, z, p)*pyth*pyth - 321*dG_r(r, z, p)*r*pyth - G(r, z, p)*(11*pyth - 572*r*r))/(16*powl(pyth, 15.0/2.0));
    return term1 + term2 + term3 + term4 + term5;
}

double ddomega_r_z(double r, double z, Params *p){
    double pyth = pythagorean(r, z, p);
    double term1 = 6*(*p).J*r*(2*z*pyth - 5*r*r*z)/powl(pyth, 7.0/2.0);
    double term2 = 16*(*p).J*(*p).M*r*(z*pyth - 3*r*r*z)/powl(pyth, 4);
    double term3 = (ddF_r_z(r, z, p)*pyth*pyth - 7*dF_r(r, z, p)*z*pyth - 7*dF_z(r, z, p)*r*pyth + 63*F(r, z, p)*r*z)/powl(pyth, 11.0/2.0);
    double term4 = (ddH_r_z(r, z, p)*pyth*pyth - 8*dH_r(r, z, p)*z*pyth - 8*dH_z(r, z, p)*r*pyth + 80*H(r, z, p)*r*z)/(2*powl(pyth, 6));
    double term5 = (ddG_r_z(r, z, p)*pyth*pyth - 11*dG_r(r, z, p)*z*pyth - 11*dG_z(r, z, p)*r*pyth + 143*F(r, z, p)*r*z)/(4*powl(pyth, 15.0/2.0));
    return term1 + term2 + term3 + term4 + term5;
}

double ddomega_z_z(double r, double z, Params *p){
    double pyth = pythagorean(r, z, p);
    double term1 = 6*r*r*(*p).J*(pyth - 5*z*z)/powl(pyth, 3.0/2.0);
    double term2 = 8*r*r*(*p).J*(*p).M*(pyth - 6*z*z)/powl(pyth, 4);
    double term3 = (ddF_z_z(r, z, p)*pyth*pyth - 14*dF_z(r, z, p)*z*pyth - 7*F(r, z, p)*(pyth - 9*z*z))/powl(pyth, 11.0/2.0);
    double term4 = (ddH_z_z(r, z, p)*pyth*pyth - 16*dH_z(r, z, p)*z*pyth - 8*H(r, z, p)*(pyth - 10*z*z))/(2*powl(pyth, 6));
    double term5 = (ddG_z_z(r, z, p)*pyth*pyth - 321*dG_z(r, z, p)*z*pyth - 11*G(r, z, p)*(pyth - 52*z*z))/(16*powl(pyth, 15.0/2.0));
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
    double term2 = - (r*pyth*(*p).M*(*p).M - 2*powl(r, 3)*(*p).M*(*p).M)/(powl(pyth, 3));
    return term1 + term2;
}

double dgamma_z(double r, double z, Params *p) {
    double pyth = pythagorean(r, z, p);
    double term1 = + (4*r*r*powl(pyth, 4)*(-(*p).J*(*p).J*16*z - (*p).M*(powl((*p).M, 3) + 3*(*p).M2)*8*z) - 32*r*r*z*powl(pyth, 3)*((*p).J*(*p).J*(r*r - 8*z*z) + (*p).M*(powl((*p).M, 3) + 3*(*p).M2)*(r*r - 4*z*z)))/(16*powl(pyth, 8));
    double term2 = + (2*r*r*z)/powl(pyth, 3);
    return term1 + term2;
}

double ddgamma_r_r(double r, double z, Params *p) {
    double pyth = pythagorean(r, z, p);
    double left = r*r;
    double dright = (*p).J*(*p).J*2*r + (*p).M*(powl((*p).M, 3) + 3*(*p).M2)*2*r;
    double right = (*p).J*(*p).J*(r*r - 8*z*z) + (*p).M*(powl((*p).M, 3) + 3*(*p).M2)*(r*r - 4*z*z);
    double dleft = 2*r;
    double denom = 4*powl(pyth, 4);
    double num = (left*dright + right*dleft)*denom - left*right*32*r*powl(pyth, 3);
    double dnum = (powl(dleft*dright, 2) + left*(2*(*p).J*(*p).J + 2*(*p).M*(powl((*p).M, 3) + 3*(*p).M2)) + 2*right)*denom;
    dnum += (left*dright + right*dleft)*32*r*powl(pyth, 3);
    dnum -= 32*(dleft*right*r*powl(pyth, 3) + left*dright*r*powl(pyth, 3) + left*right*powl(pyth, 3) + left*right*r*r*6*pyth*pyth);
    double term1 = (dnum*denom*denom - num*256*r*powl(pyth, 7))/(256*powl(pyth, 16));
    double term2 = ((4*r*r - pyth)*pyth - 6*r*r*(2*r*r - pyth))/powl(pyth, 4);
    return term1 + term2;
}

double ddgamma_r_z(double r, double z, Params *p) {
    double pyth = pythagorean(r, z, p);
    double num1 = 4*r*r*powl(pyth, 4)*(-(*p).J*(*p).J*16*z - (*p).M*(powl((*p).M, 3) + 3*(*p).M2)*8*z) - 32*r*r*z*powl(pyth, 3)*((*p).J*(*p).J*(r*r - 8*z*z) + (*p).M*(powl((*p).M, 3) + 3*(*p).M2)*(r*r - 4*z*z));
    double term1 = 4*(powl(r, 3)*8*powl(pyth, 3) + 2*r*powl(pyth, 4))*(-(*p).J*(*p).J*16*z - (*p).M*(powl((*p).M, 3) + 3*(*p).M2)*8*z);
    double term2 = (*p).J*(*p).J*(r*r - 8*z*z) + (*p).M*(powl((*p).M, 3) + 3*(*p).M2)*(r*r - 4*z*z);
    double term3 = - 32*z*(2*r*powl(pyth, 3)*term2 + powl(r, 3)*6*powl(pyth, 2)*term2 + 2*powl(r, 3)*powl(pyth, 3)*((*p).J*(*p).J + (*p).M*(powl((*p).M, 3) + 3*(*p).M2)));
    double term4 = ((term1 + term2)*powl(pyth, 8) - 256*num1*powl(pyth, 7)*r)/(256*powl(pyth, 16));
    double num2 = 8*r*r*z*pyth;
    double dnum2_r = 8*z*(2*r*pyth + 2*powl(r, 3));
    double term5 = (dnum2_r*4*powl(pyth, 4) - 32*r*pow(pyth, 3)*num2)/(16*powl(pyth, 8));
    return term4 + term5;
}

double ddgamma_z_z(double r, double z, Params *p) {
    double pyth = pythagorean(r, z, p);
    double term1 = pyth*(-2*(*p).J*(*p).J - (*p).M*(powl((*p).M, 3) + 3*(*p).M2));
    double term2 = (*p).J*(*p).J*(r*r - 8*z*z) + (*p).M*(powl((*p).M, 3) + 3*(*p).M2)*(r*r - 4*z*z);
    double dterm1 = 2*z*(-2*(*p).J*(*p).J - (*p).M*(powl((*p).M, 3) + 3*(*p).M2));
    double dterm2 = -16*z*(*p).J*(*p).J - 8*z*(*p).M*(powl((*p).M, 3) + 3*(*p).M2);
    double term3 = 2*r*r*z*(((dterm1 - dterm2)*pyth - 10*z*(term1 - term2))/powl(pyth, 6) + (term1 + term2)/powl(pyth, 5));
    double term4 = 2*r*r*(pyth - 6*z*z)/powl(pyth, 4);
    return term3 + term4;
}


double g_tt(double r, double z, Params *p) { return - f(r, z, p); }
double dg_tt_r(double r, double z, Params *p) { return - df_r(r, z, p); }
double dg_tt_z(double r, double z, Params *p) { return - df_z(r, z, p); }
double ddg_tt_r_r(double r, double z, Params *p) { return - ddf_r_r(r, z, p); }
double ddg_tt_r_z(double r, double z, Params *p) { return - ddf_r_z(r, z, p); }
double ddg_tt_z_z(double r, double z, Params *p) { return - ddf_z_z(r, z, p); }

double g_tf(double r, double z, Params *p) { return omega(r, z, p)*f(r, z, p); }
double dg_tf_r(double r, double z, Params *p) { return omega(r, z, p)*df_r(r, z, p) + domega_r(r, z, p)*f(r, z, p); }
double dg_tf_z(double r, double z, Params *p) { return omega(r, z, p)*df_z(r, z, p) + domega_z(r, z, p)*f(r, z, p); }
double ddg_tf_r_r(double r, double z, Params *p) {
    double term1 = omega(r, z, p)*ddf_r_r(r, z, p) + domega_r(r, z, p)*df_r(r, z, p);
    double term2 = domega_r(r, z, p)*df_r(r, z, p) + ddomega_r_r(r, z, p)*f(r, z, p);
    return term1 + term2;
}
double ddg_tf_r_z(double r, double z, Params *p) {
    double term1 = omega(r, z, p)*ddf_r_z(r, z, p) + domega_z(r, z, p)*df_r(r, z, p);
    double term2 = domega_r(r, z, p)*df_z(r, z, p) + ddomega_r_z(r, z, p)*f(r, z, p);
    return term1 + term2;
}
double ddg_tf_z_z(double r, double z, Params *p) { 
    double term1 = omega(r, z, p)*ddf_z_z(r, z, p) + domega_z(r, z, p)*df_z(r, z, p);
    double term2 = domega_z(r, z, p)*df_z(r, z, p) + ddomega_z_z(r, z, p)*f(r, z, p);
    return term1 + term2;
}

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
double ddg_ff_r_r(double r, double z, Params *p) {
    double f_val = f(r, z, p);
    double df_r_val = df_r(r, z, p);
    double ddf_r_r_val = ddf_r_r(r, z, p);
    double omega_val = omega(r, z, p);
    double domega_r_val = domega_r(r, z, p);

    double num = 2*r*f_val - r*r*df_r_val; 
    double dnum_r = 2*(r*df_r_val + f_val) - (2*r*df_r_val + r*r*ddf_r_r_val);
    double term1 = -2*(df_r_val*omega_val*domega_r_val + f_val*powl(domega_r_val, 2) + f_val*omega_val*ddomega_r_r(r, z, p));
    double term2 = -ddf_r_r_val*powl(omega_val, 2) - 2*df_r_val*omega_val*domega_r_val; 
    double term3 = (dnum_r*powl(f_val, 2) - 2*num*f_val*df_r_val)/powl(f_val, 4);
    return term1 + term2 + term3;
}
double ddg_ff_r_z(double r, double z, Params *p) {
    double f_val = f(r, z, p);
    double df_z_val = df_z(r, z, p);
    double omega_val = omega(r, z, p);
    double domega_r_val = domega_r(r, z, p);
    double domega_z_val = domega_z(r, z, p);

    double num = 2*r*f_val - r*r*df_r(r, z, p);
    double dnum_z = 2*r*df_z_val - r*r*ddf_r_z(r, z, p);
    double term1 = -2*(df_z_val*omega_val*domega_r_val + f_val*domega_r_val*domega_z_val + f_val*omega_val*ddomega_r_z(r, z, p));
    double term2 = -ddf_r_z(r, z, p)*powl(omega_val, 2) - 2*df_r(r, z, p)*omega_val*domega_z_val;
    double term3 = (dnum_z*powl(f_val, 2) - 2*num*f_val*df_z_val)/powl(f_val, 4);
    return term1 + term2 + term3;
}
double ddg_ff_z_z(double r, double z, Params *p) {
    double f_val = f(r, z, p);
    double df_z_val = df_z(r, z, p);
    double ddf_z_z_val = ddf_z_z(r, z, p);
    double omega_val = omega(r, z, p);
    double domega_z_val = domega_z(r, z, p);

    double num = - r*r*df_z_val;
    double dnum_z = - r*r*ddf_z_z_val;
    double term1 = -2*(df_z_val*omega_val*domega_z_val + f_val*powl(domega_z_val, 2) + f_val*omega_val*ddomega_z_z(r, z, p));
    double term2 = -ddf_z_z_val*powl(omega_val, 2) - 2*df_z_val*omega_val*domega_z_val;
    double term3 = (dnum_z*powl(f_val, 2) - 2*num*f_val*df_z_val)/powl(f_val, 4);
    return term1 + term2 + term3;
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
double ddg_rr_r_r(double r, double z, Params *p) {
    double e_2gamma = exp(2*my_gamma(r, z, p));
    double f_val = f(r, z, p);
    double dgamma_r_val = dgamma_r(r, z, p);
    double df_r_val = df_r(r, z, p);

    double num = 2*e_2gamma*dgamma_r_val*f_val - e_2gamma*df_r_val;
    double dnum_r = 2*(2*e_2gamma*powl(dgamma_r_val, 2)*f_val + e_2gamma*ddgamma_r_r(r, z, p)*f_val + e_2gamma*dgamma_r_val*df_r_val) - (2*e_2gamma*dgamma_r_val*df_r_val + e_2gamma*ddf_r_r(r, z, p));
    return (dnum_r*powl(f_val, 2) - 2*num*f_val*df_r_val)/powl(f_val, 4);
}
double ddg_rr_r_z(double r, double z, Params *p) {
    double e_2gamma = exp(2*my_gamma(r, z, p));
    double f_val = f(r, z, p);
    double df_r_val = df_r(r, z, p);
    double df_z_val = df_z(r, z, p);
    double dgamma_z_val = dgamma_z(r, z, p);
    double dgamma_r_val = dgamma_r(r, z, p);
    
    double num = 2*e_2gamma*dgamma_r_val*f_val - e_2gamma*df_r_val;
    double dnum_z = 2*(2*e_2gamma*dgamma_r_val*dgamma_z_val*f_val + e_2gamma*ddgamma_r_z(r, z, p)*f_val + e_2gamma*dgamma_r_val*df_z_val) - (2*e_2gamma*dgamma_z_val*df_r_val + e_2gamma*dff_r_z(r, z, p));
    return (dnum_z*powl(f_val, 2) - 2*num*f_val*df_z_val)/powl(f_val, 4);
}
double ddg_rr_z_z(double r, double z, Params *p) {
    double e_2gamma = exp(2*my_gamma(r, z, p));
    double f_val = f(r, z, p);
    double dgamma_z_val = dgamma_z(r, z, p);
    double df_z_val = df_z(r, z, p);

    double num = 2*e_2gamma*dgamma_z_val*f_val - e_2gamma*df_z_val;
    double dnum_z = 2*(2*e_2gamma*powl(dgamma_z_val, 2)*f_val + e_2gamma*ddgamma_z_z(r, z, p)*f_val + e_2gamma*dgamma_z_val*df_z_val) - (2*e_2gamma*dgamma_z_val*df_z_val + e_2gamma*ddf_z_z(r, z, p));
    return (dnum_z*powl(f_val, 2) - 2*num*f_val*df_z_val)/powl(f_val, 4);
}

double Det(double r, double z, Params *p) {return g_tt(r, z, p)*g_ff(r, z, p) - powl(g_tf(r, z, p), 2); }
double Det_r(double r, double z, Params *p) {
    return g_tt(r, z, p)*dg_ff_r(r, z, p) + dg_tt_r(r, z, p)*g_ff(r, z, p) - 2*g_tf(r, z, p)*dg_rf_r(r, z, p);
}
double Det_z(double r, double z, Params *p) {
    return g_tt(r, z, p)*dg_ff_z(r, z, p) + dg_tt_z(r, z, p)*g_ff(r, z, p) - 2*g_tf(r, z, p)*dg_rf_z(r, z, p);
}


double V_eff(double r, double z, double E, double L_z, Params *p) {
    return 1.0/g_rr(r, z, p)*(1 + (g_ff(r, z, p)*E*E + g_tt(r, z, p)*L_z*L_z + 2*g_tf(r, z, p)*E*L_z)/(Det(r, z, p)));
}

double line_element(double dt, double dphi, double r, double z, Params *p) {
    double term1 = -f(r, z, p)*powl(dt - omega(r, z, p)*dphi, 2);
    double term2 = (exp(2*my_gamma(r, z, p))*(r*r + z*z) + r*r*dphi)/f(r, z, p);
    return sqrt(abs(term1 + term2));
}

void update_g(double r, double z, double g[4][4], Params *p) {
    g[0][0] = g_tt(r, z, p);
    g[0][1] = g_tf(r, z, p);

    g[1][0] = g_tf(r, z, p);
    g[1][1] = g_ff(r, z, p);

    g[2][2] = g_rr(r, z, p);

    g[3][3] = g_rr(r, z, p);
}

void update_g_inv(double r, double z, double g_inv[4][4], Params *p) {
    double Det_val = Det(r, z, p);

    g_inv[0][0] = g_ff(r, z, p)/Det_val;
    g_inv[0][1] = - g_tf(r, z, p)/Det_val;

    g_inv[1][0] = g_inv[0][1];
    g_inv[1][1] = g_tt(r, z, p)/Det_val;

    g_inv[2][2] = 1/g_rr(r, z, p);

    g_inv[3][3] = g_inv[2][2];
}

void update_dg(double r, double z, double dg[4][4][4], Params *p) {
    dg[0][0][2] = dg_tt_r(r, z, p);
    dg[0][1][2] = dg_tf_r(r, z, p);
    dg[1][0][2] = dg_tf_r(r, z, p);
    dg[1][1][2] = dg_ff_r(r, z, p);
    dg[2][2][2] = dg_rr_r(r, z, p);
    dg[3][3][2] = dg_rr_r(r, z, p);
    
    dg[0][0][3] = dg_tt_z(r, z, p);
    dg[0][1][3] = dg_tf_z(r, z, p);
    dg[1][0][3] = dg_tf_z(r, z, p);
    dg[1][1][3] = dg_ff_z(r, z, p);
    dg[2][2][3] = dg_rr_z(r, z, p);
    dg[3][3][3] = dg_rr_z(r, z, p);
}

void update_dg_inv(double r, double z, double dg_inv[4][4][4], Params *p) {
    double Det_val = Det(r, z, p);
    double Det_r_val = Det_r(r, z, p);
    double Det_z_val = Det_z(r, z, p);
    dg_inv[0][0][2] = (dg_ff_r(r, z, p)*Det_val - g_ff(r, z, p)*Det_r_val)/(Det_val*Det_val);
    dg_inv[0][1][2] = - (g_tf_r(r, z, p)*Det_val - g_tf(r, z, p)*Det_r_val)/(Det_val*Det_val);

    dg_inv[1][0][2] = dg_inv[0][1][2];
    dg_inv[1][1][2] = (dg_tt_r(r, z, p)*Det_val - g_tt(r, z, p)*Det_r_val)/(Det_val*Det_val);

    dg_inv[2][2][2] = - dg_rr_r(r, z, p)/powl(g_rr(r, z, p), 2);

    dg_inv[3][3][2] = dg_inv[2][2][2];

    dg_inv[0][0][3] = (dg_ff_z(r, z, p)*Det_val - g_ff(r, z, p)*Det_z_val)/(Det_val*Det_val);
    dg_inv[0][1][3] = - (g_tf_z(r, z, p)*Det_val - g_tf(r, z, p)*Det_z_val)/(Det_val*Det_val);

    dg_inv[1][0][3] = dg_inv[0][1][3];
    dg_inv[1][1][3] = (dg_tt_z(r, z, p)*Det_val - g_tt(r, z, p)*Det_z_val)/(Det_val*Det_val);

    dg_inv[2][2][3] = - dg_rr_z(r, z, p)/powl(g_rr(r, z, p), 2);

    dg_inv[3][3][3] = dg_inv[2][2][3];
    
}

void update_ddg(double r, double z, double ddg[4][4][4][4], Params *p) {
    ddg[0][0][2][2] = ddg_tt_r_r(r, z, p);
    ddg[0][1][2][2] = ddg_tf_r_r(r, z, p);
    ddg[1][0][2][2] = ddg_tf_r_r(r, z, p);
    ddg[1][1][2][2] = ddg_ff_r_r(r, z, p);
    ddg[2][2][2][2] = ddg_rr_r_r(r, z, p);
    ddg[3][3][2][2] = ddg_rr_r_r(r, z, p);
    
    ddg[0][0][3][2] = ddg_tt_r_z(r, z, p);
    ddg[0][1][3][2] = ddg_tf_r_z(r, z, p);
    ddg[1][0][3][2] = ddg_tf_r_z(r, z, p);
    ddg[1][1][3][2] = ddg_ff_r_z(r, z, p);
    ddg[2][2][3][2] = ddg_rr_r_z(r, z, p);
    ddg[3][3][3][2] = ddg_rr_r_z(r, z, p);

    ddg[0][0][2][3] = ddg[0][0][3][2];
    ddg[0][1][2][3] = ddg[0][1][3][2];
    ddg[1][0][2][3] = ddg[1][0][3][2];
    ddg[1][1][2][3] = ddg[1][1][3][2];
    ddg[2][2][2][3] = ddg[2][2][3][2];
    ddg[3][3][2][3] = ddg[3][3][3][2];
    
    ddg[0][0][3][3] = ddg_tt_z_z(r, z, p);
    ddg[0][1][3][3] = ddg_tf_z_z(r, z, p);
    ddg[1][0][3][3] = ddg_tf_z_z(r, z, p);
    ddg[1][1][3][3] = ddg_ff_z_z(r, z, p);
    ddg[2][2][3][3] = ddg_rr_z_z(r, z, p);
    ddg[3][3][3][3] = ddg_rr_z_z(r, z, p);
}

void update_Christoffel_symbols(double r, double z, Params *p, double g[4][4], double g_inv[4][4], double dg[4][4][4], double Christoffel[4][4][4]) {
    memset(Christoffel, 0, 4*4*4*sizeof(double));

    for (int mu = 0; mu < 4; mu++){
        for (int kappa = 0; kappa < 4; kappa++){
            for (int lambda = 0; lambda < 4; lambda++){
                for (int sigma = 0; sigma < 4; sigma++){
                    Christoffel[mu][kappa][lambda] += 0.5*g_inv[mu][sigma]*(dg[sigma][kappa][lambda] + dg[lambda][sigma][kappa] - dg[kappa][lambda][sigma]);
                }
            }
        }
    }
}

void update_DChristoffel_symbols(double r, double z, Params *p, double g[4][4], double g_inv[4][4], double dg_inv[4][4][4], double dg[4][4][4], double ddg[4][4][4][4], double DChristoffel[4][4][4][4]) {
    memset(DChristoffel, 0, 4*4*4*4*sizeof(double));

    for (int mu = 0; mu < 4; mu++){
        for (int kappa = 0; kappa < 4; kappa++){
            for (int lambda = 0; lambda < 4; lambda++){
                for (int nu = 0; nu < 4; nu++) {
                    for (int sigma = 0; sigma < 4; sigma++){
                        DChristoffel[mu][kappa][lambda][nu] += 0.5*dg_inv[mu][sigma][nu]*(dg[lambda][sigma][kappa] + dg[kappa][lambda][sigma] - dg[sigma][kappa][lambda]);
                        DChristoffel[mu][kappa][lambda][nu] += 0.5*g_inv[mu][sigma]*(ddg[sigma][kappa][lambda][nu] + ddg[lambda][sigma][kappa][nu] - ddg[kappa][lambda][sigma][nu]);
                    }
                }
            }
        }
    }

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