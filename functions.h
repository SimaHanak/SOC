#ifndef FUNCTIONS_H
#define FUNCTIONS_H

typedef struct {
    double M, J, M2, S3, M4;
} Params;

void update_Christoffel_symbols(double r, double z, Params *p, double g[4][4], double g_inv[4][4], double dg[4][4][4], double Christoffel[4][4][4]);

void update_DChristoffel_symbols(double r, double z, Params *p, double g[4][4], double g_inv[4][4], double dg_inv[4][4][4], double dg[4][4][4], double ddg[4][4][4][4], double DChristoffel[4][4][4][4]);

void update_ddg(double r, double z, double ddg[4][4][4][4], Params *p);

void update_dg(double r, double z, double dg[4][4][4], Params* p);

void update_dg_inv(double r, double z, double dg_inv[4][4][4], Params* p);

void update_g(double r, double z, double g[4][4], Params *p);

void update_g_inv(double r, double z, double g_inv[4][4], Params *p);

double V_eff(double r, double z, double E, double L_Z, Params *p);

double line_element(double dt, double dphi, double r, double z, Params *p);

#endif
