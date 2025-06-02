#ifndef FUNCTIONS_H
#define FUNCTIONS_H

typedef struct {
    double M, J, M2, S3, M4;
} Params;

double*** generate_Christoffel_symbols(double r, double z, Params* p, double** g, double** g_inv, double*** dg);
void free_Christoffel(double*** Christoffel);

double*** make_dg();
void free_dg(double*** dg);

double** make_g();
void update_g(double r, double z, double** g, Params* p);
void free_g(double** g);

double** make_g_inv();
void update_g_inv(double r, double z, double** g_inv, Params *p);
void free_g_inv(double** g_inv);

#endif
