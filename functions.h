#ifndef FUNCTIONS_H
#define FUNCTIONS_H

typedef struct {
    double M, J, M2, S3, M4;
} Params;

double*** generate_Christoffel_symbols(double r, double z, const Params* p);
void free_Christoffel(double*** Christoffel);
double** make_g(double r, double z, const Params* p);
void free_g(double** g);
double** make_g_inv(double r, double z, const Params* p);
void free_g_inv(double** g_inv);

#endif
