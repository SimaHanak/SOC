#ifndef FUNCTIONS_H
#define FUNCTIONS_H

typedef struct {
    long double M, J, M2, S3, M4;
} Params;

long double*** generate_Christoffel_symbols(long double r, long double z, Params* p, long double** g, long double** g_inv, long double*** dg);
void free_Christoffel(long double*** Christoffel);

long double*** make_dg();
void update_dg(long double r, long double z, long double*** dg, Params* p);
void free_dg(long double*** dg);

long double** make_g();
void update_g(long double r, long double z, long double** g, Params* p);
void free_g(long double** g);

long double** make_g_inv();
void update_g_inv(long double r, long double z, long double** g_inv, Params *p);
void free_g_inv(long double** g_inv);

#endif
