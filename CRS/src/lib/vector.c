#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include "alm.h"
#include "vector.h"

Vector *vec_init(int l)
{
    Vector *self = malloc(sizeof *self);
    if(self == NULL) {exit(1);}

    self->len = l;
    self->val = dvec(l);

    return self;
}

Vector *vec_read(FILE *fp)
{
    int i;
    int l;
    double m;

    fscanf(fp,"%d\n",&l);
    Vector *self = vec_init(l);

    for(i=0;i<self->len;i++)
    {
        fscanf(fp,"%lf\n",&m);
        self->val[i] = m;
    }

    return self;
}

double vec_norm(Vector *a, Vector *b)
{
    int l;
    double sum = 0;
    #pragma omp parallel
    {
    #pragma omp for reduction (+:sum)
    for(l=0;l<(a->len);l++)
    {
        sum += a->val[l] * b->val[l];
    }
    }
    return sum;
}

Vector *vec_normCrs(Vector *self, CRS *A, Vector *x)
{
    if(A->N != x->len) { printf("error!\tA->N != a->len\n"); exit(1); }
    if(self == NULL)
    {
        Vector *self = vec_init(x->len);
    }

    int i,j;
    #pragma omp parallel
    {
    #pragma omp for private(j)
    for(i=0;i<(A->N);i++)//並列するなら，行ごとに
    {
        self->val[i] = 0.0;
        for(j=(A->row_ptr[i]);j<(A->row_ptr[i+1]);j++)
        {
            self->val[i] += A->val[j]*x->val[A->col_ind[j]];
        }
    }
    }

    return self;
}

Vector *vec_plus(Vector *self, Vector *a, Vector *b)
{
    if(a->len != b->len) { printf("error!\tA->N != a->len\n"); exit(1); }

    int i;
    #pragma omp parallel
    {
    #pragma omp for
    for(i=0;i<a->len;i++)
    {
        self->val[i] = a->val[i] + b->val[i];
    }
    }

    return self;
}

Vector *vec_minus(Vector *self, Vector *a, Vector *b)
{
    if(a->len != b->len) { printf("error!\tA->N != a->len\n"); exit(1); }

    int i;
    #pragma omp parallel
    {
    #pragma omp for
    for(i=0;i<a->len;i++)
    {
        self->val[i] = a->val[i] - b->val[i];
    }
    }
    return self;
}

Vector *vec_copy(Vector *self, Vector *a)
{
    if(self == NULL)
    {
        Vector *self = vec_init(a->len);
    } 

    int i;
    #pragma omp parallel
    {
    #pragma omp for
    for(i=0;i<a->len;i++)
    {
        self->val[i] = a->val[i];
    }
    }

    return self;
}

Vector *vec_scala(Vector *self, Vector *x, double a)
{
    int i;
    #pragma omp parallel
    {
    #pragma omp for
    for(i=0;i<(self->len);i++)
    {
        self->val[i] = a * x->val[i]; 
    }
    }

    return self;
}

void vec_print(Vector *a)
{
    int i;

    printf("val:\n[\n");

    for(i=0;i<(a->len);i++)
    {
       printf("%lf\n",a->val[i]);
    }
    printf("]\n");
}
 
Vector *vec_makeZero(Vector *self)
{
    int i;

    for(i=0;i<self->len;i++)
    {
        self->val[i] = 0.0;
    }

    return self;
}

Vector *vec_makeRandom(Vector *self, double min, double max)
{
    int i;

    for(i=0;i<self->len;i++)
    {
        self->val[i] = min + (rand() * (max - min + 1.0) / (1.0 + RAND_MAX));
    }

    return self;
}