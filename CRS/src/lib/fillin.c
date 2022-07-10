#include <stdlib.h>
#include <stdio.h>
#include "fillin.h"

FILLIN *fillin_init(FILLIN *self)
{
    self->len = 0;
    self->point = NULL;

    return self;
}

FILLIN *fillin_addPoint(FILLIN *self, POINT a)
{
    int m;

    if(self->len == 0)
    {
        self->point = calloc(1, sizeof(POINT));
        self->point[0] = a;
        self->len = 1;
    }
    else if(self->len == 1)
    {   
        if(point_comparison(self->point[0],a) == 1)
        {
            fillin_realloc(self, 1);
            self->point[1] = a;
        }
        else
        {
            fillin_shift(self, 1, 0);
            self->point[0] = a; 
        }
    }
    else
    {   
        m = fillin_bsrshPoint(self, a);
        fillin_shift(self, 1, m);
        self->point[m] = a;
    }

    return self;
}

int fillin_realloc(FILLIN *x, int n)
{
    int l;

    l = x->len + n;
    if(l < 0) { printf("error!\tminus\n");exit(1); }
    x->point = realloc(x->point, l*sizeof(POINT));
    x->len = l;

    return 0;
}

int fillin_shift(FILLIN *x, int n, int m)
{
    int i;

    fillin_realloc(x,n);
    if(m > (x->len-1) || m < 0) return 0;
    for(i=(x->len-1);i>=m;i--)
    {
        x->point[i] = x->point[i-n];
    }

    return 0;
}

int fillin_free(FILLIN *x)
{
    free(x -> point);

    x -> point = NULL;
}

int fillin_bsrshPoint(FILLIN *self, POINT a)
{
    int mid;
    int left, right;

    left = 0;
    right = self->len - 1;
    if(point_comparison(a,self->point[left]) == 1) return left;
    if(point_comparison(self->point[right],a) == 1) return right + 1;
    while((right - left) != 1)
    {   
        mid = (left + right) / 2;
        if(point_comparison(self->point[mid],a) == 1)
        {
            left = mid;
        }
        else
        {
            right = mid;
        }
    }
    return right;
}

int point_comparison(POINT a, POINT b)
{
    if(a.row < b.row)
    {
        return 1;
    }
    else if (a.row == b.row)
    {
        if(a.col < a.col)
        {
            return 1;
        }   
    }
    return 0;
}

void fillin_print(FILLIN *r)
{
    int i;
    printf("len = %d\n",r->len);
    for(i=0;i<r->len;i++)
    {
        printf("point[%d]\t| row=%d col=%d\n",i,r->point[i].row, r->point[i].col);   
    }
}