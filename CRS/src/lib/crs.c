#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <omp.h>
#include "alm.h"
#include "crs.h"

/* 構造体のメンバーにconst型を用いるための，特別な処理．*/
CRS *crs_make(int n, int nnz)
{
    CRS init = {.N = n};
    CRS *self = malloc(sizeof *self);
    if(self == NULL) {exit(1);}

    memcpy(self, &init, sizeof *self);
    self->nnz = nnz;
    self->row_ptr=ivec(n+1);
    self->diag_ptr=ivec(n+1);
    self->col_ind = ivec(nnz);
    self->val=dvec(nnz);

    return self;
}

CRS *crs_copy(CRS *A)
{
    int i, ptr;
    int n, nnz;

    n = A->N;
    nnz = A->nnz;

    CRS *B = crs_make(n,nnz);

    //#pragma omp parallel
    //{
    //#pragma omp for
    for(i=0;i<n+1;i++)
    {
        B->row_ptr[i] = A->row_ptr[i];
        B->diag_ptr[i] = A->diag_ptr[i];
    }
    //#pragma omp for
    for(ptr=0;ptr<nnz;ptr++)
    {
        B->col_ind[ptr] = A->col_ind[ptr];
        B->val[ptr] = A->val[ptr];
    }
    //}

    return B;
}

void crs_free(CRS* A)
{
    free(A -> row_ptr);
    free(A -> diag_ptr);
    free(A -> col_ind);
    free(A -> val);

    A -> row_ptr = NULL;
    A -> diag_ptr = NULL;
    A -> col_ind = NULL;
    A -> val = NULL;
}

int crs_realloc(CRS* A, int n)
{
    int l;

    l = A->nnz + n;
    if(l < 0) { printf("error!\tminus\n");exit(1); }
    A->nnz += n;
    A->col_ind = realloc(A->col_ind, l*sizeof(CRS));
    A->val = realloc(A->val, l*sizeof(CRS));

    return 0;
}

int crs_shift(CRS* A, int n, int m)
{
    int i;
    int l;

    crs_realloc(A, n);//できる
    
    
    for(i=A->nnz-1;i>=m;i--)
    {
        A->col_ind[i] = A->col_ind[i-n];
        A->val[i] = A->val[i-n];
    }
    #pragma omp parallel 
    {
    #pragma omp for
    for(i=0;i<=A->N;i++)
    {
        if(A->diag_ptr[i]>m)
        {
            A->diag_ptr[i] += n;
        }
        if(A->row_ptr[i]>m)
        {
            A->row_ptr[i] += n;
        }
    }
    }

    return 0;
}

/* 検索結果がない場合，-1を返す．*/
int find_ptr(int row, int col, CRS* A)
{
    int i;
    int ptr = -1;

    //行列の範囲外を検索した場合に異常終了する
    if((row >= A->N)||(col >= A->N))
    { 
        printf("error!\trow=%d, col=%d, A->N=%d\n",row,col,A->N); 
        exit(1);
    }

    for(i=A->row_ptr[row];i<A->row_ptr[row+1];i++)
    {
        if(A->col_ind[i] == col)
        {
            ptr = i;
            break;
        }
    }

    return ptr;
}

int find_rightPtr(int row, int col, CRS* A)
{
    int i;
    int ptr = -1;

    //行列の範囲外を検索した場合に異常終了する
    if((row >= A->N)||(col >= A->N))
    { 
        printf("error!right_ptr\trow=%d, col=%d, A->N=%d\n",row,col,A->N); 
        exit(1);
    }

    for(i=A->row_ptr[row];i<A->row_ptr[row+1];i++)
    {
        if(A->col_ind[i] > col)
        {
            ptr = i;
            break;
        }
    }

    return ptr;
}

int find_leftPtr(int row, int col, CRS* A)
{
    int i;
    int ptr = -1;

    //行列の範囲外を検索した場合に異常終了する
    if((row >= A->N)||(col >= A->N))
    { 
        printf("error!left_ptr\trow=%d, col=%d, A->N=%d\n",row,col,A->N); 
        exit(1);
    }

    for(i=A->row_ptr[row];i<A->row_ptr[row+1];i++)
    {
        if(A->col_ind[i+1] > col)
        {
            ptr = i;
            break;
        }
    }

    return ptr;
}

int find_downPtr(int row, int col, CRS *A)
{
    int i;
    int ptr = -1;
    int ptr2;

    for(ptr2=A->row_ptr[row+1];ptr2<(A->nnz);ptr2++)
    {
        if(A->col_ind[ptr2] == col)
        {
            ptr = ptr2;
                break;
        }
    }
    
    return ptr; 
}

int find_upPtr(int row, int col, CRS *A)
{
    int i;
    int ptr = -1;
    int ptr2;

    for(ptr2=A->row_ptr[row]-1;ptr2>=0;ptr2--)
    {
        if(A->col_ind[ptr2] == col)
        {
            ptr = ptr2;
                break;
        }
    }
    
    return ptr; 
}

CRS *crs_read(FILE *fp)
{
    int i, j, k;
    int n, nnz;
    int flg;
    double m;

    fscanf(fp,"%d\t%d\n",&n,&nnz);
    CRS *p = crs_make(n, nnz);

    k = 0;
    flg = 0;
    for(i=0;i<n;i++)
    {
        for(j=0;j<n;j++)
        {
            fscanf(fp,"%lf\t",&m);
            if(m!=0)
            {
                p->val[k] = m;
                p->col_ind[k] = j;
                if(i == j)
                {
                    p->diag_ptr[i] = k;
                }
                if(flg == 0)
                {
                    flg = 1;
                    p->row_ptr[i] = k;
                }
                k += 1;
            }
        }
        flg = 0;
        fscanf(fp,"\n");
    }
    p->diag_ptr[n] = k;
    p->row_ptr[n] = k; 

    return p;
}

void crs_write(FILE* fp, CRS* A)
{
	int i, j, k=0;

	fprintf(fp,"%d\n",A->N);
	for (i=0;i<A->N;i++)
    {
        for(j=0;j<A->N;j++)
        {
            if(A->col_ind[k] == j)
            {
                fprintf(fp,"%lf\t",A->val[k]);
                k += 1;
            }else
            {
		    fprintf(fp,"0.0\t");
            }
        }
        fprintf(fp,"\n");
	}
}

CRS *crs_makePositiveDefiniteMatrix(int n, int m, double val_cen, double val_nei)
{  
    int i;
    int nnz;
    int flg = 0, ptr = 0;
    double val_c = val_cen, val_s = val_nei, val_n = val_nei, val_e = val_nei, val_w = val_nei;

    if(m<3||m>n){ printf("error!\tm must be more than 2.\t"); exit(1); }

    m += 1;
    nnz = 5*n - 2*m;
    CRS *self = crs_make(n,nnz);
    
    for(i=0;i<n;i++)
    {
        // s
        if(i>=(m-1))
        {
            self->val[ptr] = val_s;
            self->col_ind[ptr] = i - (m - 1);
            if(flg == 0)
            {
                self->row_ptr[i] = ptr;
                flg = 1;
            }
            ptr += 1;
        }
        // w
        if(i % (m-1) != 0)
        {
            if(i>=1)
            {
                self->val[ptr] = val_w;
                self->col_ind[ptr] = i - 1;
                if(flg == 0)
                {
                    self->row_ptr[i] = ptr;
                    flg = 1;
                }
                ptr += 1;
            }
        }
        // c
        self->val[ptr] = val_c;
        self->col_ind[ptr] = i;
        self->diag_ptr[i] = ptr;
        if(flg == 0)
        {
            self->row_ptr[i] = ptr;
            flg = 1;
        } 
        ptr += 1;
        // e
        if((i+1) % (m-1) != 0)
        {
            if(i<=n-2)
            {
                self->val[ptr] = val_e;
                self->col_ind[ptr] = i + 1;
                ptr += 1;
            }
        }
        // n
        if(i<= n - m)
        {
            self->val[ptr] = val_n;
            self->col_ind[ptr] = i + (m-1);
            ptr += 1;
        }
        flg = 0;
    }
    self->row_ptr[n] = ptr;
    self->diag_ptr[n] = ptr;

    return self;
}

void crs_print(CRS* A)
{
    int mode = 15015;
    int i, j, k=0;
    
    if(mode % 3 == 0)
    {
        printf("val: [\n");
        for (i=0;i<A->N;i++)
        {
            for(j=0;j<A->N;j++)
            {
                if(A->col_ind[k] == j)
                {
                    printf("%2.2lf\t",A->val[k]);
                    k += 1;
                }else
                {
                printf("0.00\t");
                }
            }
            printf("\n");
        }
        printf("]\n");
    }
    if(mode % 5 == 0)
    {
        printf("ptr:[\n");
        for(i=0;i<A->N;i++)
        {
            for(j=0;j<A->N;j++)
            {
                printf("%d\t",find_ptr(i,j,A));
            } 
            printf("\n");
        }
        printf("]\n");
    } 
    if(mode % 7 == 0)
    {
        printf("col_ind:[\n");
        for(i=0;i<A->N;i++)
        {
            for(j=0;j<A->N;j++)
            {
                if(find_ptr(i,j,A)!=-1)
                {
                    printf("%d\t",A->col_ind[find_ptr(i,j,A)]);
                }else
                {
                    printf("-\t");
                }
            } 
            printf("\n");
        }
        printf("]\n"); 
    }
    if(mode % 11 == 0)
    {
        printf("diag_ptr: [\n");
        for(i=0;i<A->N+1;i++)
        {
            printf("%d\n",A->diag_ptr[i]);
        }
        printf("\n]\n");
    }
    if(mode % 13  == 0)
    {
        printf("row_ptr: [\n");
        for(i=0;i<A->N+1;i++)
        {
            printf("%d\n",A->row_ptr[i]);
        }
        printf("\n]\n");
    }
}
