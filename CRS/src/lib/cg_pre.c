#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cg_pre.h"

Vector *ILU0_forward_substitution(Vector *self, CRS *L, Vector *r)
{
    int i, ptr;
    double sum = 0;

    if(self == NULL)
    {
        self = vec_init(r->len);
    }

    for(i=0;i<L->N;i++)
    {
        for(ptr=L->row_ptr[i];ptr<L->diag_ptr[i];ptr++)
        {
            sum += L->val[ptr] * self->val[L->col_ind[ptr]];
        }    
        self->val[i] = r->val[i] - sum;
        self->val[i] /= L->val[L->diag_ptr[i]]; //<- 上三角行列の対角成分を1とするなら，コメントアウトを解除．
        sum = 0;
    }

    return self;
}

Vector *ILU0_backward_substitution(Vector *self, CRS *U, Vector *r)
{
    int i, ptr;
    double sum = 0;

    if(self == NULL)
    {
        self = vec_init(r->len);
    }

    for(i=(U->N-1);i>=0;i--)
    {
        for(ptr=(U->diag_ptr[i]+1);ptr<(U->row_ptr[i+1]);ptr++)
        {
            sum += U->val[ptr] * self->val[U->col_ind[ptr]];
        }    
        self->val[i] = r->val[i] - sum;
        self->val[i] /= U->val[U->diag_ptr[i]]; //<- 上三角行列の対角成分を1とするなら，コメントアウト．
        sum = 0;
    }

    return self;
}

/* Pの下三角行列と上三角行列が，それぞれLとUに当たる．
上三角行列の対角成分が，1．
*/
CRS *ILU0(CRS *A) //not parallel
{
    CRS *P = crs_copy(A);

    int i, j, k;
    int ptrij, ptrik, ptrkj;

    for(i=1;i<(P->N);i++)
    {
        for(ptrik=P->row_ptr[i];ptrik<P->diag_ptr[i];ptrik++)
        {   
            k = P->col_ind[ptrik];
            P->val[ptrik] /= P->val[P->diag_ptr[k]];
            for(ptrkj=(P->diag_ptr[k]+1);ptrkj<(P->row_ptr[k+1]);ptrkj++)
            {   
                j = P->col_ind[ptrkj];
                ptrij = find_ptr(i,j,P);
                if(ptrij != -1)
                {
                    P->val[ptrij] -= P->val[ptrik] * P->val[ptrkj];
                }
            }
        }
    }

    return P;
}

CRS *ILU0_residual(CRS *A, CRS *P, FILLIN *r)
{
    CRS *R = crs_copy(A);

    int i, j, k;
    int ptri, ptrik, ptrkj;
    double sum = 0.0;

    for(i=0;i<(A->N);i++)
    {
        for(ptri=(P->row_ptr[i]);ptri<(P->row_ptr[i+1]);ptri++)
        {
            j = P->col_ind[ptri];
            if(j>=i)
            {
                for(ptrik=(P->row_ptr[i]);ptrik<=(P->diag_ptr[i]);ptrik++)
                {
                    k = P->col_ind[ptrik];
                    ptrkj = find_ptr(k, j, P);
                    if(ptrkj != -1)
                    {
                        if(k == i)
                        {
                            sum += P->val[ptrkj];
                        }
                        else
                        {
                            sum += P->val[ptrik] * P->val[ptrkj];
                        }
                    } 
                }
            }
            else
            {
                for(ptrik=(P->row_ptr[i]);ptrik<=ptri;ptrik++)
                {
                    k = P->col_ind[ptrik];
                    ptrkj = find_ptr(k, j, P);
                    if(ptrkj != -1)
                    {
                        sum += P->val[ptrik] * P->val[ptrkj];
                    } 
                } 
            }
            
            R->val[ptri] -= sum;
            sum = 0;
        } 
    }
    //fill-inの確認
    R = crs_fillin(R, P, r);
    return R;
}

Vector *IC0_forward_substitution(Vector *self, CRS *L, Vector *r)
{
    int i, ptr;
    double sum = 0;

    if(self == NULL)
    {
        self = vec_init(r->len);
    }

    for(i=0;i<L->N;i++)
    {
        for(ptr=L->row_ptr[i];ptr<L->diag_ptr[i];ptr++)
        {
            sum += L->val[ptr] * self->val[L->col_ind[ptr]];
        }    
        self->val[i] = r->val[i] - sum;
        self->val[i] /= L->val[L->diag_ptr[i]];
        sum = 0;
    }

    return self;
}

Vector *IC0_backward_substitution(Vector *self, CRS *U, Vector *r)
{
    int i, ptr;
    double sum = 0;
    
    if(self == NULL)
    {
        self = vec_init(r->len);
    }

    for(i=(U->N-1);i>=0;i--)
    {
        for(ptr=(U->diag_ptr[i]+1);ptr<(U->row_ptr[i+1]);ptr++)
        {
            sum += U->val[ptr] * self->val[U->col_ind[ptr]];
        }    
        self->val[i] = r->val[i] - sum;
        self->val[i] /= U->val[U->diag_ptr[i]];
        sum = 0;
    }

    return self;
}

void IC0_transpose(CRS *P, int mode)
{   
    int i, j;
    int ptr, ptrji;

    if(mode == 0)
    {
        /* Upper to Lower */
        for(i=1;i<(P->N);i++)
        {
            for(ptr=(P->row_ptr[i]);ptr<(P->diag_ptr[i]);ptr++)
            {
                j = P->col_ind[ptr];
                ptrji = find_ptr(j,i,P);
                P->val[ptr] = P->val[ptrji];
            }
        }
    }
    else if(mode == 1)
    {
        // Lower to Upper (未検証)
        for(i=0;i<(P->N);i++)
        {
            for(ptr=(P->diag_ptr[i]+1);ptr<(P->row_ptr[i+1]);ptr++)
            {
                j = P->col_ind[ptr];
                ptrji = find_ptr(j,i,P);
                P->val[ptr] = P->val[ptrji];
            }
        }
    }
    else
    {
        printf("error!\tMode is uncorected!\n");
        exit(1);
    }
}

/* CRS形式上，上三角行列の計算が早い．*/
CRS *IC0(CRS *A)
{
   CRS *P = crs_copy(A);

    int i, j, k;
    int ptrij, ptrki, ptrkj;
    
    for(k=0;k<(P->N);k++)
    {
        P->val[P->diag_ptr[k]] = sqrt(P->val[P->diag_ptr[k]]);
        for(ptrkj=(P->diag_ptr[k]+1);ptrkj<(P->row_ptr[k+1]);ptrkj++)
        {
            j = P->col_ind[ptrkj];
            P->val[ptrkj] /= P->val[P->diag_ptr[k]];
        }
        for(ptrki=(P->diag_ptr[k]+1);ptrki<(P->row_ptr[k+1]);ptrki++)
        {
            i = P->col_ind[ptrki];
            for(ptrkj=ptrki;ptrkj<(P->row_ptr[k+1]);ptrkj++)//ここ大丈夫か
            {
                j = P->col_ind[ptrkj];
                ptrij = find_ptr(i,j,P);
                if(ptrij != -1)
                {
                    P->val[ptrij] -= P->val[ptrki] * P->val[ptrkj];
                }
            }
        } 
    }

    IC0_transpose(P,UT);

    return P;
}

//　上三角行列の対角成分を非１にしたため，変更が必要かも？
CRS *IC0_residual(CRS *A, CRS *P, FILLIN *r)
{
    CRS *R = crs_copy(A);

    int i, j, k;
    int ptri, ptri2, ptrk, ptrij;
    double sum = 0.0;

    for(i=0;i<(P->N);i++)
    {
        for(ptri=(P->row_ptr[i]);ptri<=(P->diag_ptr[i]);ptri++)
        {
            //対角成分
            if(ptri == P->diag_ptr[i])
            {
                for(ptri2=(P->row_ptr[i]);ptri2<=(P->diag_ptr[i]);ptri2++)
                {
                    sum += P->val[ptri2] * P->val[ptri2];
                }
            }
            else
            {
                k = P->col_ind[ptri];
                for(ptrk=(P->row_ptr[k]);ptrk<=(P->diag_ptr[k]);ptrk++)
                {
                    j = P->col_ind[ptrk];
                    ptrij = find_ptr(i,j,P);
                    if(ptrij != -1)
                    {
                        sum += P->val[ptrk] * P->val[ptrij];
                    }
                }
            }
            R->val[ptri] -= sum;
            sum = 0.0;
        }
    }

    IC0_transpose(R,LT);

    // fillinの計算
    R = crs_fillin(R, P, r);
    return R;
}

CRS *crs_fillin(CRS *self, CRS *P, FILLIN *r)
{
    int l;
    int ptr;

    if(r->len == 0)
    {   
        r = fillin_getList(r,P);
        for(l=0;l<(r->len);l++)
        {   
            ptr = find_rightPtr(r->point[l].row, r->point[l].col, self); 
            crs_shift(self, 1, ptr);
            self->val[ptr] = pre_calVal(P, r->point[l].row, r->point[l].col);
            self->col_ind[ptr] = r->point[l].col;
        }
    }
    else
    {
        for(l=0;l<(r->len);l++)
        {
            ptr = find_ptr(r->point[l].row,r->point[l].col,P);
            self->val[ptr] = pre_calVal(P, r->point[l].row, r->point[l].col);
            self->col_ind[ptr] = r->point[l].col;
        }
    }

    return self;
}

FILLIN *fillin_getList(FILLIN *self, CRS *P)
{
    int i, j, k;
    int ptr, ptr2, ptrkj;
    POINT a1,a2;

    for(i=0;i<(P->N);i++)
    {
        for(ptr=(P->row_ptr[i]);ptr<(P->diag_ptr[i]);ptr++)
        {
            for(j=(P->col_ind[ptr]+1);j<(P->col_ind[ptr+1]);j++)
            {
                for(ptr2=(P->row_ptr[i]);ptr2<=ptr;ptr2++)
                {
                    k = P->col_ind[ptr2];
                    ptrkj = find_ptr(k, j, P);
                    if(ptrkj != -1)
                    {
                        a1.row = i;
                        a1.col = j;
                        self = fillin_addPoint(self,a1);
                        a2.row = P->N - i - 1;
                        a2.col = P->N - j - 1;
                        self = fillin_addPoint(self,a2);
                        break;
                    }
                }
            }
        }
    }
    return self;
}

double pre_calVal(CRS *P, int row, int col)
{
    int i, j ,k;
    int ptr, ptrik, ptrkj;
    double m = 0;

    i = row;
    j = col;

    if(i == j) { printf("error!\tdiag\n"); exit(1); }
    else if(i < j)
    {
        ptrik = find_rightPtr(i,j,P);
        for(ptr=ptrik;ptr<(P->row_ptr[i+1]);ptr++)
        {
            k = P->col_ind[ptr];
            ptrkj = find_ptr(k,j,P);
            m +=  P->val[ptrik] * P->val[ptrkj];
        }
    }
    else
    {
        ptrik = find_leftPtr(i,j,P);
        for(ptr=ptrik;ptr>=(P->row_ptr[i]);ptr--)
        {
            k = P->col_ind[ptr];
            ptrkj = find_ptr(k,j,P);
            m += P->val[ptrik] * P->val[ptrkj];
        }
    }

    return m; 
}

double fillin_res(FILLIN *r, CRS *P)
{
    int i;
    double m;
    double res = 0;

    for(i=0;i<r->len;i++)
    {
        m = pre_calVal(P, r->point[i].row, r->point[i].col);
        m *= m;

        res += m;
    }

    return res;
}

double crs_res(CRS *R)
{
    int ptr;
    double res = 0;

    for(ptr=0;ptr<(R->nnz);ptr++)
    {
        res += R->val[ptr] * R->val[ptr];
    }
    res = sqrt(res);

    return res;
}
