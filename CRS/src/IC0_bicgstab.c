/*
 *	このプログラムは，Bi-CGstab法のテストプログラムです．
 *	解くべき行列は，正定値対称行列の二次元の楕円方程式（ex. 二次元平板の熱伝熱方程式）を想定しています．
 *	前処理行列には，不完全コレスキー分解を用いています．
 *
 */
#include <stdio.h>
#include "cg_solver.h"

static int N = 10000;//行×積
static int M = 100;//行の要素数
static int Kmax = 100000;
static double EPS = 1.0e-12;
static double EPSmax = 1.0e50;
static double val_cen = 6.0;
static double val_nei = -1.0;
Vector *b, *x0, *r0, *x, *s, *r, *d;
static CRS *A, *P, *R;
static Vector *xA_t, *res_t;
       
void alm()
{   
    b = vec_init(N);
    x0 = vec_init(N);
    r0 = vec_init(N);
    x = vec_init(N);
    s = vec_init(N);
    r = vec_init(N);
    d = vec_init(N);
    xA_t = vec_init(N);
    res_t = vec_init(N);
}

void Bi_CGstab_init()
{
    double c_res;
    //FILE *fp;

    //fp = fopen("matrix.out","w");
    FILLIN l = *fillin_init(&l);

    A = crs_makePositiveDefiniteMatrix(N,M,val_cen,val_nei);
    P = IC0(A);
    R = IC0_residual(A,P,&l);

    /*
    crs_write(fp,R);

    c_res = crs_res(R);
    printf("ave.c_res = %lf\n",c_res/l.len);
    */

    if(N<30)
    {
        printf("A.");
        crs_print(A);

        printf("P.");
        crs_print(P);

        printf("R.");
        crs_print(R);
    }

    alm();
}

void Bi_CGstab_solver()
{
    int k;
    double alpha, beta, omega;
    double normr0rk, normr0rk1;
    double bb;
    double err;

    Vector *Ax0 = vec_init(N);
    Vector *P_1r = vec_init(N);
    Vector *Ad = vec_init(N);
    Vector *aAd = vec_init(N);
    Vector *P_1Ad = vec_init(N);
    Vector *aP_1Ad = vec_init(N);
    Vector *P_1s = vec_init(N);
    Vector *AP_1s = vec_init(N);
    Vector *ad = vec_init(N);
    Vector *oP_1s = vec_init(N);
    Vector *oAP_1s = vec_init(N);
    Vector *oP_1Ad = vec_init(N);

    k = 0;

    vec_makeZero(x0);
    vec_makeRandom(b,0,5);
    vec_normCrs(Ax0, A, x0);
    vec_minus(r, b, Ax0);
    IC0_forward_substitution(P_1r,P,r);
    vec_copy(r0,P_1r);
    vec_copy(d,r0);
    bb = vec_norm(b,b);
    normr0rk1 = vec_norm(r0, P_1r);

    printf("start iteration\n");
    do
    {
        k += 1;
        // ①　
        normr0rk = normr0rk1;
        IC0_forward_substitution(P_1Ad,P,vec_normCrs(Ad,A,d));
        alpha = normr0rk / vec_norm(r0, P_1Ad);
        // ②
        vec_minus(s,r,vec_scala(aAd, Ad, alpha));
        // ③
        vec_minus(P_1s, P_1r, vec_scala(aP_1Ad,P_1Ad,alpha));
        vec_normCrs(AP_1s, A, P_1s);
        omega = vec_norm(AP_1s, s) / vec_norm(AP_1s, AP_1s);
        // ④
        vec_scala(ad, d, alpha);
        vec_scala(oP_1s, P_1s, omega);
        vec_plus(x, x, ad);
        vec_plus(x, x, oP_1s);
        // ⑤
        vec_minus(r, s, vec_scala(oAP_1s, AP_1s, omega));
        IC0_forward_substitution(P_1r, P, r);
        normr0rk1 = vec_norm(r0, P_1r);
        // ⑦
        beta = (alpha * normr0rk1) / (omega * normr0rk);
        // ⑧
        vec_scala(oP_1Ad, P_1Ad, omega);
        vec_minus(oP_1Ad, d, oP_1Ad);
        vec_plus(d, P_1r, vec_scala(oP_1Ad, oP_1Ad, beta));

        err = vec_norm(r,r) / bb;
        if(k % 100 == 0)
        {
            printf("k = %6d, err = %le\n", k, err);
        }

    } while (err>EPS && k < Kmax && err < EPSmax);
    printf("iteration end!\n");

    if(err > EPS)
    {
        printf("Not Converged!!\n");
    }
    printf("k = %6d, err = %le\n", k, err); 
}

void check()
{
    double res;

    xA_t = vec_normCrs(xA_t,A,x);
    res_t = vec_minus(res_t,b,xA_t);
    res = vec_norm(res_t,res_t);

    printf("error = %lf\n",res);
}

int main(void)
{
    Bi_CGstab_init();
    Bi_CGstab_solver();
    check();
}
