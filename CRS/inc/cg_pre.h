#ifndef CG_PRE_H
#define CG_PRE_H

#include "crs.h"
#include "fillin.h"
#include "vector.h"

/* 前処理行列
ーーーーーーーーーーーーーーーーーーーーーーーーーーーー */
/* 不完全LU分解 */
CRS *ILU0(CRS *A);
// x = L^(-1)r (ILU0用前進代入)
Vector *ILU0_forward_substitution(Vector *self, CRS *L, Vector *r);
// x = U^(-1)r (ILU0用後退代入)
Vector *ILU0_backward_substitution(Vector *self, CRS *U, Vector *r);

/* 不完全コレスキー分解 */
CRS *IC0(CRS *A);
// 転置行列の代入
#define UT 0 //上三角行列の転置
#define LT 1 //下三角行列の転置
void IC0_transpose(CRS *P, int mode);
// x = L^(-1)r (IC0用前進代入)
Vector *IC0_forward_substitution(Vector *x, CRS *L, Vector *r);
// x = U^(-1)r (IC0用後退代入)
Vector *IC0_backward_substitution(Vector *x, CRS *U, Vector *r);

/* 残差およびfillinの計算 
ーーーーーーーーーーーーーーーーーーーーーーーーーーーー */
// R = A - LU (IC0用残差計算)
CRS *ILU0_residual(CRS *A, CRS *P, FILLIN *r);
// R = A - LL^(T) (ILU0用残差計算)
CRS *IC0_residual(CRS *A, CRS *P, FILLIN *r);
// fillinの計算
CRS *crs_fillin(CRS *self, CRS *P, FILLIN *r);
// 前処理行列Pから，fillinの要素を取得
FILLIN *fillin_getList(FILLIN *self, CRS *P);
// fillinの値の計算(対角成分の計算は仮定していない)
double pre_calVal(CRS *P, int row, int col);

/* デバッグ用 
ーーーーーーーーーーーーーーーーーーーーーーーーーーーー */
// fillinの評価
double fillin_res(FILLIN *r, CRS *P);
// 残差の評価
double crs_res(CRS *R);

/* テスト中
ーーーーーーーーーーーーーーーーーーーーーーーーーーーー */

#endif