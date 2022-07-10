#ifndef VEC_H
#define VEC_H

#include "crs.h"

/* Vectorの配列の構造体
ーーーーーーーーーーーーーーーーーーーーーーーーーーーー */
typedef struct 
{
    int len; //次元数
    double* val; //各要素の値
}Vector;

// ベクター配列のメモリ確保
Vector *vec_init(int l);
// ベクター配列の値のみコピー
Vector *vec_copy(Vector *self, Vector *a);
// ベクター配列の内積
double vec_norm(Vector *a, Vector *b);
// self = A・x (ベクター配列とCRS形式の積)
Vector *vec_normCrs(Vector *self, CRS *A, Vector *x);
// self = a + b (ベクター配列同士の和)
Vector *vec_plus(Vector *self, Vector *a, Vector *b);
// self = a - b (ベクター配列同士の差)
Vector *vec_minus(Vector *self, Vector *a, Vector *b);
// self = a * self (ベクター配列のスカラー倍)
Vector *vec_scala(Vector *self, Vector *x, double a);

/* デバッグ用 
ーーーーーーーーーーーーーーーーーーーーーーーーーーーー */
// コンソールに出力
void vec_print(Vector *a);
// 零ベクトルを作成
Vector *vec_makeZero(Vector *self);
// ランダムでベクトルを作成
Vector *vec_makeRandom(Vector *self, double min, double max);

#endif