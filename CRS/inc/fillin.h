#ifndef FILLIN_H
#define FILLIN_H

/* 残差計算に用いる構造体
ーーーーーーーーーーーーーーーーーーーーーーーーーーーー */
typedef struct 
{
   int row;
   int col;
}POINT;

typedef struct
{
    POINT* point;
    int len;
}FILLIN;

/* FILLIN配列の基本操作 
ーーーーーーーーーーーーーーーーーーーーーーーーーーーー */
// FILLIN配列の初期化
FILLIN *fillin_init(FILLIN *self);
// FILLIN配列をn個拡張
int fillin_realloc(FILLIN *x, int n);
// FILLIN配列をm番目からn個シフトする
int fillin_shift(FILLIN *x, int n, int m);
// FILLIN配列のメモリ領域を解放
int fillin_free(FILLIN *x);
//　コンソールへfillinを出力する
void fillin_print(FILLIN *r);
// POINT追加位置を二分探索法
int fillin_bsrshPoint(FILLIN *self, POINT a);
// aよりbの方が後に来る場合，TRUE(=1)を返す
int point_comparison(POINT a, POINT b);
// FILLIN配列にPOINTを追加する
FILLIN *fillin_addPoint(FILLIN *self, POINT a);

#endif