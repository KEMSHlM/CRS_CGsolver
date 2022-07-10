#ifndef CRS_H
#define CRS_H

#include<stdio.h>

/* CRS形式の行列の構造体
ーーーーーーーーーーーーーーーーーーーーーーーーーーーー */
typedef struct 
{
    const int N; //次元数
    int nnz; //要素数
    int* row_ptr; //各列の最初の要素番号
    int* diag_ptr; //対角成分の要素番号
    int* col_ind; //各要素の列番号
    double* val; //各要素の値
}CRS;

/* CRS形式の初期化
ーーーーーーーーーーーーーーーーーーーーーーーーーーーー */
// 引数に，次元数と要素数を受け取り，CRS形式の行列を作成する
CRS *crs_make(int n, int nnz);
// 行列データから，CRS形式の行列を作成する
CRS *crs_read(FILE *fp);
// CRS形式の値のみのコピー
CRS *crs_copy(CRS *A);
// メモリの解放
void crs_free(CRS* A);
// メモリ領域を変更する
int crs_realloc(CRS* A, int n);
// (n,n+m)番目を空白にし，右にシフトする
int crs_shift(CRS* A, int n, int m);

/*　CRS形式の基本演算
ーーーーーーーーーーーーーーーーーーーーーーーーーーーー */
// 行番号と列番号を指定し，CRS形式での要素番号を返す
int find_ptr(int row, int col, CRS* A);
// 行番号と列番号を指定し，右方向に一番近い要素番号を返す
int find_rightPtr(int row, int col, CRS *A);
// 行番号と列番号を指定し，左方向に一番近い要素番号を返す
int find_leftPtr(int row, int col, CRS *A);

/* デバッグ用
ーーーーーーーーーーーーーーーーーーーーーーーーーーーー */
// データの書き出し
void crs_write(FILE* fp, CRS* A);
// コンソール画面への書き出し
void crs_print(CRS* A);
// 正定値対称行列を生成
CRS *crs_makePositiveDefiniteMatrix(int n, int m, double val_cen, double val_nei);

/* テスト中
ーーーーーーーーーーーーーーーーーーーーーーーーーーーー */
// 行番号と列番号を指定し，下方向に一番近い要素番号を返す
int find_downPtr(int row, int col, CRS *A);
// 行番号と列番号を指定し，上方向に一番近い要素番号を返す
int find_upPtr(int row, int col, CRS *A);


#endif 
