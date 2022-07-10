# CRS_CGsolver
This is a conjugate gradient method program for solving simultaneous linear equations.

# 1. 概要
流体解析(CFD)や構造解析(CAE)で解くべき行列は，大規模かつ疎(スパース)である．そこでCRS(CSR)形式でデータを格納することでメモリ消費を節約できる．CRS形式の行列は，その特殊なデータ格納形態から行列の演算が複雑になる．そこで，本プログラムではCRS形式の行列の基本演算に加え，不完全LU分解(ILU(0))と不完全コレスキー分解(IC(0))の前処理行列生成を実装した．これら基本演算と前処理行列を組み合わせることで，共役勾配法を構築することが可能である．  
使用可能な関数の説明は，各ヘッダーファイルに記載している．  

>[URL : CRSについて](https://zenn.dev/hishinuma_t/books/sparse-matrix-and-vector-product/viewer/crs)  
>[URL : 不完全コレスキー分解](https://cattech-lab.com/science-tools/lecture-mini-preconditioned-matrix/#%E4%B8%8D%E5%AE%8C%E5%85%A8%E3%82%B3%E3%83%AC%E3%82%B9%E3%82%AD%E3%83%BC%E5%88%86%E8%A7%A3)
>[URL : 不完全LU分解](https://cattech-lab.com/science-tools/lecture-mini-preconditioned-matrix/#%E4%B8%8D%E5%AE%8C%E5%85%A8LU%E5%88%86%E8%A7%A3) 

# 2.　コンパイル
コンパイルには，[Makefile](./CRS/src/makefile)を用いる．Makefileには，デバッグ，リリース，並列のオプションを実装した．

# 3. テストプログラム
テストプログラムに，Bi-CGstab法を実装した．従来のCG法は，正定値対称行列のみ適用可能であった．そこで，広い適用範囲と，より強固に収束する性質を持つBi-CGstab法が開発された．  
テストプログラムは，前処理行列に不完全コレスキー分解を用いた[IC0_bicgstab.c](./CRS/src/IC0_bicgstab.c)と前処理行列に不完全LU分解を用いた[ILU0_bicgstab.c](./CRS/src/ILU0_bicgstab.c)の二つを用意した.使い分けには，[Makefile](./CRS/src/makefile) の``$(TARGET)``を変更してほしい．
どちらも初期設定では，二次元の楕円方程式(ex. 二次元平板の熱伝導方程式)を解くことを仮定している．

>[URL : Bi-CGstab法](https://www2.ccs.tsukuba.ac.jp/workshop/HPCseminar/2011/material/2011-04-linear-system.pdf)

# 4. 今後の実装について
本プログラムは，流体解析に用いるために開発した．今後は，fillinを考慮した前処理行列や，不完全LU分解や不完全コレスキー分解の並列化を実装したく思う．