void almerror(char error_text[]);
int *ivec(int nv);
double *dvec(int nv);
double **dmat(int nm,int nv);
int **imat(int nm,int nv);
double ***dten(int nt,int nm,int nv);
void freedvec(double *v,int nv);
void freedmat(double **m,int nm,int nv);
void freedten(double ***t,int nt,int nm,int nv);
void dmcpvec(double *d,double *o,int nv);
void dmcpmat(double **d,double **o,int nm,int nv);
void dmcpten(double ***d,double ***o,int nt,int nm,int nv);

