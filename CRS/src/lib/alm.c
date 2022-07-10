#include<stdio.h>
#include<stdlib.h>
#include<memory.h>

void almerror(char error_text[])
{
	fprintf(stderr,"%s\n",error_text);
	exit(1);
}

int *ivec(int nv)
{
	int *v;

	v=(int *) calloc(nv,sizeof(int));
	if (!v) almerror("allocation failure in dvec()");
	return v;
}

double *dvec(int nv)
{
	double *v;

	v=(double *) calloc(nv,sizeof(double));
	if (!v) almerror("allocation failure in dvec()");
	return v;
}

double **dmat(int nm,int nv)
{
	int i;
	double **m;

	m=(double **) malloc((size_t)((nm)*sizeof(double*)));
	if (!m) almerror("allocation failure 1 in dmat()");

	for (i=0;i<nm;i++){
		m[i]=(double *) calloc(nv,sizeof(double));
		if (!m[i]) almerror("allocation failure 2 in dten()");
	}

	return m;
}

int **imat(int nm,int nv)
{
	int i;
	int **m;

	m=(int **) malloc((size_t)((nm)*sizeof(int*)));
	if (!m) almerror("allocation failure 1 in dmat()");

	for (i=0;i<nm;i++){
		m[i]=(int *) calloc(nv,sizeof(int));
		if (!m[i]) almerror("allocation failure 2 in dten()");
	}

	return m;
}

double ***dten(int nt,int nm,int nv)
{
	int i,j;
	double ***t;

	t=(double ***) malloc((size_t)((nt)*sizeof(double**)));
	if (!t) almerror("allocation failure 1 in dten()");

	t[0]=(double **) malloc((size_t)((nt*nm)*sizeof(double*)));
	if (!t[0]) almerror("allocation failure 2 in dten()");

	for(i=1;i<nt;i++)
		t[i]=t[i-1]+nm;

	for(i=0;i<nt;i++){
		for(j=0;j<nm;j++){
			t[i][j]=(double *) calloc(nv,sizeof(double));
			if (!t[i][j]) almerror("allocation failure 3 in dten()");
		}
	}
	return t;
}

void freedvec(double *v,int nv)
{
	if (nv!=0)	free(v);
}

void freedmat(double **m,int nm,int nv)
{
	int i;
	if (nv!=0){
		for (i=0;i<nm;i++)	free((double *)m[i]);
		free((double **)m);
	}
}

void freedten(double ***t,int nt,int nm,int nv)
{
	int i,j;
	if (nv!=0){
		for (i=0;i<nt;i++){
			for (j=0;j<nm;j++)	free((double *)t[i][j]);
			free((double **)t[i]);
		}
		free((double ***)t);
	}
}

void dmcpvec(double *d,double *o,int nv)
{
	memcpy(d,o,nv*sizeof(double));
}

void dmcpmat(double **d,double **o,int nm,int nv)
{
	int i;
	for (i=0;i<nm;i++) memcpy(d[i],o[i],nv*sizeof(double));
}

void dmcpten(double ***d,double ***o,int nt,int nm,int nv)
{
	int i,j;
	for (i=0;i<nt;i++)
		for (j=0;j<nm;j++) memcpy(d[i][j],o[i][j],nv*sizeof(double));
}

