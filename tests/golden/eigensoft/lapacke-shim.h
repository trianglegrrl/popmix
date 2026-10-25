/* Prototypes for the four LAPACKE routines ksrc/kjg_gsl.c calls (fastmode
   only); the system has liblapacke.so.3 but not its development header. */
#ifndef LAPACKE_SHIM_H
#define LAPACKE_SHIM_H
typedef int lapack_int;
#define LAPACK_ROW_MAJOR 101
#define LAPACK_COL_MAJOR 102
double LAPACKE_dlange(int matrix_layout, char norm, lapack_int m, lapack_int n, const double *a, lapack_int lda);
lapack_int LAPACKE_dgeqrf(int matrix_layout, lapack_int m, lapack_int n, double *a, lapack_int lda, double *tau);
lapack_int LAPACKE_dorgqr(int matrix_layout, lapack_int m, lapack_int n, lapack_int k, double *a, lapack_int lda, const double *tau);
lapack_int LAPACKE_dgesvd(int matrix_layout, char jobu, char jobvt, lapack_int m, lapack_int n, double *a, lapack_int lda, double *s, double *u, lapack_int ldu, double *vt, lapack_int ldvt, double *superb);
#endif
