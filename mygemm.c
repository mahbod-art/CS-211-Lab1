#include "mygemm.h"

/**
 * 
 * Implement all functions here in this file.
 * Do NOT change input parameters and return type.
 * 
 **/

//Register Reuse part 1
void dgemm0(const double* A, const double* B, double* C, const int n)
{
    /*int i = 0, j = 0, k = 0;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            for (k = 0; k < n; k++)
            {
                C[i * n + j] += A[i * n + k] * B[k * n + j];
            }
        }
    }*/
}

void dgemm1(const double* A, const double* B, double* C, const int n)
{
    /*int i = 0, j = 0, k = 0;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            register double C_Temp = C[i * n + j];
            for (k = 0; k < n; k++)
            {
                C_Temp += A[i * n + k] * B[k * n + j];
            }
            C[i * n + j] = C_Temp;
        }
    }*/
}
//Register Reuse part 1 End

//Register Reuse part 2
void dgemm2(const double* A, const double* B, double* C, const int n)
{
    int i = 0, j = 0, k = 0;
    for (i = 0; i < n; i += 2)
    {
        for (j = 0; j < n; j += 2)
        {
            register double C1 = C[i * n + j];
            register double C2 = i < (n - 1) ? C[(i + 1) * n + j] : 0;
            register double C3 = j < (n - 1) ? C[i * n + (j + 1)] : 0;
            register double C4 = (i < (n - 1)) && (j < (n - 1))? C[(i + 1) * n + (j + 1)] : 0;
            for (k = 0; k < n; k += 2)
            {
                register double A1 = A[i * n + k];
                register double A2 = i < (n - 1) ? A[(i + 1) * n + k] : 0;
                register double A3 = k < (n - 1) ? A[i * n + (k + 1)] : 0;
                register double A4 = (i < (n - 1)) && (k < (n - 1)) ? A[(i + 1) * n + (k + 1)] : 0;

                register double B1 = B[k * n + j];
                register double B2 = k < (n - 1) ? B[(k + 1) * n + j] : 0;
                register double B3 = j < (n - 1) ? B[k * n + (j + 1)] : 0;
                register double B4 = (k < (n - 1)) && (j < (n - 1)) ? B[(k + 1) * n + (j + 1)] : 0;

                C1 = C1 + A1 * B1 + A3 * B2;
                C2 = C2 + A2 * B1 + A4 * B2;
                C3 = C3 + A1 * B3 + A3 * B4;
                C4 = C4 + A2 * B3 + A4 * B4;
            }

            C[i * n + j] = C1;
            C[(i + 1) * n + j] = C2;
            C[i * n + (j + 1)] = C3;
            C[(i + 1) * n + (j + 1)] = C4;
        }
    }
}

//Register Reuse part 2 End

//Register Reuse part 3
void dgemm3(const double* A, const double* B, double* C, const int n)
{
   /* int i = 0;
    int j = 0;
    int k = 0;
    for (i = 0; i < n; i += 3)
    {
        for (j = 0; j < n; j += 4)
        {
            register double C1 = C[i * n + j];
            register double C2 = i < (n - 1) ? C[(i + 1) * n + j] : 0;
            register double C3 = i < (n - 2) ? C[(i + 2) * n + j] : 0;
            register double C4 = j < (n - 1) ? C[i * n + (j + 1)] : 0;
            register double C5 = i < (n - 1) && j < (n - 1) ? C[(i + 1) * n + (j + 1)] : 0;
            register double C6 = i < (n - 2) && j < (n - 1) ? C[(i + 2) * n + (j + 1)] : 0;
            register double C7 = j < (n - 2) ? C[i * n + (j + 2)] : 0;
            register double C8 = i < (n - 1) && j < (n - 2) ? C[(i + 1) * n + (j + 2)] : 0;
            register double C9 = i < (n - 2) && j < (n - 2) ? C[(i + 2) * n + (j + 2)] : 0;
            register double C10 = j < (n - 3) ? C[i * n + (j + 3)] : 0;
            register double C11 = i < (n - 1) && j < (n - 3) ? C[(i + 1) * n + (j + 3)] : 0;
            register double C12 = i < (n - 2) && j < (n - 3) ? C[(i + 2) * n + (j + 3)] : 0;
            for (k = 0; k < n; k++)
            {
                register double A1 = A[i * n + k];
                register double A2 = i < (n - 1) ? A[(i + 1) * n + k] : 0;
                register double A3 = i < (n - 2) ? A[(i + 2) * n + k] : 0;
                register double B1 = B[k * n + j];

                C1 += A1 * B1;
                C2 += A2 * B1;
                C3 += A3 * B1;

                B1 = j < (n - 1) ? B[k * n + (j + 1)] : 0;
                C4 += A1 * B1;
                C5 += A2 * B1;
                C6 += A3 * B1;

                B1 = j < (n - 2) ? B[k * n + (j + 2)] : 0;
                C7 += A1 * B1;
                C8 += A2 * B1;
                C9 += A3 * B1;

                B1 = j < (n - 3) ? B[k * n + (j + 3)] : 0;
                C10 += A1 * B1;
                C11 += A2 * B1;
                C12 += A3 * B1;
            }

            C[i * n + j] = C1;
            if (i < (n - 1)) C[(i + 1) * n + j] = C2;
            if (i < (n - 2)) C[(i + 2) * n + j] = C3;

            if (j < (n - 1)) C[i * n + (j + 1)] = C4;
            if (i < (n - 1) && j < (n - 1)) C[(i + 1) * n + (j + 1)] = C5;
            if (i < (n - 2) && j < (n - 1)) C[(i + 2) * n + (j + 1)] = C6;

            if (j < (n - 2)) C[i * n + (j + 2)] = C7;
            if (i < (n - 1) && j < (n - 2)) C[(i + 1) * n + (j + 2)] = C8;
            if (i < (n - 2) && j < (n - 2)) C[(i + 2) * n + (j + 2)] = C9;

            if (j < (n - 3)) C[i * n + (j + 3)] = C10;
            if (i < (n - 1) && j < (n - 3)) C[(i + 1) * n + (j + 3)] = C11;
            if (i < (n - 2) && j < (n - 3)) C[(i + 2) * n + (j + 3)] = C12;
        }
    }*/
}
//Register Reuse part 3 End

//Cache Reuse part 3
void ijk(const double *A, const double *B, double *C, const int n) 
{
   /* int i = 0;
    int j = 0;
    int k = 0;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            register double Cij = C[i * n + j];
            for (k = 0; k < n; k++)
            {
                Cij += A[i * n + k] * B[k * n + j];
            }
            C[i * n + j] = Cij;
        }
    }*/
}

void bijk(const double *A, const double *B, double *C, const int n, const int b) 
{
   /* int k=0;
    int j=0;
    int i=0;
    int kB=0;
    int jB=0;
    int iB=0;
    for( i=0;i<n;i+=b)
		for( j=0;j<n;j+=b)	
			for( k=0;k<n;k+=b)
			{
				for( iB=i;iB<i+b && iB<n;iB++)
					for( jB=j;jB<j+b && jB<n;jB++)
					{
						register double res=C[iB*n+jB];
						for( kB=k;kB<k+b && kB<n;kB++)
							res+=A[iB*n+kB]*B[kB*n+jB];
						C[iB*n+jB]=res;
					}
			}*/
}

void jik(const double *A, const double *B, double *C, const int n) 
{
   /* int i = 0;
    int j = 0;
    int k = 0;
    for (j = 0; j < n; j++)
    {
        for (i = 0; i < n; i++)
        {
            register double Cij = C[i * n + j];
            for (k = 0; k < n; k++)
            {
                Cij += A[i * n + k] * B[k * n + j];
            }
            C[i * n + j] = Cij;
        }
    }*/
}

void bjik(const double *A, const double *B, double *C, const int n, const int b) 
{
  /*  int k=0;
    int j=0;
    int i=0;
    int kB=0;
    int jB=0;
    int iB=0;
    for( j=0;j<n;j+=b)
		for( i=0;i<n;i+=b)	
			for( k=0;k<n;k+=b)
			{
				for( jB=j;jB<j+b && jB<n;jB++)
					for( iB=i;iB<i+b && iB<n;iB++)
					{
						register double res=C[iB*n+jB];
						for( kB=k;kB<k+b && kB<n;kB++)
							res+=A[iB*n+kB]*B[kB*n+jB];
						C[iB*n+jB]=res;
					}
			}*/
}

void kij(const double *A, const double *B, double *C, const int n) 
{
  /*  int i = 0;
    int j = 0;
    int k = 0;
    for (k = 0; k < n; k++)
    {
        for (i = 0; i < n; i++)
        {
            register double Aik = A[i * n + k];
            for (j = 0; j < n; j++)
            {
                C[i * n + j] += Aik * B[k * n + j];
            }
        }
    }*/
}

void bkij(const double *A, const double *B, double *C, const int n, const int b) 
{
  /*  int k=0;
    int j=0;
    int i=0;
    int kB=0;
    int jB=0;
    int iB=0;
    for( k=0;k<n;k+=b)
		for( i=0;i<n;i+=b)	
			for( j=0;j<n;j+=b)
			{
				for( kB=k;kB<k+b && kB<n;kB++)
					for( iB=i;iB<i+b && iB<n;iB++)
					{
						register double res=A[iB*n+kB];
						for( jB=j;jB<j+b && jB<n;jB++)
							C[iB*n+jB]+=res*B[kB*n+jB];
					}
			}*/
}


void ikj(const double *A, const double *B, double *C, const int n) 
{
   /* int i = 0;
    int j = 0;
    int k = 0;
    for (i = 0; i < n; i++)
    {
        for (k = 0; k < n; k++)
        {
            register double Aik = A[i * n + k];
            for (j = 0; j < n; j++)
            {
                C[i * n + j] += Aik * B[k * n + j];
            }
        }
    }*/
}

void bikj(const double *A, const double *B, double *C, const int n, const int b) 
{
  /*  int k=0;
    int j=0;
    int i=0;
    int kB=0;
    int jB=0;
    int iB=0;
    for( i=0;i<n;i+=b)
		for( k=0;k<n;k+=b)	
			for( j=0;j<n;j+=b)
			{
				for( iB=i;iB<i+b && iB<n;iB++)
					for( kB=k;kB<k+b && kB<n;kB++)
					{
						register double res=A[iB*n+kB];
						for( jB=j;jB<j+b && jB<n;jB++)
							C[iB*n+jB]+=res*B[kB*n+jB];
					}
			}*/
}

void jki(const double *A, const double *B, double *C, const int n) 
{
  /*  int i = 0;
    int j = 0;
    int k = 0;
    for (j = 0; j < n; j++)
    {
        for (k = 0; k < n; k++)
        {
            register double Bkj = B[k * n + j];
            for (i = 0; i < n; i++)
            {
                C[i * n + j] += A[i * n + k] * Bkj;
            }
        }
    }*/
}

void bjki(const double *A, const double *B, double *C, const int n, const int b) 
{
 /*   int k=0;
    int j=0;
    int i=0;
    int kB=0;
    int jB=0;
    int iB=0;
    for( j=0;j<n;j+=b)
		for( k=0;k<n;k+=b)	
			for( i=0;i<n;i+=b)
			{
				for( jB=j;jB<j+b && jB<n;jB++)
					for( kB=k;kB<k+b && kB<n;kB++)
					{
						register double res=B[kB*n+jB];
						for( iB=i;iB<i+b && iB<n;iB++)
							C[iB*n+jB]+=res*A[iB*n+kB];
					}
			}*/
}

void kji(const double *A, const double *B, double *C, const int n) 
{
 /*   int i = 0;
    int j = 0;
    int k = 0;
    for (k = 0; k < n; k++)
    {
        for (j = 0; j < n; j++)
        {
            register double Bkj = B[k * n + j];
            for (i = 0; i < n; i++)
            {
                C[i * n + j] += A[i * n + k] * Bkj;
            }
        }
    }*/
}

void bkji(const double *A, const double *B, double *C, const int n, const int b) 
{
 /*   int k=0;
    int j=0;
    int i=0;
    int kB=0;
    int jB=0;
    int iB=0;
    for(k=0;k<n;k+=b)
		for(j=0;j<n;j+=b)	
			for(i=0;i<n;i+=b)
			{
				for(kB=k;kB<k+b && kB<n;kB++)
					for(jB=j;jB<j+b && jB<n;jB++)
					{
						register double res=B[kB*n+jB];
						for(iB=i;iB<i+b && iB<n;iB++)
							C[iB*n+jB]+=res*A[iB*n+kB];
					}
			}*/
}
//Cache Reuse part 3 End 

//Cache Reuse part 4
void optimal(const double* A, const double* B, double *C, const int n, const int b)
{
 /*   int i = 0;
    for (i = 0; i < n; i += b)
    {
        int j = 0;
        for (j = 0; j < n; j += b)
        {
            int k = 0;
            for (k = 0; k < n; k += b)
            {
                int i1 = 0;
                for (i1 = i; i1 < (i + b > n? n : (i + b)); i1 += 3)
                {
                    int j1 = 0;
                    for (j1 = j; j1 < (j + b > n? n : (j + b)); j1 += 3)
                    {
                        register double C_0_0 = C[i1 * n + j1];
                        register double C_1_0 = C[(i1 + 1) * n + j1];
                        register double C_2_0 = C[(i1 + 2) * n + j1];

                        register double C_0_1 = C[i1 * n + (j1 + 1)];
                        register double C_1_1 = C[(i1 + 1) * n + (j1 + 1)];
                        register double C_2_1 = C[(i1 + 2) * n + (j1 + 1)];

                        register double C_0_2 = C[i1 * n + (j1 + 2)];
                        register double C_1_2 = C[(i1 + 1) * n + (j1 + 2)];
                        register double C_2_2 = C[(i1 + 2) * n + (j1 + 2)];

                        int k1 = 0;
                        for (k1 = k; k1 < (k + b > n? n : (k + b)); k1++)
                        {
                            register double A_0_M = A[i1 * n + k1];
                            register double A_1_M = A[(i1 + 1) * n + k1];
                            register double A_2_M = A[(i1 + 2) * n + k1];

                            register double B_M =  B[k1 * n + j1];
                            C_0_0 += A_0_M * B_M;
                            C_1_0 += A_1_M * B_M;
                            C_2_0 += A_2_M * B_M;

                            B_M = B[k1 * n + (j1 + 1)];
                            C_0_1 += A_0_M * B_M;
                            C_1_1 += A_1_M * B_M;
                            C_2_1 += A_2_M * B_M;

                            B_M = B[k1 * n + (j1 + 2)];
                            C_0_2 += A_0_M * B_M;
                            C_1_2 += A_1_M * B_M;
                            C_2_2 += A_2_M * B_M;

                        }
                        C[i1 * n + j1] = C_0_0;
                        C[(i1 + 1) * n + j1] = C_1_0;
                        C[(i1 + 2) * n + j1] = C_2_0;

                        C[i1 * n + (j1 + 1)] = C_0_1;
                        C[(i1 + 1) * n + (j1 + 1)] = C_1_1;
                        C[(i1 + 2) * n + (j1 + 1)] = C_2_1;

                        C[i1 * n + (j1 + 2)] = C_0_2;
                        C[(i1 + 1) * n + (j1 + 2)] = C_1_2;
                        C[(i1 + 2) * n + (j1 + 2)] = C_2_2;
                    
                    }
                }
            }
        }
    }*/
}