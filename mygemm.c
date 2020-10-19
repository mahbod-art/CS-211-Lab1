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
    int i = 0;
    int j = 0;
    int k = 0;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            for (k = 0; k < n; k++)
            {
                C[i * n + j] += A[i * n + k] * B[k * n + j];
            }
        }
    }
}

void dgemm1(const double* A, const double* B, double* C, const int n)
{
    int i = 0;
    int j = 0;
    int k = 0;
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
    }
}
//Register Reuse part 1 End

//Register Reuse part 2
void dgemm2(const double* A, const double* B, double* C, const int n)
{
    int i = 0;
    int j = 0;
    int k = 0;
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

                C1 += A1 * B1 + A3 * B2;
                C2 += A2 * B1 + A4 * B2;
                C3 += A1 * B3 + A3 * B4;
                C4 += A2 * B3 + A4 * B4;
            }

            C[i * n + j] = C1;
            if (i < (n - 1)) C[(i + 1) * n + j] = C2;
            if (j < (n - 1)) C[i * n + (j + 1)] = C3;
            if (i < (n - 1) && j < (n - 1)) C[(i + 1) * n + (j + 1)] = C4;
        }
    }
}

//Register Reuse part 2 End

//Register Reuse part 3
void dgemm3(const double* A, const double* B, double* C, const int n)
{
    int i = 0;
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
    }
}
//Register Reuse part 3 End

//Cache Reuse part 3
void ijk(const double *A, const double *B, double *C, const int n) 
{
    int i = 0;
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
    }
}

void bijk(const double *A, const double *B, double *C, const int n, const int b) 
{
    int i = 0;
    int j = 0;
    int k = 0;
    int l = 0;
    int m = 0;
    for (i = 0; i < n; i += b)
    {
        for (j = 0; j < n; j += b)
        {
            for (k = 0; k < n; k += b)
            {
                for (l = i; l < i + b && l < n; l++)
                {
                    for (m = j; m < j + b && m < n; m++)
                    {
                        register double Clm = C[l * n + m];
                        int n = 0;
                        for (n = k; n < k + b && n < n; n++)
                        {
                            Clm += A[l * n + n] * B[n * n + m];
                        }
                        C[l * n + m] = Clm;
                    }
                }
            }
        }
    }
}

void jik(const double *A, const double *B, double *C, const int n) 
{
    int i = 0;
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
    }
}

void bjik(const double *A, const double *B, double *C, const int n, const int b) 
{
    int i = 0;
    int j = 0;
    int k = 0;
    int l = 0;
    int m = 0;
    for (j = 0; j < n; j += b)
    {
        for (i = 0; i < n; i += b)
        {
            for (k = 0; k < n; k += b)
            {
                for (m = j; m < j + b && m < n; m++)
                {
                    for (l = i; l < i + b && l < n; l++)
                    {
                        register double Clm = C[l * n + m];
                        int n = 0;
                        for (n = k; n < k + b && n < n; n++)
                        {
                            Clm += A[l * n + n] * B[n * n + m];
                        }
                        C[l * n + m] = Clm;
                    }
                }
            }
        }
    }
}

void kij(const double *A, const double *B, double *C, const int n) 
{
    int i = 0;
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
    }
}

void bkij(const double *A, const double *B, double *C, const int n, const int b) 
{
    int i = 0;
    int j = 0;
    int k = 0;
    int l = 0;
    int m = 0;
    for (j = 0; j < n; j += b)
    {
        for (i = 0; i < n; i += b)
        {
            for (k = 0; k < n; k += b)
            {
                for (m = j; m < j + b && m < n; m++)
                {
                    for (l = i; l < i + b && l < n; l++)
                    {
                        register double Clm = C[l * n + m];
                        int n = 0;
                        for (n = k; n < k + b && n < n; n++)
                        {
                            Clm += A[l * n + n] * B[n * n + m];
                        }
                        C[l * n + m] = Clm;
                    }
                }
            }
        }
    }
}


void ikj(const double *A, const double *B, double *C, const int n) 
{
    int i = 0;
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
    }
}

void bikj(const double *A, const double *B, double *C, const int n, const int b) 
{
    int i = 0;
    int j = 0;
    int k = 0;
    int l = 0;
    int m = 0;
    for (k = 0; k < n; k += b)
    {
        for (i = 0; i < n; i += b)
        {
            for (j = 0; j < n; j += b)
            {
                int n = 0;
                for (n = k; n < k + b && n < n; n++)
                {
                    for (l = i; l < i + b && l < n; l++)
                    {
                        register double Aln = A[l * n + n];
                        for (m = j; m < j + b && m < n; m++)
                        {
                            C[l * n + m] += Aln * B[n * n + m];
                        }
                    }
                }
            }
        }
    }
}

void jki(const double *A, const double *B, double *C, const int n) 
{
    int i = 0;
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
    }
}

void bjki(const double *A, const double *B, double *C, const int n, const int b) 
{
    int i = 0;
    int j = 0;
    int k = 0;
    int l = 0;
    int m = 0;
    for (i = 0; i < n; i += b)
    {
        for (k = 0; k < n; k += b)
        {
            for (j = 0; j < n; j += b)
            {
                for (l = i; l < i + b && l < n; l++)
                {
                    int n = 0;
                    for (n = k; n < k + b && n < n; n++)
                    {
                        register double Aln = A[l * n + n];
                        for (m = j; m < j + b && m < n; m++)
                        {
                            C[l * n + m] += Aln * B[n * n + m];
                        }
                    }
                }
            }
        }
    }
}

void kji(const double *A, const double *B, double *C, const int n) 
{
    int i = 0;
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
    }
}

void bkji(const double *A, const double *B, double *C, const int n, const int b) 
{
    int k = 0;
    for (k = 0; k < n; k += b)
    {
        int j = 0;
        for (j = 0; j < n; j += b)
        {
            int i = 0;
            for (i = 0; i < n; i += b)
            {
                int k1 = 0;
                for (k1 = k; k1 < k + b && k1 < n; k1++)
                {
                    int j1 = 0;
                    for (j1 = j; j1 < j + b && j1 < n; j1++)
                    {
                        register double B_k1_j1 = B[k1 * n + j1];
                        int i1 = 0;
                        for (i1 = i; i1 < i + b && i1 < n; i1++)
                        {
                            C[i1 * n + j1] += A[i1 * n + k1] * B_k1_j1;
                        }
                    }
                }
            }
        }
    }
}
//Cache Reuse part 3 End 

//Cache Reuse part 4
void optimal(const double* A, const double* B, double *C, const int n, const int b)
{

}