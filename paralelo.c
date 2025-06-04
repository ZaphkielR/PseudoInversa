#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

// Función para imprimir una matriz en la consola
void printMatrix(double** matriz, int m, int n) {
    #pragma omp parallel for
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            printf("%lf ", matriz[i][j]);
        }
        printf("\n");
    }
}

// Función para leer una matriz desde un archivo
double** leer_matriz(const char* nombre_archivo, int* m, int* n) {
    FILE* archivo = fopen(nombre_archivo, "r");
    if (!archivo) {
        perror("Error al abrir el archivo");
        return NULL;
    }

    if (fscanf(archivo, "%d %d", m, n) != 2) {
        fprintf(stderr, "Error al leer las dimensiones del archivo\n");
        fclose(archivo);
        return NULL;
    }

    double** matriz = (double**)malloc(*m * sizeof(double*));
    if (!matriz) {
        fprintf(stderr, "Error al asignar memoria para la matriz\n");
        fclose(archivo);
        return NULL;
    }

    #pragma omp parallel for
    for (int i = 0; i < *m; i++) {
        matriz[i] = (double*)malloc(*n * sizeof(double));
        if (!matriz[i]) {
            fprintf(stderr, "Error al asignar memoria para la fila %d\n", i);
            for (int j = 0; j < i; j++)
                free(matriz[j]);
            free(matriz);
            fclose(archivo);
            return NULL;
        }
    }

    #pragma omp parallel for
    for (int i = 0; i < *m; i++) {
        for (int j = 0; j < *n; j++) {
            if (fscanf(archivo, "%lf", &matriz[i][j]) != 1) {
                fprintf(stderr, "Error al leer el elemento (%d, %d)\n", i, j);
                for (int k = 0; k <= i; k++)
                    free(matriz[k]);
                free(matriz);
                fclose(archivo);
                return NULL;
            }
        }
    }

    fclose(archivo);
    return matriz;
}

// Función para calcular el rango de una matriz usando eliminación gaussiana
int rango(int m, int n, double** A) {
    double temp[m][n];
    int rank = 0;

    #pragma omp parallel for
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            temp[i][j] = A[i][j];
        }
    }

    for (int i = 0; i < n; i++) {
        int pivotRow = -1;
        double maxVal = 0;

        #pragma omp parallel for reduction(max:maxVal)
        for (int j = rank; j < m; j++) {
            double absVal = fabs(temp[j][i]);
            if (absVal > 1e-10 && absVal > maxVal) {
                #pragma omp critical
                {
                    if (absVal > maxVal) {
                        maxVal = absVal;
                        pivotRow = j;
                    }
                }
            }
        }

        if (pivotRow != -1) {
            #pragma omp parallel for
            for (int k = 0; k < n; k++) {
                double tmp = temp[rank][k];
                temp[rank][k] = temp[pivotRow][k];
                temp[pivotRow][k] = tmp;
            }

            #pragma omp parallel for
            for (int j = 0; j < m; j++) {
                if (j != rank) {
                    double factor = temp[j][i] / temp[rank][i];
                    #pragma omp parallel for
                    for (int k = 0; k < n; k++)
                        temp[j][k] -= factor * temp[rank][k];
                }
            }

            rank++;
        }
    }

    return rank;
}

// Función para calcular la transpuesta de una matriz
double** transpuesta(int m, int n, double** A) {
    double** B = (double**)malloc(n * sizeof(double*));

    #pragma omp parallel for
    for (int i = 0; i < n; i++) {
        B[i] = (double*)malloc(m * sizeof(double));
    }

    #pragma omp parallel for
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            B[j][i] = A[i][j];
        }
    }

    return B;
}

// Función para multiplicar dos matrices
double** multiplicar(int m, int n, int p, double** A, double** B) {
    double** C = (double**)malloc(m * sizeof(double*));
    
    #pragma omp parallel for
    for (int i = 0; i < m; i++) {
        C[i] = (double*)malloc(p * sizeof(double));
    }

    #pragma omp parallel for collapse(2)
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < p; j++) {
            double suma = 0;
            for (int k = 0; k < n; k++) {
                suma += A[i][k] * B[k][j];
            }
            C[i][j] = suma;
        }
    }

    return C;
}

// Función para calcular la inversa de una matriz cuadrada usando Gauss-Jordan
double** inverse(int m, double** A) {
    const double EPS = 1e-12;
    double** aug = (double**)malloc(m * sizeof(double*));

    #pragma omp parallel for
    for (int i = 0; i < m; i++) {
        aug[i] = (double*)malloc(2 * m * sizeof(double));
        if (!aug[i]) {
            fprintf(stderr, "Error al reservar memoria (fila aug %d)\n", i);
            for (int k = 0; k < i; k++) free(aug[k]);
            free(aug);
            return NULL;
        }

        #pragma omp parallel for
        for (int j = 0; j < m; j++) {
            aug[i][j] = A[i][j];
            aug[i][j + m] = (i == j) ? 1.0 : 0.0;
        }
    }

    for (int col = 0, row = 0; col < m && row < m; col++, row++) {
        int piv = row;
        double maxVal = fabs(aug[row][col]);

        #pragma omp parallel for reduction(max:maxVal)
        for (int k = row + 1; k < m; k++) {
            double absVal = fabs(aug[k][col]);
            if (absVal > maxVal) {
                #pragma omp critical
                {
                    if (absVal > maxVal) {
                        maxVal = absVal;
                        piv = k;
                    }
                }
            }
        }

        if (fabs(aug[piv][col]) < EPS) {
            for (int i = 0; i < m; i++) free(aug[i]);
            free(aug);
            fprintf(stderr, "Error: matriz no invertible (pivote cero)\n");
            return NULL;
        }

        if (piv != row) {
            #pragma omp parallel for
            for (int k = 0; k < 2 * m; k++) {
                double tmp = aug[piv][k];
                aug[piv][k] = aug[row][k];
                aug[row][k] = tmp;
            }
        }

        double piv_val = aug[row][col];
        #pragma omp parallel for
        for (int j = 0; j < 2 * m; j++)
            aug[row][j] /= piv_val;

        #pragma omp parallel for
        for (int i = 0; i < m; i++) {
            if (i == row) continue;
            double factor = aug[i][col];
            #pragma omp parallel for
            for (int j = 0; j < 2 * m; j++)
                aug[i][j] -= factor * aug[row][j];
        }
    }

    double** inv = (double**)malloc(m * sizeof(double*));
    #pragma omp parallel for
    for (int i = 0; i < m; i++) {
        inv[i] = (double*)malloc(m * sizeof(double));
        if (!inv[i]) {
            fprintf(stderr, "Error al reservar memoria (fila inv %d)\n", i);
            for (int k = 0; k < i; k++) free(inv[k]);
            free(inv);
            for (int k = 0; k < m; k++) free(aug[k]);
            free(aug);
            return NULL;
        }
        #pragma omp parallel for
        for (int j = 0; j < m; j++)
            inv[i][j] = aug[i][j + m];
    }

    #pragma omp parallel for
    for (int i = 0; i < m; i++)
        free(aug[i]);
    free(aug);

    return inv;
}

int main() {
    int m, n;
    double** matriz = leer_matriz("entrada.ent", &m, &n);

    if (!matriz) {
        printf("No se pudo leer la matriz.\n");
        return 1;
    }

    printMatrix(matriz, m, n);
    int r = rango(m, n, matriz);

    if (r < (m < n ? m : n)) {
        printf("-1\n");
        return 0;
    }

    if (r == m) {
        double** AT = transpuesta(m, n, matriz);
        printMatrix(AT, n, m);

        double** ATA = multiplicar(m, n, m, matriz, AT);
        printMatrix(ATA, m, m);

        double** invATA = inverse(m, ATA);
        printMatrix(invATA, m, m);
        
        double** P = multiplicar(n, m, m, AT, invATA);
        printMatrix(P, n, m);
    }

    if (r == n) {
        double** AT = transpuesta(m, n, matriz);
        printMatrix(AT, n, m);

        double** ATA = multiplicar(n, m, n, AT, matriz);
        printMatrix(ATA, n, n);

        double** invATA = inverse(n, ATA);
        printMatrix(invATA, n, n);
        
        double** P = multiplicar(n, n, m, invATA, AT);
        printMatrix(P, n, m);
    }

    return 0;
}

