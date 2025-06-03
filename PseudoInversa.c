#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Definición de constante para el tamaño máximo de la matriz
// MAX: número máximo de filas o columnas que puede tener la matriz
#define MAX 100

// -----------------------------------------------------------------------------
// Función para imprimir cualquier matriz en la consola
// Parámetros:
// - rows: número de filas de la matriz
// - cols: número de columnas de la matriz
// - A: matriz a imprimir (constante para evitar modificaciones)
// -----------------------------------------------------------------------------
void printMatrix(int rows, int cols, const double A[rows][cols]) {
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            printf("%10.6f ", A[i][j]);
        }
        printf("\n");
    }
}

// -----------------------------------------------------------------------------
// Función para transponer una matriz
// Parámetros:
// - rows: número de filas de la matriz original
// - cols: número de columnas de la matriz original
// - A: matriz original a transponer
// - B: matriz resultado donde se guardará la transpuesta
// -----------------------------------------------------------------------------
void transpose(int rows, int cols, const double A[rows][cols], double B[cols][rows]) {
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            B[j][i] = A[i][j];
}

// -----------------------------------------------------------------------------
// Función para multiplicar dos matrices
// Parámetros:
// - m: número de filas de A
// - n: número de columnas de A y filas de B
// - p: número de columnas de B
// - A: primera matriz
// - B: segunda matriz
// - C: matriz resultado de la multiplicación
// -----------------------------------------------------------------------------
void multiplyMatrix(int m, int n, int p, const double A[m][n], const double B[n][p], double C[m][p]) {
    for (int i = 0; i < m; i++)
        for (int j = 0; j < p; j++) {
            C[i][j] = 0;
            for (int k = 0; k < n; k++)
                C[i][j] += A[i][k] * B[k][j];
        }
}

// -----------------------------------------------------------------------------
// Función para calcular la inversa de una matriz cuadrada usando Gauss-Jordan
// Parámetros:
// - n: dimensión de la matriz cuadrada
// - A: matriz original
// - inv: matriz donde se guardará la inversa
// -----------------------------------------------------------------------------
void inverseMatrix(int n, double A[n][n], double inv[n][n]) {
    double aug[n][2 * n];  // Matriz aumentada [A | I]

    // Construcción de la matriz aumentada
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) {
            aug[i][j] = A[i][j];  // Parte izquierda A
            aug[i][j + n] = (i == j) ? 1.0 : 0.0;  // Parte derecha I
        }

    // Aplicar el método de Gauss-Jordan
    for (int i = 0; i < n; i++) {
        double pivot = aug[i][i];

        // Validación de pivote distinto de cero
        if (fabs(pivot) < 1e-10) {
            fprintf(stderr, "Error: matriz no invertible (pivote cero)\n");
            exit(1);
        }

        // Normalizar fila pivote
        for (int j = 0; j < 2 * n; j++)
            aug[i][j] /= pivot;

        // Eliminar elementos en otras filas
        for (int k = 0; k < n; k++) {
            if (k != i) {
                double factor = aug[k][i];
                for (int j = 0; j < 2 * n; j++)
                    aug[k][j] -= factor * aug[i][j];
            }
        }
    }

    // Extraer la matriz inversa desde la parte derecha de la matriz aumentada
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            inv[i][j] = aug[i][j + n];
}

// -----------------------------------------------------------------------------
// Función para calcular el rango de una matriz usando eliminación de Gauss
// Parámetros:
// - m: número de filas
// - n: número de columnas
// - A: matriz original
// Retorna:
// - rango de la matriz (número de filas no nulas en forma escalonada)
// -----------------------------------------------------------------------------
int rango(int m, int n, double A[m][n]) {
    double temp[m][n];
    int rank = 0;

    // Copiar la matriz original a una temporal
    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
            temp[i][j] = A[i][j];

    // Recorrer columnas buscando pivotes
    for (int i = 0; i < n; i++) {
        int pivotRow = -1;

        // Buscar fila con elemento distinto de 0 en la columna i
        for (int j = rank; j < m; j++) {
            if (fabs(temp[j][i]) > 1e-10) {
                pivotRow = j;
                break;
            }
        }

        if (pivotRow != -1) {
            // Intercambiar filas
            for (int k = 0; k < n; k++) {
                double tmp = temp[rank][k];
                temp[rank][k] = temp[pivotRow][k];
                temp[pivotRow][k] = tmp;
            }

            // Eliminar elementos de otras filas
            for (int j = 0; j < m; j++) {
                if (j != rank) {
                    double factor = temp[j][i] / temp[rank][i];
                    for (int k = 0; k < n; k++)
                        temp[j][k] -= factor * temp[rank][k];
                }
            }

            rank++;  // Aumentar el conteo de filas linealmente independientes
        }
    }

    return rank;
}

// -----------------------------------------------------------------------------
// Función principal del programa
// Lee la matriz desde "entrada.ent", calcula la pseudoinversa si es posible,
// y guarda el resultado en "salida.sal".
// -----------------------------------------------------------------------------
int main() {
    int m, n;
    double A[MAX][MAX];

    // Abrir archivo de entrada
    FILE *fin = fopen("entrada.ent", "r");
    if (!fin) {
        perror("Error al abrir entrada.ent");
        return 1;
    }

    // Leer dimensiones de la matriz
    fscanf(fin, "%d %d", &m, &n);

    // Leer elementos de la matriz A
    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
            fscanf(fin, "%lf", &A[i][j]);

    fclose(fin);

    // Calcular el rango de la matriz
    int r = rango(m, n, A);

    // Abrir archivo de salida
    FILE *fout = fopen("salida.sal", "w");
    if (!fout) {
        perror("Error al crear salida.sal");
        return 1;
    }

    // Verificar si tiene rango completo
    if (r < (m < n ? m : n)) {
        // No tiene pseudoinversa
        fprintf(fout, "-1\n");
        fclose(fout);
        return 0;
    }

    // Caso: pseudoinversa por la izquierda
    if (r == n) {
        double AT[n][m];
        transpose(m, n, A, AT);

        double ATA[n][n];
        multiplyMatrix(n, m, n, AT, A, ATA);

        double invATA[n][n];
        inverseMatrix(n, ATA, invATA);

        double P[n][m];
        multiplyMatrix(n, n, m, invATA, AT, P);

        fprintf(fout, "L\n");
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                fprintf(fout, "%.6f ", P[i][j]);
            }
            fprintf(fout, "\n");
        }

    // Caso: pseudoinversa por la derecha
    } else if (r == m) {
        double AT[n][m];
        transpose(m, n, A, AT);

        double AAT[m][m];
        multiplyMatrix(m, n, m, A, AT, AAT);

        double invAAT[m][m];
        inverseMatrix(m, AAT, invAAT);

        double P[n][m];
        multiplyMatrix(n, m, m, AT, invAAT, P);

        fprintf(fout, "R\n");
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                fprintf(fout, "%.6f ", P[i][j]);
            }
            fprintf(fout, "\n");
        }

    } else {
        // No tiene pseudoinversa válida
        fprintf(fout, "-1\n");
    }

    fclose(fout);
    return 0;
}
