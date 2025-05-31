#include <stdio.h>
#include <stdlib.h>

// Definición de constantes para las dimensiones de la matriz
// ROWS: número de filas (3)
// COLS: número de columnas (2)
#define ROWS 3
#define COLS 2

// Matriz original A de 3x2 con valores constantes
// Esta matriz será la base para calcular su pseudo-inversa
const double A[ROWS][COLS] = {
    {1, 2},  // Fila 1
    {3, 4},  // Fila 2
    {5, 6}   // Fila 3
};


// Función para imprimir cualquier matriz en la consola
// Parámetros:
// - rows: número de filas de la matriz
// - cols: número de columnas de la matriz
// - A: matriz a imprimir (constante para evitar modificaciones)
void printMatrix(int rows, int cols, const double A[rows][cols]) {
    // Bucle externo para filas
    for (int i = 0; i < rows; i++) {
        // Bucle interno para columnas
        for (int j = 0; j < cols; j++) {
            // Imprime cada elemento con formato: 8 dígitos totales y 4 decimales
            printf("%8.4f ", A[i][j]);
        }
        // Nueva línea después de cada fila
        printf("\n");
    }
}

// Función para transponer una matriz
// Parámetros:
// - rows: número de filas de la matriz original
// - cols: número de columnas de la matriz original
// - A: matriz original a transponer
// - B: matriz resultado donde se guardará la transpuesta
void transpose(int rows, int cols, const double A[rows][cols], double B[cols][rows]) {
    // Bucle para recorrer cada elemento de la matriz
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            // Intercambia filas por columnas
            B[j][i] = A[i][j];
        }
    }
}

// Función para multiplicar dos matrices
// Parámetros:
// - m: número de filas de A
// - n: número de columnas de A y filas de B
// - p: número de columnas de B
// - A: primera matriz
// - B: segunda matriz
// - C: matriz resultado de la multiplicación
void multiplyMatrix(int m, int n, int p, const double A[m][n], const double B[n][p], double C[m][p]) {
    // Bucle externo para filas de C
    for (int i = 0; i < m; i++) {
        // Bucle interno para columnas de C
        for (int j = 0; j < p; j++) {
            // Inicializa cada elemento en 0
            C[i][j] = 0;
            // Bucle para sumar productos de elementos correspondientes
            for (int k = 0; k < n; k++) {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
}

// Función para calcular la inversa de una matriz cuadrada usando Gauss-Jordan
// Parámetros:
// - n: dimensión de la matriz cuadrada
// - A: matriz original
// - inv: matriz donde se guardará la inversa
void inverseMatrix(int n, double A[n][n], double inv[n][n]) {
    // Crear matriz aumentada de tamaño n x 2n
    double aug[n][2 * n];

    // Construir matriz aumentada [A | I]
    // Donde I es la matriz identidad
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            // Copiar elementos de A a la parte izquierda
            aug[i][j] = A[i][j];
            // Crear matriz identidad en la parte derecha
            aug[i][j + n] = (i == j) ? 1.0 : 0.0;
        }
    }

    // Método de Gauss-Jordan para encontrar la inversa
    for (int i = 0; i < n; i++) {
        // Obtener el pivote (elemento diagonal)
        double pivot = aug[i][i];

        // Normalizar la fila i dividiendo por el pivote
        for (int j = 0; j < 2 * n; j++) {
            aug[i][j] /= pivot;
        }

        // Eliminar elementos en otras filas usando la fila normalizada
        for (int k = 0; k < n; k++) {
            if (k != i) {
                // Calcular factor para eliminar elemento
                double factor = aug[k][i];
                for (int j = 0; j < 2 * n; j++) {
                    // Eliminar elemento usando la fila normalizada
                    aug[k][j] -= factor * aug[i][j];
                }
            }
        }
    }

    // Extraer la matriz inversa de la parte derecha de la matriz aumentada
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            // La inversa está en la parte derecha de la matriz aumentada
            inv[i][j] = aug[i][j + n];
        }
    }
}

int main() {
    // Imprimir la matriz original A
    printf("Matriz original A:\n");
    printMatrix(ROWS, COLS, A);
    printf("\n");

    // Calcular y mostrar la transpuesta de A
    printf("Transpuesta de A (AT):\n");
    double AT[COLS][ROWS];
    transpose(ROWS, COLS, A, AT);
    printMatrix(COLS, ROWS, AT);
    printf("\n");

    // Calcular y mostrar el producto AT * A
    printf("Producto AT * A:\n");
    double ATA[COLS][COLS];
    multiplyMatrix(COLS, ROWS, COLS, AT, A, ATA);
    printMatrix(COLS, COLS, ATA);
    printf("\n");

    // Calcular y mostrar la inversa de ATA
    printf("Inversa de ATA:\n");
    double invATA[COLS][COLS];
    inverseMatrix(COLS, ATA, invATA);
    printMatrix(COLS, COLS, invATA);
    printf("\n");

    // Calcular y mostrar la pseudo-inversa final P
    printf("Pseudo-inversa de A:\n");
    double P[COLS][ROWS];
    multiplyMatrix(COLS, COLS, ROWS, invATA, AT, P);
    printMatrix(COLS, ROWS, P);
    
    return 0;
}