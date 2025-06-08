//Alumnos: Christian Delgado, Rafael Morales
//Fecha: 10-06-2025
//Profesor: Sergio Antonio Baltierra Valenzuela
//Carrera: Ingenieria Civil Informatica
//Ramo: Computación de Alto Rendimiento

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Función para imprimir una matriz en la consola
void printMatrix(double** matriz, char* label, int m, int n) {    
    FILE* archivo = fopen("salida_sec.sal", "w");

    if (!archivo) {
        perror("Error al abrir el archivo");
        return;
    }

    fprintf(archivo, "%s\n", label);
    
    // Imprime cada elemento de la matriz en formato de tabla
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            fprintf(archivo, "%lf ", matriz[i][j]);
        }
        fprintf(archivo, "\n");
    } 

    fclose(archivo);
}

// Función para leer una matriz desde un archivo
// Retorna NULL si hay algún error durante la lectura
// Parámetros:
// - nombre_archivo: nombre del archivo a leer
// - m: número de filas (salida)
// - n: número de columnas (salida)
double** leer_matriz(const char* nombre_archivo, int* m, int* n) {
    FILE* archivo = fopen(nombre_archivo, "r");
    if (!archivo) {
        perror("Error al abrir el archivo");
        return NULL;
    }

    // Leer dimensiones de la matriz (primera línea del archivo)
    if (fscanf(archivo, "%d %d", m, n) != 2) {
        fprintf(stderr, "Error al leer las dimensiones del archivo\n");
        fclose(archivo);
        return NULL;
    }

    // Asignar memoria para la matriz (matriz de punteros)
    double** matriz = (double**)malloc(*m * sizeof(double*));
    if (!matriz) {
        fprintf(stderr, "Error al asignar memoria para la matriz\n");
        fclose(archivo);
        return NULL;
    }

    // Asignar memoria para cada fila
    for (int i = 0; i < *m; i++) {
        matriz[i] = (double*)malloc(*n * sizeof(double));
        if (!matriz[i]) {
            fprintf(stderr, "Error al asignar memoria para la fila %d\n", i);
            fclose(archivo);
            return NULL;
        }
    }

    // Leer elementos de la matriz del archivo
    for (int i = 0; i < *m; i++) {
        for (int j = 0; j < *n; j++) {
            if (fscanf(archivo, "%lf", &matriz[i][j]) != 1) {
                fprintf(stderr, "Error al leer el elemento (%d, %d)\n", i, j);
                fclose(archivo);
                return NULL;
            }
        }
    }

    fclose(archivo);
    return matriz;
}

// Función para calcular el rango de una matriz usando eliminación gaussiana
// Retorna el rango de la matriz
double rango(int m, int n, double** A) {
    // Crear una copia temporal de la matriz para no modificar la original
    double temp[m][n];
    int rank = 0;

    // Copiar la matriz original a la temporal
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            temp[i][j] = A[i][j];
        }
    }

    // Algoritmo de eliminación gaussiana
    for (int i = 0; i < n; i++) {
        int pivotRow = -1;

        // Buscar el pivote (primer elemento no nulo)
        for (int j = rank; j < m; j++) {
            if (fabs(temp[j][i]) > 1e-10) {
                pivotRow = j;
                break;
            }
        }

        // Si encontramos un pivote
        if (pivotRow != -1) {
            // Intercambiar filas si es necesario
            for (int k = 0; k < n; k++) {
                int tmp = temp[rank][k];
                temp[rank][k] = temp[pivotRow][k];
                temp[pivotRow][k] = tmp;
            }

            // Eliminación gaussiana: hacer cero los elementos debajo del pivote
            for (int j = 0; j < m; j++) {
                if (j != rank) {
                    int factor = temp[j][i] / temp[rank][i];
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
    // Crear matriz transpuesta (n x m)
    double** B = (double**)malloc(n * sizeof(double*));

    // Asignar memoria para cada fila
    for (int i = 0; i < n; i++) {
        B[i] = (double*)malloc(m * sizeof(double));
    }

    // Calcular la transpuesta: B[i][j] = A[j][i]
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            B[j][i] = A[i][j];
        }
    }

    return B;
}

// Función para multiplicar dos matrices
double** multiplicar(int m, int n, int p, double** A, double** B) {
    // Crear matriz resultado (m x p)
    double** C = (double**)malloc(m * sizeof(double*));
    
    // Asignar memoria para cada fila
    for (int i = 0; i < m; i++) {
        C[i] = (double*)malloc(p * sizeof(double));
    }

    // Realizar la multiplicación de matrices
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < p; j++) {
            C[i][j] = 0;
            // Sumatorio de la multiplicación
            for (int k = 0; k < n; k++) {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }

    return C;
}

// Función para calcular la inversa de una matriz cuadrada usando Gauss-Jordan
double** inverse(int m, double** A) {
    const double EPS = 1e-12;  // Tolerancia para considerar un número como cero
    
    // Crear matriz aumentada [A | I] donde I es la matriz identidad
    double** aug = (double**)malloc(m * sizeof(double*));
    if (!aug) {
        fprintf(stderr, "Error al reservar memoria (aug)\n");
        return NULL;
    }

    // Asignar memoria para cada fila
    for (int i = 0; i < m; i++) {
        aug[i] = (double*)malloc(2 * m * sizeof(double));
        if (!aug[i]) {
            fprintf(stderr, "Error al reservar memoria (fila aug %d)\n", i);
            return NULL;
        }

        // Llenar la matriz aumentada
        // Primera mitad: copia de A
        for (int j = 0; j < m; j++) {         
            aug[i][j] = A[i][j];
        }
        // Segunda mitad: matriz identidad
        for (int j = 0; j < m; j++) {         
            aug[i][j + m] = (i == j) ? 1.0 : 0.0;
        }
    }

    /* Algoritmo de Gauss-Jordan para calcular la inversa */
    for (int col = 0, row = 0; col < m && row < m; col++, row++) {

        /* 1. Buscar pivote (elemento máximo en la columna) */
        int piv = row;
        for (int k = row + 1; k < m; k++)
            if (fabs(aug[k][col]) > fabs(aug[piv][col]))
                piv = k;

        // Si no se encuentra un pivote válido, la matriz no es invertible
        if (fabs(aug[piv][col]) < EPS) {
            fprintf(stderr, "Error: matriz no invertible (pivote cero)\n");
            return NULL;
        }

        /* 2. Intercambiar filas si es necesario para poner el pivote en la posición correcta */
        if (piv != row) {
            double* tmp = aug[piv];
            aug[piv] = aug[row];
            aug[row] = tmp;
        }

        /* 3. Normalizar la fila pivote (dividir por el pivote) */
        double piv_val = aug[row][col];
        for (int j = 0; j < 2 * m; j++)
            aug[row][j] /= piv_val;

        /* 4. Eliminar la columna en las demás filas (hacer ceros debajo del pivote) */
        for (int i = 0; i < m; i++) {
            if (i == row) continue;
            double factor = aug[i][col];
            for (int j = 0; j < 2 * m; j++)
                aug[i][j] -= factor * aug[row][j];
        }
    }

    /* Extraer la parte derecha como matriz inversa */
    double** inv = (double**)malloc(m * sizeof(double*));
    if (!inv) {
        fprintf(stderr, "Error al reservar memoria (inv)\n");
        return NULL;
    }

    // Copiar la parte derecha de la matriz aumentada (que ahora contiene la inversa)
    for (int i = 0; i < m; i++) {
        inv[i] = (double*)malloc(m * sizeof(double));
        if (!inv[i]) {
            fprintf(stderr, "Error al reservar memoria (fila inv %d)\n", i);
            return NULL;
        }
        // Copiar la parte derecha de la matriz aumentada
        for (int j = 0; j < m; j++)
            inv[i][j] = aug[i][j + m];
    }

    return inv;
}

int main(int argc, char* argv[]) {
    // Variables para las dimensiones de la matriz
    int m, n;

    // Leer la matriz del archivo de entrada entregado por terminal
    const char *entrada = argv[1];  // Guardamos la cadena en una variable
    double** matriz = leer_matriz(entrada, &m, &n);

    // Verificar si hubo error al leer la matriz
    if (!matriz) {
        printf("No se pudo leer la matriz.\n");
        return 1;
    }

    // Calcular el rango de la matriz
    int r = rango(m, n, matriz);

    // Si el rango es menor que el mínimo entre filas y columnas,
    // la matriz no tiene pseudo-inversa
    if (r < (m < n ? m : n)) {
        printf("-1\n");
        return 0;
    }

    // Caso 1: rango = número de filas (m)
    if (r == m) {
        // Calcular la transpuesta de A
        double** AT = transpuesta(m, n, matriz);

        // Calcular A * A^T
        double** ATA = multiplicar(m, n, m, matriz, AT);

        // Calcular la inversa de A * A^T
        double** invATA = inverse(m, ATA);
        
        // Calcular la pseudo-inversa: P = A^T * (A * A^T)^-1
        double** P = multiplicar(n, m, m, AT, invATA);

        printMatrix(P, "R", n, m);
    }

    // Caso 2: rango = número de columnas (n)
    if (r == n) {
        // Calcular la transpuesta de A
        double** AT = transpuesta(m, n, matriz);

        // Calcular A^T * A
        double** ATA = multiplicar(n, m, n, AT, matriz);

        // Calcular la inversa de A^T * A
        double** invATA = inverse(n, ATA);
        
        // Calcular la pseudo-inversa: P = (A^T * A)^-1 * A^T
        double** P = multiplicar(n, n, m, invATA, AT);

        printMatrix(P, "L", n, m);
    }

    return 0;
}