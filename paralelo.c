/*
 * Programa para calcular la pseudo-inversa de una matriz usando paralelización OpenMP
 * Autores: Christian Delgado, Rafael Morales
 * Fecha: 10-06-2025
 * Ramo: Computación de Alto Rendimiento
 * Docente: Sergio Antonio Baltierra Valenzuela
 * 
 * Compilar: gcc -fopenmp paralelo.c -o paralelo -lm [entrada] [num_threads]
 */
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

/**
 * Función que imprime una matriz en un archivo de salida
 * @param label Etiqueta para identificar la matriz en el archivo
 * @param m Número de filas de la matriz
 * @param n Número de columnas de la matriz
 * @param matriz Matriz a imprimir
 */
void printMatrix(char* label, int m, int n, double matriz[m][n]) {    
    FILE* archivo = fopen("salida.sal", "w");

    if (!archivo) {
        perror("Error al abrir el archivo");
        return;
    }

    fprintf(archivo, "%s\n", label);
    
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            fprintf(archivo, "%.6f ", matriz[i][j]);
        }
        fprintf(archivo, "\n");
    } 

    fclose(archivo);
}



/**
 * Función que calcula el rango de una matriz usando eliminación gaussiana
 * @param m Número de filas de la matriz
 * @param n Número de columnas de la matriz
 * @param A Matriz de entrada
 * @return El rango de la matriz
 * 
 * La función crea una copia de la matriz y aplica eliminación gaussiana
 * para encontrar el rango. Utiliza paralelización OpenMP para la operación
 * de reducción de filas.
 */
int rango(int m, int n, double A[m][n]) {
    const double EPS = 1e-16;

    // Crear copia de la matriz original
    double temp[m][n];
    int rank = 0;

    // Copiar matriz original a matriz temporal
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            temp[i][j] = A[i][j];
        }
    }

    // Para cada columna
    for (int i = 0; i < n; i++) {
        int pivotRow = -1;

        // Buscar pivote en la columna actual
        for (int j = rank; j < m; j++) {
            if (fabs(temp[j][i]) > EPS) {
                pivotRow = j;
                break;
            }
        }

        // Si se encontró pivote
        if (pivotRow != -1) {
            // Intercambiar fila pivote con fila actual
            for (int k = 0; k < n; k++) {
                int tmp = temp[rank][k];
                temp[rank][k] = temp[pivotRow][k];
                temp[pivotRow][k] = tmp;
            }

            // Reducir las demás filas usando la fila pivote (paralelo)
            #pragma omp parallel for
            for (int j = 0; j < m; j++) {
                if (j != rank) {
                    int factor = temp[j][i] / temp[rank][i];
                    for (int k = 0; k < n; k++)
                        temp[j][k] -= factor * temp[rank][k];
                }
            }

            // Incrementar rango encontrado
            rank++;
        }
    }

    // Retornar rango calculado
    return rank;
}

/**
 * Función que calcula la transpuesta de una matriz
 * @param m Número de filas de la matriz de entrada
 * @param n Número de columnas de la matriz de entrada
 * @param A Matriz de entrada
 * @param B Matriz de salida que contendrá la transpuesta
 * 
 * La función utiliza paralelización OpenMP con collapse(2) para
 * optimizar el cálculo de la transpuesta.
 */
void transpuesta(int m, int n, double A[m][n], double B[n][m]) {
    #pragma omp parallel for collapse(2)
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            B[j][i] = A[i][j];
        }
    }

    return;
}

/**
 * Función que multiplica dos matrices
 * @param m Número de filas de la primera matriz
 * @param n Número de columnas de la primera matriz y filas de la segunda
 * @param p Número de columnas de la segunda matriz
 * @param A Primera matriz
 * @param B Segunda matriz
 * @param C Matriz de salida que contendrá el resultado
 * 
 * La función utiliza paralelización OpenMP con collapse(2) para
 * optimizar la multiplicación de matrices.
 */
void multiplicar(int m, int n, int p, double A[m][n], double B[n][p], double C[m][p]) {
    #pragma omp parallel for collapse(2)
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < p; j++) {
            C[i][j] = 0;

            for (int k = 0; k < n; k++) {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }

    return;
}

/**
 * Función que calcula la inversa de una matriz usando eliminación gaussiana
 * @param m Dimensión de la matriz cuadrada
 * @param A Matriz de entrada
 * @param inv Matriz donde se almacenará la inversa
 * @return 0 si se calcula la inversa exitosamente, 1 si la matriz no es invertible
 * 
 * La función utiliza paralelización OpenMP para optimizar las operaciones
 * de eliminación gaussiana y cálculo de la inversa.
 */
int inverse(int m, double A[m][m], double inv[m][m]) {
    const double EPS = 1e-16; 
    
    // Crear matriz aumentada [A | I]
    double aug[m][2 * m];

    #pragma omp parallel for
    for (int i = 0; i < m; i++) {
        // Copiar matriz original a la parte izquierda
        for (int j = 0; j < m; j++) {         
            aug[i][j] = A[i][j];
        }

        // Crear matriz identidad en la parte derecha
        for (int j = 0; j < m; j++) {         
            aug[i][j + m] = (i == j) ? 1.0 : 0.0;
        }
    }

    // Para cada columna
    for (int col = 0, row = 0; col < m && row < m; col++, row++) {
        // Buscar el mejor pivote (pivoteo parcial)
        int piv = row;
        for (int k = row + 1; k < m; k++)
            if (fabs(aug[k][col]) > fabs(aug[piv][col]))
                piv = k;

        // Verificar si hay pivote válido
        if (fabs(aug[piv][col]) < EPS) {
            fprintf(stderr, "Error: matriz no invertible (pivote cero)\n");
            return 1;
        }

        // Intercambiar filas si el pivote no está en la posición correcta
        if (piv != row) {
            for (int j = 0; j < 2 * m; j++) {
                double tmp = aug[piv][j];
                aug[piv][j] = aug[row][j];
                aug[row][j] = tmp;
            }
        }

        // Normalizar fila pivote
        double piv_val = aug[row][col];
        for (int j = 0; j < 2 * m; j++)
            aug[row][j] /= piv_val;

        // Reducir las demás filas usando la fila pivote
        #pragma omp parallel for
        for (int i = 0; i < m; i++) {
            if (i == row) continue;
            double factor = aug[i][col];
            for (int j = 0; j < 2 * m; j++)
                aug[i][j] -= factor * aug[row][j];
        }
    }

    for (int i = 0; i < m; i++) {
        for (int j = 0; j < m; j++)
            inv[i][j] = aug[i][j + m];
    }

    return 0;
}

int main(int argc, char* argv[]) {
    // Verificar argumentos de línea de comandos
    if (argc < 2) {
        printf("Uso: %s archivo_entrada [num_threads]\n", argv[0]);
        return 1;
    }

    // Configurar número de hilos (por defecto: núcleos disponibles)
    int num_threads = omp_get_num_procs();

    // Si se especificó número de hilos, usar ese valor
    if (argc > 2) {
        int temp_threads = atoi(argv[2]);
        if (temp_threads > 0) {
            num_threads = temp_threads;
        } else {
            printf("Número de hilos debe ser positivo\n");
            return 1;
        }
    }

    omp_set_num_threads(num_threads);

    // Leer dimensiones y datos de la matriz
    int m, n;
    const char *entrada = argv[1];
    FILE* archivo = fopen(entrada, "r");

    // Verificar apertura exitosa del archivo
    if (!archivo) {
        perror("Error al abrir el archivo");
        return 1;
    }

    // Leer dimensiones de la matriz
    if (fscanf(archivo, "%d %d", &m, &n) != 2) {
        fprintf(stderr, "Error al leer las dimensiones del archivo\n");
        fclose(archivo);
        return 1;
    }

    // Leer elementos de la matriz
    double matriz[m][n];
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            if (fscanf(archivo, "%lf", &matriz[i][j]) != 1) {
                fprintf(stderr, "Error al leer el elemento (%d, %d)\n", i, j);
                fclose(archivo);
                return 1;
            }
        }
    }

    // Cerrar archivo después de leer datos
    fclose(archivo);

    // Calcular rango de la matriz
    int r = rango(m, n, matriz);

    // Verificar si la matriz tiene rango completo
    if (r < (m < n ? m : n)) {
        printf("-1\n");
        return 0;
    }

    // Caso 1: Rango igual al número de filas
    if (r == m) {
        // Calcular transpuesta
        double AT[n][m];
        transpuesta(m, n, matriz, AT);

        // Calcular A^T * A
        double AAT[m][m];
        multiplicar(m, n, m, matriz, AT, AAT);
        
        // Invertir A^T * A
        double invATA[m][m];
        int error = inverse(m, AAT, invATA);

        // Verificar si hubo error en la inversión
        if (error == 1) {
            printf("-1\n");
            return 1;
        }

        // Calcular pseudo-inversa final: (A^T * A)^-1 * A^T
        double P[n][m];
        multiplicar(n, m, m, AT, invATA, P);

        // Imprimir resultado
        printMatrix("R", n, m, P);
    }

    // Caso 2: Rango igual al número de columnas
    if (r == n) {
        // Calcular transpuesta
        double AT[n][m];
        transpuesta(m, n, matriz, AT);

        // Calcular A * A^T
        double ATA[n][n];
        multiplicar(n, m, n, AT, matriz, ATA);

        // Invertir A * A^T
        double invATA[n][n];
        int error = inverse(n, ATA, invATA);

        // Verificar si hubo error en la inversión
        if (error == 1) {
            printf("-1\n");
            return 1;
        }

        double P[n][m];
        multiplicar(n, n, m, invATA, AT, P);

        // Imprimir resultado
        printMatrix("L", n, m, P);
    }

    return 0;
}
