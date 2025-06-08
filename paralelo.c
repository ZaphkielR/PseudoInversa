#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

// Imprimir matriz a archivo
void printMatrix(double** matriz, char* label, int m, int n) {    
    FILE* archivo = fopen("salida.sal", "w");

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

double rango(int m, int n, double** A) {
    // Matriz de trabajo en pila (si m·n no es demasiado grande)
    double temp[m][n];
    int rank = 0;

    // 1) Copiar A → temp en paralelo
    #pragma omp parallel shared(A, temp, m, n)
    {
        int num_threads = omp_get_num_threads();
        int tid = omp_get_thread_num();
        for (int i = tid; i < m; i += num_threads) {
            for (int j = 0; j < n; j++) {
                temp[i][j] = A[i][j];
            }
        }
    }

    // 2) Eliminación gaussiana con pivoteo
    for (int col = 0; col < n; col++) {
        int pivotRow = -1;
        double maxVal = 0.0;

        // Búsqueda del mayor valor absoluto en la columna 'col'
        #pragma omp parallel shared(temp, pivotRow, maxVal, m, col, rank)
        {
            int num_threads = omp_get_num_threads();
            int tid = omp_get_thread_num();
            double local_max = 0.0;
            int local_pivot = -1;

            for (int i = rank + tid; i < m; i += num_threads) {
                double v = fabs(temp[i][col]);
                if (v > local_max) {
                    local_max = v;
                    local_pivot = i;
                }
            }

            #pragma omp critical
            {
                if (local_max > maxVal) {
                    maxVal = local_max;
                    pivotRow = local_pivot;
                }
            }
        }

        if (pivotRow == -1 || maxVal < 1e-12) {
            // Todos los elementos son (casi) cero → no hay más pivotes
            break;
        }

        // 3) Intercambio de filas (rank) ↔ (pivotRow)
        if (pivotRow != rank) {
            for (int j = 0; j < n; j++) {
                double tmp = temp[rank][j];
                temp[rank][j] = temp[pivotRow][j];
                temp[pivotRow][j] = tmp;
            }
        }

        // 4) Eliminación de todas las otras filas en paralelo
        #pragma omp parallel shared(temp, m, n, rank, col)
        {
            int num_threads = omp_get_num_threads();
            int tid = omp_get_thread_num();

            for (int i = tid; i < m; i += num_threads) {
                if (i == rank) continue;
                double factor = temp[i][col] / temp[rank][col];
                for (int j = col; j < n; j++) {
                    temp[i][j] -= factor * temp[rank][j];
                }
            }
        }

        rank++;
    }

    return rank;
}

// Transpuesta paralela
double** transpuesta(int m, int n, double** A) {
    double** B = (double**)malloc(n * sizeof(double*));
    for (int i = 0; i < n; i++) {
        B[i] = (double*)malloc(m * sizeof(double));
    }

    #pragma omp parallel for collapse(2)
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            B[j][i] = A[i][j];
        }
    }

    return B;
}



// Multiplicación paralela
double** multiplicar(int m, int n, int p, double** A, double** B) {
    // Reserva de C y puesta a cero
    double** C = malloc(m * sizeof(double*));
    for (int i = 0; i < m; i++) {
        C[i] = calloc(p, sizeof(double));
    }

    // Zona paralela: reparto intercalado de filas entre hilos
    #pragma omp parallel shared(A, B, C, m, n, p)
    {
        int num_threads = omp_get_num_threads();
        int tid = omp_get_thread_num();

        for (int i = tid; i < m; i += num_threads) {
            for (int k = 0; k < n; k++) {
                double a_ik = A[i][k];
                for (int j = 0; j < p; j++) {
                    C[i][j] += a_ik * B[k][j];
                }
            }
        }
    }

    return C;
}

// Inversa Gauss-Jordan paralela
double** inverse(int m, double** A) {
    const double EPS = 1e-12;
    // Matriz aumentada [A | I]
    double** aug = malloc(m * sizeof(double*));
    for (int i = 0; i < m; i++) {
        aug[i] = malloc(2 * m * sizeof(double));
    }

    // Inicializar en paralelo
    #pragma omp parallel shared(aug, A, m)
    {
        int num_threads = omp_get_num_threads();
        int tid = omp_get_thread_num();
        for (int i = tid; i < m; i += num_threads) {
            for (int j = 0; j < m; j++) {
                aug[i][j] = A[i][j];
                aug[i][j + m] = (i == j) ? 1.0 : 0.0;
            }
        }
    }

    // Proceso de eliminación Gauss–Jordan
    for (int col = 0, row = 0; col < m && row < m; col++, row++) {
        int piv = row;
        double maxVal = fabs(aug[row][col]);

        // Búsqueda de pivote en columna 'col'
        #pragma omp parallel shared(aug, piv, maxVal, m, row, col)
        {
            int num_threads = omp_get_num_threads();
            int tid = omp_get_thread_num();
            double local_max = maxVal;
            int local_piv = piv;

            for (int i = row + tid; i < m; i += num_threads) {
                double v = fabs(aug[i][col]);
                if (v > local_max) {
                    local_max = v;
                    local_piv = i;
                }
            }

            #pragma omp critical
            {
                if (local_max > maxVal) {
                    maxVal = local_max;
                    piv = local_piv;
                }
            }
        }

        if (fabs(aug[piv][col]) < EPS) {
            fprintf(stderr, "Error: matriz no invertible en columna %d\n", col);
            return NULL;
        }

        // Swap de filas piv ↔ row
        if (piv != row) {
            double* tmp = aug[piv];
            aug[piv] = aug[row];
            aug[row] = tmp;
        }

        // Normalizar fila 'row'
        double diag = aug[row][col];
        #pragma omp parallel for shared(aug, m, row, col)
        for (int j = 0; j < 2*m; j++) {
            aug[row][j] /= diag;
        }

        // Eliminar resto de filas
        #pragma omp parallel shared(aug, m, row, col)
        {
            int num_threads = omp_get_num_threads();
            int tid = omp_get_thread_num();

            for (int i = tid; i < m; i += num_threads) {
                if (i == row) continue;
                double factor = aug[i][col];
                for (int j = 0; j < 2*m; j++) {
                    aug[i][j] -= factor * aug[row][j];
                }
            }
        }
    }

    // Extraer la inversa de la parte derecha
    double** inv = malloc(m * sizeof(double*));
    #pragma omp parallel shared(inv, aug, m)
    {
        int num_threads = omp_get_num_threads();
        int tid = omp_get_thread_num();
        for (int i = tid; i < m; i += num_threads) {
            inv[i] = malloc(m * sizeof(double));
            for (int j = 0; j < m; j++) {
                inv[i][j] = aug[i][j + m];
            }
        }
    }

    return inv;
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        printf("Uso: %s archivo_entrada [num_threads]\n", argv[0]);
        return 1;
    }

    int num_threads = omp_get_num_procs();

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

    int m, n;

    const char *entrada = argv[1];
    double** matriz = leer_matriz(entrada, &m, &n);

    if (!matriz) {
        printf("No se pudo leer la matriz.\n");
        return 1;
    }

    int r = rango(m, n, matriz);
    if (r < (m < n ? m : n)) {
        printf("-1\n");
        return 0;
    }

    if (r == m) {
        double** AT = transpuesta(m, n, matriz);
        double** ATA = multiplicar(m, n, m, matriz, AT);
        double** invATA = inverse(m, ATA);
        double** P = multiplicar(n, m, m, AT, invATA);

        printMatrix(P, "R", n, m);
    }

    if (r == n) {
        double** AT = transpuesta(m, n, matriz);
        double** ATA = multiplicar(n, m, n, AT, matriz);
        double** invATA = inverse(n, ATA);
        double** P = multiplicar(n, n, m, invATA, AT);

        printMatrix(P, "L", n, m);
    }

    return 0;
}
