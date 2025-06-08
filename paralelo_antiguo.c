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
    double temp[m][n];
    int rank = 0;

    // Copiar matriz (puede paralelizarse si m y n son grandes)
    #pragma omp parallel for if(m * n > 10000)
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            temp[i][j] = A[i][j];
        }
    }

    for (int i = 0; i < n; i++) {
        int pivotRow = -1;
        double maxVal = 0.0;

        // Paralelizar la búsqueda del pivote
        #pragma omp parallel
        {
            int local_pivot = -1;
            double local_max = 0.0;

            #pragma omp for
            for (int j = rank; j < m; j++) {
                if (fabs(temp[j][i]) > local_max) {
                    local_max = fabs(temp[j][i]);
                    local_pivot = j;
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

        if (pivotRow != -1) {
            // Intercambio de filas (secuencial, pero rápido)
            for (int k = 0; k < n; k++) {
                double tmp = temp[rank][k];
                temp[rank][k] = temp[pivotRow][k];
                temp[pivotRow][k] = tmp;
            }

            // Eliminación gaussiana paralela
            #pragma omp parallel for
            for (int j = 0; j < m; j++) {
                if (j != rank) {
                    double factor = temp[j][i] / temp[rank][i];
                    for (int k = 0; k < n; k++) {
                        temp[j][k] -= factor * temp[rank][k];
                    }
                }
            }
            rank++;
        }
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
    // Asignación de memoria para la matriz resultante C
    double** C = (double**)malloc(m * sizeof(double*));
    for (int i = 0; i < m; i++) {
        C[i] = (double*)malloc(p * sizeof(double));
        for (int j = 0; j < p; j++) {
            C[i][j] = 0; // Inicializar C para evitar basura
        }
    }

    // Región paralela con variables compartidas
    #pragma omp parallel shared(A, B, C, m, n, p)
    {
        // Obtener el número total de hilos y el ID del hilo actual
        int num_threads = omp_get_num_threads();
        int thread_id = omp_get_thread_num();

        // Distribución intercalada: cada hilo procesa filas i donde i = thread_id + k*num_threads
        for (int i = thread_id; i < m; i += num_threads) {
            for (int j = 0; j < p; j++) {
                for (int k = 0; k < n; k++) {
                    C[i][j] += A[i][k] * B[k][j];
                }
            }
        }
    }

    return C;
}

// Inversa Gauss-Jordan paralela
double** inverse(int m, double** A) {
    const double EPS = 1e-12;
    double** aug = (double**)malloc(m * sizeof(double*));
    for (int i = 0; i < m; i++) {
        aug[i] = (double*)malloc(2 * m * sizeof(double));
    }

    #pragma omp parallel for
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < m; j++) {
            aug[i][j] = A[i][j];
        }
        for (int j = 0; j < m; j++) {
            aug[i][j + m] = (i == j) ? 1.0 : 0.0;
        }
    }

    for (int col = 0, row = 0; col < m && row < m; col++, row++) {
        int piv = row;
        double maxVal = fabs(aug[row][col]);

        // Paralelizar la búsqueda del pivote
        #pragma omp parallel shared(aug, piv, maxVal, row, col, m)
        {
            int local_piv = row;
            double local_max = maxVal;

            // Distribución intercalada para buscar el pivote
            int num_threads = omp_get_num_threads();
            int thread_id = omp_get_thread_num();
            for (int k = row + thread_id; k < m; k += num_threads) {
                if (fabs(aug[k][col]) > local_max) {
                    local_max = fabs(aug[k][col]);
                    local_piv = k;
                }
            }

            // Reducir el máximo local al global
            #pragma omp critical
            {
                if (local_max > maxVal) {
                    maxVal = local_max;
                    piv = local_piv;
                }
            }
        }

        if (fabs(aug[piv][col]) < EPS) {
            fprintf(stderr, "Error: matriz no invertible\n");
            return NULL;
        }

        if (piv != row) {
            double* tmp = aug[piv];
            aug[piv] = aug[row];
            aug[row] = tmp;
        }

        double piv_val = aug[row][col];
        #pragma omp parallel for
        for (int j = 0; j < 2 * m; j++) {
            aug[row][j] /= piv_val;
        }

        // Eliminación de filas en paralelo con distribución intercalada
        #pragma omp parallel shared(aug, row, col, m)
        {
            int num_threads = omp_get_num_threads();
            int thread_id = omp_get_thread_num();

            // Cada hilo procesa filas de forma intercalada
            for (int i = thread_id; i < m; i += num_threads) {
                if (i != row) {
                    double factor = aug[i][col];
                    for (int j = 0; j < 2 * m; j++) {
                        aug[i][j] -= factor * aug[row][j];
                    }
                }
            }
        }
    }

    double** inv = (double**)malloc(m * sizeof(double*));
    #pragma omp parallel for
    for (int i = 0; i < m; i++) {
        inv[i] = (double*)malloc(m * sizeof(double));
        for (int j = 0; j < m; j++) {
            inv[i][j] = aug[i][j + m];
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
