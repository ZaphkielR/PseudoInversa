//Alumnos: Christian Delgado, Rafael Morales
//Fecha: 10-06-2025
//Profesor: Sergio Antonio Baltierra Valenzuela
//Carrera: Ingenieria Civil Informatica
//Ramo: Computación de Alto Rendimiento
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


void printMatrix(char* label, int m, int n, double matriz[m][n]) {    
    FILE* archivo = fopen("salida_sec.sal", "w");

    if (!archivo) {
        perror("Error al abrir el archivo");
        return;
    }

    fprintf(archivo, "%s\n", label);
    
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            fprintf(archivo, "%lf ", matriz[i][j]);
        }
        fprintf(archivo, "\n");
    } 

    fclose(archivo);
}

int rango(int m, int n, double A[m][n]) {
    double temp[m][n];
    int rank = 0;

    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            temp[i][j] = A[i][j];
        }
    }

    for (int i = 0; i < n; i++) {
        int pivotRow = -1;

        for (int j = rank; j < m; j++) {
            if (fabs(temp[j][i]) > 1e-10) {
                pivotRow = j;
                break;
            }
        }

        if (pivotRow != -1) {
            for (int k = 0; k < n; k++) {
                int tmp = temp[rank][k];
                temp[rank][k] = temp[pivotRow][k];
                temp[pivotRow][k] = tmp;
            }

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

void transpuesta(int m, int n, double A[m][n], double B[n][m]) {
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            B[j][i] = A[i][j];
        }
    }
    return;
}

void multiplicar(int m, int n, int p, double A[m][n], double B[n][p], double C[m][p]) {
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

int inverse(int m, double A[m][m], double inv[m][m]) {
    const double EPS = 1e-12;  

    double aug[m][2 * m];

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

        for (int k = row + 1; k < m; k++)
            if (fabs(aug[k][col]) > fabs(aug[piv][col]))
                piv = k;

        if (fabs(aug[piv][col]) < EPS) {
            fprintf(stderr, "Error: matriz no invertible (pivote cero)\n");
            return 1;
        }

        if (piv != row) {
            for (int j = 0; j < 2 * m; j++) {
                double tmp = aug[piv][j];
                aug[piv][j] = aug[row][j];
                aug[row][j] = tmp;
            }
        }   

        double piv_val = aug[row][col];
        for (int j = 0; j < 2 * m; j++)
            aug[row][j] /= piv_val;

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

    int m, n;

    const char *entrada = argv[1]; 

    FILE* archivo = fopen(entrada, "r");

    if (!archivo) {
        perror("Error al abrir el archivo");
        return 1;
    }

    if (fscanf(archivo, "%d %d", &m, &n) != 2) {
        fprintf(stderr, "Error al leer las dimensiones del archivo\n");
        fclose(archivo);
        return 1;
    }

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

    fclose(archivo);


    int r = rango(m, n, matriz);

    if (r < (m < n ? m : n)) {
        printf("-1\n");
        return 0;
    }

    if (r == m) {
        double AT[n][m];
        transpuesta(m, n, matriz, AT);

        double ATA[m][m];
        multiplicar(m, n, m, matriz, AT, ATA);
        
        double invATA[m][m];
        int error = inverse(m, ATA, invATA);

        if (error == 1) {
            printf("-1\n");
            return 1;
        }

        double P[n][m];
        multiplicar(n, m, m, AT, invATA, P);

        printMatrix("R", n, m, P);
    }

    if (r == n) {
        double AT[n][m];
        transpuesta(m, n, matriz, AT);

        double ATA[n][n];
        multiplicar(n, m, n, AT, matriz, ATA);
        
        double invATA[n][n];
        inverse(n, ATA, invATA);

        double P[n][m];
        multiplicar(n, n, m, invATA, AT, P);

        printMatrix("L", n, m, P);
    }


    return 0;
}