#include <stdio.h>
#include <stdlib.h>
#include <time.h>

int main() {
    // Definir dimensiones de la matriz
    int rows = 800;  // Número de filas
    int cols = 1000;  // Número de columnas
    
    // Abrir archivo para escritura
    FILE *file = fopen("entrada.ent", "w");
    if (file == NULL) {
        printf("Error al abrir el archivo!\n");
        return 1;
    }
    
    // Inicializar generador de números aleatorios
    srand(time(NULL));
    
    // Escribir dimensiones en la primera línea
    fprintf(file, "%d %d\n", rows, cols);
    
    // Generar y escribir la matriz
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            // Generar número double aleatorio entre 0 y 100
            double value = (double)rand() / RAND_MAX * 100.0;
            fprintf(file, "%.6f", value);
            
            // Añadir espacio o salto de línea según corresponda
            if (j < cols - 1) {
                fprintf(file, " ");
            } else {
                fprintf(file, "\n");
            }
        }
    }
    
    // Cerrar el archivo
    fclose(file);
    printf("Archivo generado exitosamente.\n");
    
    return 0;
}