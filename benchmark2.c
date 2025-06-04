#include <stdio.h>
#include <stdlib.h>
#include <sys/wait.h>
#include <unistd.h>
#include <time.h>

#define MAX_THREADS 12

int main() {
    FILE *output = fopen("tiempos.txt", "w");
    if (!output) {
        perror("Error al abrir archivo de salida");
        return 1;
    }

    // Compilar versión mono-hilo
    printf("Compilando versión mono-hilo...\n");
    system("gcc mono-hilo.c -o mono-hilo -lm");


    // Ejecutar versión mono-hilo
    printf("\nEjecutando versión mono-hilo:\n");
    fprintf(output, "Version mono-hilo:\n");
    
    clock_t start = clock();
    system("./mono-hilo entrada.ent");
    clock_t end = clock();
    
    double mono_time = (double)(end - start) / CLOCKS_PER_SEC;
    printf("Tiempo de ejecución: %.6f segundos\n", mono_time);
    fprintf(output, "Tiempo de ejecución: %.6f segundos\n", mono_time);


    fclose(output);
    printf("\nResultados guardados en tiempos.txt\n");
    return 0;
}
