#include <stdio.h>
#include <stdlib.h>
#include <sys/wait.h>
#include <unistd.h>
#include <time.h>
#include <math.h>

int main() {
    
    char ENTRADA[] = "entrada.ent";

    // Compilar versión secuencial
    printf("Compilando versión secuencial...\n");
    system("gcc secuencial.c -o secuencial -lm");

    // =============================
    // MÉTRICAS: SECUENCIAL
    // =============================
    char comando_sec[128];
    sprintf(comando_sec, "./secuencial %s", ENTRADA);

    printf("\nEjecutando secuencial\n");
    clock_t start = clock();
    system(comando_sec);
    clock_t end = clock();
    
    double sec_time = (double)(end - start) / CLOCKS_PER_SEC;
    printf("Tiempo de ejecución secuencial: %.6f segundos\n", sec_time);

    // =============================
    // MÉTRICAS: PARALELO
    // =============================
    
    double* par_time;

    int N = 12;
    par_time = (double*)malloc(N * sizeof(double));

    for (int i = 1; i <= N; i++) {
        printf("\nCompilando versión paralela con %d hilos...\n", (int)pow(2, i));
        system("gcc -fopenmp paralelo.c -o paralelo -lm");

        char comando[128];
        sprintf(comando, "./paralelo %s %d", ENTRADA, (int)pow(2, i));
        printf("Ejecutando con %d hilos...\n", (int)pow(2, i));

        start = clock();
        system(comando);
        end = clock();

        par_time[i - 1] = (double)(end - start) / CLOCKS_PER_SEC;
        printf("Tiempo de ejecución paralela con %d hilos: %.6f segundos\n", (int)pow(2, i), par_time[i - 1]);
    }


    FILE *output = fopen("metricas.met", "w");

    if (!output) {
        perror("Error al abrir archivo de salida");
        return 1;
    }

    fprintf(output, "Ensayo\t Hilos\t\t Speedup\t\t\t Eficiencia\n");
    for (int i = 0; i < N; i++) {
        fprintf(output, "%d\t\t 2 ^ (%d)\t %.15f\t %.15f\n", i+1, i+1, sec_time / par_time[i], (sec_time / par_time[i]) / pow(2, i+1));
    }    
    fclose(output);
    printf("\n| Resultados guardados en metricas.met |\n");
    return 0;
}
