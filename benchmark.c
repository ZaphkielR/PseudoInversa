#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <unistd.h>
#include <math.h>

double get_time() {
    struct timeval t;
    gettimeofday(&t, NULL);
    return t.tv_sec + t.tv_usec / 1e6;
}

int main() {
    
    char ENTRADA[] = "entrada_1.ent";

    // Compilar versión secuencial
    printf("Compilando version secuencial...\n");
    system("gcc secuencial.c -o secuencial -lm");

    // =============================
    // MÉTRICAS: SECUENCIAL
    // =============================
    char comando_sec[128];
    sprintf(comando_sec, "secuencial %s", ENTRADA);

    printf("\nEjecutando secuencial\n");
    double start = get_time();
    system(comando_sec);
    double end = get_time();
    
    double sec_time = end - start;
    printf("Tiempo de ejecucion secuencial: %.6f segundos\n", sec_time);

    // =============================
    // MÉTRICAS: PARALELO
    // =============================
    
    double* par_time;

    int N = 10;
    par_time = (double*)malloc(N * sizeof(double));

    system("gcc -fopenmp paralelo.c -o paralelo -lm");
    printf("\nCompilando version paralela...\n");


    for (int i = 1; i <= N; i++) {
        int num_threads = (int)pow(2, i);

        char comando[128];

        sprintf(comando, "paralelo %s %d", ENTRADA, num_threads);
        printf("\nEjecutando con %d hilos...\n", num_threads);

        start = get_time();
        system(comando);
        end = get_time();

        par_time[i - 1] = end - start;
        printf("Tiempo de ejecucion paralela con %d hilos: %.6f segundos\n", num_threads, par_time[i - 1]);
    }


    FILE *output = fopen("metricas.met", "w");

    if (!output) {
        perror("Error al abrir archivo de salida");
        return 1;
    }

    fprintf(output, "Ensayo\t Hilos\t\t Speedup\t\t\t Eficiencia\n");
    for (int i = 0; i < N; i++) {
        fprintf(output, "%d\t\t 2 ^ (%d)\t %.40f\t %.40f\n", i+1, i+1, sec_time / par_time[i], (sec_time / par_time[i]) / pow(2, i+1));
    }    
    fclose(output);
    printf("\n| Resultados guardados en metricas.met |\n");
    return 0;
}
