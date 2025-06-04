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


    // =============================
    // MÉTRICAS: MONO-HILO
    // =============================
    printf("\nEjecutando versión mono-hilo:\n");
    fprintf(output, "Ensayo Hilos Speedup Eficiencia\n");
    
    clock_t start = clock();
    system("./mono-hilo entrada.ent");
    clock_t end = clock();
    
    double mono_time = (double)(end - start) / CLOCKS_PER_SEC;
    printf("Tiempo de ejecución mono-hilo: %.6f segundos\n", mono_time);
    // Primer ensayo: mono-hilo, speedup y eficiencia son 1
    fprintf(output, "%d %d %.15f %.15f\n", 1, 1, 1.0, 1.0);

    // =============================
    // MÉTRICAS: PARALELO
    // =============================
    // Compilar versión paralela
    printf("\nCompilando versión paralela...\n");
    system("gcc paralelo.c -o paralelo -fopenmp -lm");

    // Ejecutar versión paralela con diferentes hilos (potencias de 2)
    for (int ensayo = 1; ensayo < 10; ensayo++) {
        int hilos = 1 << ensayo; // 2, 4, 8, ... 512
        char cmd[128];
        snprintf(cmd, sizeof(cmd), "OMP_NUM_THREADS=%d ./paralelo entrada.ent", hilos);

        clock_t p_start = clock();
        system(cmd);
        clock_t p_end = clock();

        double par_time = (double)(p_end - p_start) / CLOCKS_PER_SEC;
        double speedup = mono_time / par_time;
        double eficiencia = speedup / hilos;

        fprintf(output, "%d %d %.15f %.15f\n", ensayo + 1, hilos, speedup, eficiencia);
        printf("Ensayo %d, Hilos: %d, Speedup: %.6f, Eficiencia: %.6f\n", ensayo + 1, hilos, speedup, eficiencia);
    }

    fclose(output);
    printf("\nResultados guardados en tiempos.txt\n");
    return 0;
}
