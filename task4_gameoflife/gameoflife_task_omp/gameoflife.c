#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <omp.h>

#define for_x for (int x = 0; x < w; x++)
#define for_y for (int y = 0; y < h; y++)
#define for_xy for_x for_y

void init(void *u, int w, int h) {

    int(*univ)[w] = u;
    
    for_xy univ[y][x] = rand() < RAND_MAX / 10 ? 1 : 0;

}

void show(void *u, int w, int h) {
    int(*univ)[w] = u;
    printf("\033[H");
    for_y {
        for_x printf(univ[y][x] ? "\033[07m  \033[m" : "  ");
        printf("\033[E");
    }
    fflush(stdout);
}

void evolve(void *u, int w, int h) {

    unsigned(*univ)[w] = u;
    unsigned new[h][w];

#pragma omp parallel 
{ 
#pragma omp single
{
    for_y 
#pragma omp task
{
        for_x 
        {
            int n = 0;
            for (int y1 = y - 1; y1 <= y + 1; y1++)
                for (int x1 = x - 1; x1 <= x + 1; x1++)
                    if (univ[(y1 + h) % h][(x1 + w) % w]) n++;

            if (univ[y][x]) n--;
            new[y][x] = (n == 3 || (n == 2 && univ[y][x]));
        }
} // task

} // single
#pragma omp taskwait

#pragma omp single
{
    for_y 
#pragma omp task
{    
        for_x 
            univ[y][x] = new[y][x];
} // task

} // single

} // region



}

void game(unsigned *u, int w, int h, int iter) {
    double time1_local, time2_local, elapsed_local;

    for (int i = 0; i < iter; i++) {
#ifdef LIFE_VISUAL
        show(u, w, h);
#endif

    if(i == iter / 2) {
        time1_local = omp_get_wtime();
    }
        evolve(u, w, h);

    if(i == iter / 2) {

        time2_local = omp_get_wtime();
        elapsed_local = time2_local - time1_local;
        printf("Evolve function time > %f \n", elapsed_local);
    }

#ifdef LIFE_VISUAL
        usleep(200000);
#endif
    }
}

int main(int c, char *v[]) {

    double time1_global, time2_global, elapsed_global;

    time1_global = omp_get_wtime();

    int w = 0, h = 0, iter = 0;
    unsigned *u;

    if (c > 1) w = atoi(v[1]);
    if (c > 2) h = atoi(v[2]);
    if (c > 3) iter = atoi(v[3]);
    if (w <= 0) w = 30;
    if (h <= 0) h = 30;
    if (iter <= 0) iter = 1000;
    
    printf("Number of iteration %d \n", iter);

    u = (unsigned *)malloc(w * h * sizeof(unsigned));
    if (!u) {
        printf("No memory!\n");
        exit(1);
    }

    init(u, w, h);
    game(u, w, h, iter);

    free(u);

    time2_global = omp_get_wtime();
    elapsed_global = time2_global - time1_global;
    printf("Total simulation time > %f \n\n", elapsed_global);


}
