#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "gsl/gsl_eigen.h"
#include "gsl/gsl_matrix.h"
#include "gsl/gsl_vector.h"
#include <math.h>

const double precision = 1e-6;

typedef struct{
    int u;
    int o;
    double score;
}Edge_score;

void make_tree_random(int **A, int v, int *edge_count) {
    for (int i = 1; i < v; i++) {
        int parent = rand() % i;   // losowy wierzcholek z [0, i-1]
        A[i][parent] = A[parent][i] = 1;
        edge_count++;
    }
}

void get_eigenvalues(int **A, int v, double *eigenvalues) {
    gsl_matrix *gsl_A = gsl_matrix_alloc(v, v);
    for (int i = 0; i < v; i++) {
        for (int j = 0; j < v; j++) {
            gsl_matrix_set(gsl_A, i, j, A[i][j]);
        }
    }

    gsl_vector *eval = gsl_vector_alloc(v);
    gsl_matrix *evec = gsl_matrix_alloc(v, v);
    gsl_eigen_symmv_workspace *w = gsl_eigen_symmv_alloc(v);

    gsl_eigen_symmv(gsl_A, eval, evec, w);
    gsl_eigen_symmv_free(w);

    for (int i = 0; i < v; i++) {
        eigenvalues[i] = gsl_vector_get(eval, i);
    }

    gsl_vector_free(eval);
    gsl_matrix_free(evec);
    gsl_matrix_free(gsl_A);
}

void get_score(int **A, int v, double *score) { // Mniejsza wartosc lepsza
    double *eigenvalues = malloc(v * sizeof(double));
    get_eigenvalues(A, v, eigenvalues);

    for (int i = 0; i < v; i++) {
        *score += fabs(eigenvalues[i] - round(eigenvalues[i]));
    }

    free(eigenvalues);

}

int is_integral_graph(int **A, int v) {
    double *eigenvalues = malloc(v * sizeof(double));
    get_eigenvalues(A, v, eigenvalues);

    for (int i = 0; i < v; i++) {
        if (fabs(eigenvalues[i] - round(eigenvalues[i])) > precision) {
            free(eigenvalues);
            return 0; 
        }
    }

    free(eigenvalues);
    return 1; 
}



int main(int argc, char *argv[]) {
    int v = atoi(argv[1]);
    int e = atoi(argv[2]);
    int greed_number = atoi(argv[3]);
    int number_of_trials = atoi(argv[4]);

    if (v < 1) {
        printf("Blad: v musi byc >= 1\n");
        return 1;
    }

    int e_min = v - 1;
    int e_max = v * (v - 1) / 2;

    if (e < e_min || e > e_max) {
        printf("Nie da sie stworzyc spojnego grafu o %d wierzcholkach i %d krawedziach.\n", v, e);
        return 1;
    }

    srand(time(NULL));
    int licznik_grafow = 0;
    int ***Graphs = malloc(number_of_trials * sizeof(int **)); // Tablica na grafy

    // Macierz sąsiedztwa 
    for (int g = 0; g < number_of_trials; g++) {
        int **A = malloc(v * sizeof(int *));
        for (int i = 0; i < v; i++) {
            A[i] = calloc(v, sizeof(int));
        }

        int edge_count = 0;

        make_tree_random(A, v, &edge_count);

        while (edge_count < e) {
            Edge_score *Score = malloc(greed_number * sizeof(Edge_score));
            for (int i = 0; i < greed_number; i++) {
                int potential_u = rand() % v;
                int potential_o = rand() % v;
                while (potential_u == potential_o || A[potential_u][potential_o] == 1) {
                    potential_u = rand() % v;
                    potential_o = rand() % v;
                }
                A[potential_u][potential_o] = A[potential_o][potential_u] = 1;
                double score = 0.0;
                get_score(A, v, &score);
                Score[i].u = potential_u;
                Score[i].o = potential_o;
                Score[i].score = score;
                A[potential_u][potential_o] = A[potential_o][potential_u] = 0;
            }

            // Wybierz najlepsza krawedz
            int best_index = 0;
            for (int i = 1; i < greed_number; i++) {
                if (Score[i].score < Score[best_index].score) {
                    best_index = i;
                }
            }
            A[Score[best_index].u][Score[best_index].o] = A[Score[best_index].o][Score[best_index].u] = 1;
            edge_count++;

        }
        if(is_integral_graph(A, v)) {
            licznik_grafow++;
            printf("Wykopałeś coś w: %d\n", g);
        }
        Graphs[g] = A;
    }
    // Wyniki 
    printf("n=%d\n", v);
    for (int k = 0; k < number_of_trials; k++) {
        printf("Graf %d:\n", k + 1);
        for (int i = 0; i < v; i++) {
            for (int j = 0; j < v; j++) {
                printf("%d ", Graphs[k][i][j]);
            }
            printf("\n");
        }
    }
    printf("Liczba wykopanych grafow calkowitoliczbowych: %d\n", licznik_grafow);

    // Czyszczenie
    for (int i = 0; i < number_of_trials; i++) free(Graphs[i]);
    free(Graphs);

    return 0;
}
