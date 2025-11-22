#include <stdio.h>
#include <stdlib.h>
#include <time.h>

typedef struct {
    int u, v;
} Edge;

void make_tree_random(int **A, int v, Edge *edges, int *edge_count) {
    for (int i = 1; i < v; i++) {
        int parent = rand() % i;   // losowy wierzcholek z [0, i-1]
        A[i][parent] = A[parent][i] = 1;
        edges[(*edge_count)++] = (Edge){parent, i};
    }
}

int main(int argc, char *argv[]) {
    int v = atoi(argv[1]);
    int e = atoi(argv[2]);

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

    // ----- Macierz sąsiedztwa -----
    int **A = malloc(v * sizeof(int *));
    for (int i = 0; i < v; i++) {
        A[i] = calloc(v, sizeof(int));
    }

    Edge *edges = malloc(e * sizeof(Edge));
    int edge_count = 0;

    // ---------------------------------------------------
    // 1) Generowanie losowego drzewa (spojność gwarantowana)
    // ---------------------------------------------------
    make_tree_random(A, v, edges, &edge_count);
    // ---------------------------------------------------
    // 2) Dodawanie losowych brakujących krawędzi aż do e
    // ---------------------------------------------------
    while (edge_count < e) {
        int a = rand() % v;
        int b = rand() % v;

        if (a == b) continue;
        if (A[a][b]) continue;

        A[a][b] = A[b][a] = 1;
        edges[edge_count++] = (Edge){a, b};
    }

    // ----- Wyniki -----
    printf("n=%d\n", v);
    for (int i = 0; i < v; i++) {
        for (int j = 0; j < v; j++) {
            printf("%d ", A[i][j]);
        }
        printf("\n");
    }

    // Czyszczenie
    for (int i = 0; i < v; i++) free(A[i]);
    free(A);
    free(edges);

    return 0;
}
