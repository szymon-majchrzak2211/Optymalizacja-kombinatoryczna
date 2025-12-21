/**
 * @file greedy_random.c
 * @brief Algorytm zachłanny ulosowiony generowania grafów całkowitych (integral graphs)
 *
 * Program znajduje spójny graf całkowity o zadanej liczbie
 * wierzchołków i krawędzi, wykorzystując algorytm zachłanny ulosowiony.
 *
 * @par Argumenty linii poleceń:
 * - <b>v</b> – liczba wierzchołków grafu
 * - <b>e</b> – docelowa liczba krawędzi
 * - <b>greed_number</b> – liczba potencjalnych krawędzi rozważanych na każdym kroku
 * - <b>number_of_trials</b> – liczba wygenerowanych grafów
 * - <b>outfile</b> – plik wyjściowy do zapisu wygenerowanych grafów w formacie graph6
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "gsl/gsl_eigen.h"
#include "gsl/gsl_matrix.h"
#include "gsl/gsl_vector.h"

const double EPS = 1e-6;

/**
 * @struct Edge_score
 * @brief Struktura przechowująca ocenę potencjalnej krawędzi
 */
typedef struct {
    int u;          /**< pierwszy wierzchołek krawędzi */
    int o;          /**< drugi wierzchołek krawędzi */
    double score;   /**< wartość funkcji celu po dodaniu krawędzi */
} Edge_score;

/**
 * @brief Generuje losowe drzewo rozpinające
 *
 * Zapewnia spójność grafu poprzez utworzenie losowego drzewa.
 *
 * @param A macierz sąsiedztwa
 * @param v liczba wierzchołków
 * @param edge_count wskaźnik na licznik krawędzi
 */
void make_tree_random(int **A, int v, int *edge_count) {
    for (int i = 1; i < v; i++) {
        int parent = rand() % i;
        A[i][parent] = A[parent][i] = 1;
        (*edge_count)++;
    }
}

/**
 * @brief Oblicza wartości własne macierzy sąsiedztwa
 *
 * @param A macierz sąsiedztwa
 * @param v liczba wierzchołków
 * @param eigenvalues tablica wynikowa wartości własnych
 */
void get_eigenvalues(int **A, int v, double *eigenvalues) {
    gsl_matrix *gsl_A = gsl_matrix_alloc(v, v);

    for (int i = 0; i < v; i++)
        for (int j = 0; j < v; j++)
            gsl_matrix_set(gsl_A, i, j, A[i][j]);

    gsl_vector *eval = gsl_vector_alloc(v);
    gsl_matrix *evec = gsl_matrix_alloc(v, v);
    gsl_eigen_symmv_workspace *w = gsl_eigen_symmv_alloc(v);

    gsl_eigen_symmv(gsl_A, eval, evec, w);
    gsl_eigen_symmv_free(w);

    for (int i = 0; i < v; i++)
        eigenvalues[i] = gsl_vector_get(eval, i);

    gsl_vector_free(eval);
    gsl_matrix_free(evec);
    gsl_matrix_free(gsl_A);
}

/**
 * @brief Oblicza funkcję celu grafu
 *
 * Funkcja celu to suma odchyleń wartości własnych od liczb całkowitych.
 * Mniejsza wartość oznacza lepszy graf.
 *
 * @param A macierz sąsiedztwa
 * @param v liczba wierzchołków
 * @param score wskaźnik na wartość funkcji celu
 */
void get_score(int **A, int v, double *score) {
    double *eigenvalues = malloc(v * sizeof(double));
    get_eigenvalues(A, v, eigenvalues);

    for (int i = 0; i < v; i++)
        *score += fabs(eigenvalues[i] - round(eigenvalues[i]));

    free(eigenvalues);
}

/**
 * @brief Sprawdza czy graf jest całkowitoliczbowy (integral graph)
 *
 * @param A macierz sąsiedztwa
 * @param v liczba wierzchołków
 * @return 1 jeśli graf jest całkowitoliczbowy, 0 w przeciwnym razie
 */
int is_integral_graph(int **A, int v) {
    double *eigenvalues = malloc(v * sizeof(double));
    get_eigenvalues(A, v, eigenvalues);

    for (int i = 0; i < v; i++) {
        if (fabs(eigenvalues[i] - round(eigenvalues[i])) > EPS) {
            free(eigenvalues);
            return 0;
        }
    }

    free(eigenvalues);
    return 1;
}

/**
 * @brief Konwertuje macierz sąsiedztwa do formatu graph6
 *
 * @param A macierz sąsiedztwa
 * @param n liczba wierzchołków
 * @param output bufor wyjściowy (string)
 */
void adjacency_to_graph6(int **A, int n, char *output) {
    int bit_count = 0;
    int value = 0;
    int out_index = 0;

    output[out_index++] = (char)(n + 63);

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < i; j++) {
            value = (value << 1) | (A[i][j] ? 1 : 0);
            bit_count++;

            if (bit_count == 6) {
                output[out_index++] = (char)(value + 63);
                bit_count = 0;
                value = 0;
            }
        }
    }

    if (bit_count > 0) {
        value <<= (6 - bit_count);
        output[out_index++] = (char)(value + 63);
    }

    output[out_index] = '\0';
}

/**
 * @brief Funkcja główna programu
 *
 * @param argc liczba argumentów
 * @param argv tablica argumentów
 * @return kod zakończenia programu
 */
int main(int argc, char *argv[]) {
    int v = atoi(argv[1]);
    int e = atoi(argv[2]);
    int greed_number = atoi(argv[3]);
    int number_of_trials = atoi(argv[4]);
    char *outfile = argv[5];

    char graph6_output[1000];

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

    int ***Graphs = malloc(number_of_trials * sizeof(int **));

    for (int g = 0; g < number_of_trials; g++) {
        int **A = malloc(v * sizeof(int *));
        for (int i = 0; i < v; i++)
            A[i] = calloc(v, sizeof(int));

        int edge_count = 0;
        make_tree_random(A, v, &edge_count);

        while (edge_count < e) {
            Edge_score *Score = malloc(greed_number * sizeof(Edge_score));

            for (int i = 0; i < greed_number; i++) {
                int potential_u, potential_o;
                do {
                    potential_u = rand() % v;
                    potential_o = rand() % v;
                } while (potential_u == potential_o || A[potential_u][potential_o]);

                A[potential_u][potential_o] = A[potential_o][potential_u] = 1;

                double score = 0.0;
                get_score(A, v, &score);

                Score[i] = (Edge_score){potential_u, potential_o, score};

                A[potential_u][potential_o] = A[potential_o][potential_u] = 0;
            }

            int best = 0;
            for (int i = 1; i < greed_number; i++)
                if (Score[i].score < Score[best].score)
                    best = i;

            A[Score[best].u][Score[best].o] = A[Score[best].o][Score[best].u] = 1;
            edge_count++;

            free(Score);
        }

        if (is_integral_graph(A, v)) {
            licznik_grafow++;
            adjacency_to_graph6(A, v, graph6_output);
            printf("%s\n", graph6_output);
        }

        Graphs[g] = A;
    }

    FILE *fp = fopen(outfile, "w");
    for (int i = 0; i < number_of_trials; i++) {
        adjacency_to_graph6(Graphs[i], v, graph6_output);
        fprintf(fp, "%s\n", graph6_output);
        free(Graphs[i]);
    }

    fclose(fp);
    free(Graphs);

    printf("Liczba wykopanych grafow calkowitoliczbowych: %d\n", licznik_grafow);
    return 0;
}
