/**
 * @file genetic.c
 * @brief Algorytm genetyczny wyszukiwania grafów całkowitych (integral graphs)
 *
 * Program znajduje spójny graf całkowity o zadanej liczbie
 * wierzchołków i krawędzi, wykorzystując algorytm genetyczny.
 *
 * @par Argumenty linii poleceń:
 * - <b>V</b> – liczba wierzchołków grafu
 * - <b>E</b> – docelowa liczba krawędzi
 * - <b>POP</b> – rozmiar populacji
 * - <b>GEN</b> – maksymalna liczba generacji
 * - <b>MUT_RATE</b> – liczba mutacji na osobnika
 * - <b>ELITES</b> – liczba najlepszych osobników kopiowanych bez zmian
 * - <b>in_file</b> – plik wejściowy z grafami w formacie graph6
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_eigen.h>

#define EPS 1e-6

/**
 * @struct Graph
 * @brief Reprezentuje graf w algorytmie genetycznym
 */
typedef struct {
    int **adj;        /**< macierz sąsiedztwa (0/1) */
    double fitness;   /**< wartość funkcji przystosowania */
} Graph;

/**
 * @brief Alokuje graf o zadanej liczbie wierzchołków
 *
 * @param V liczba wierzchołków
 * @return nowy graf
 */
Graph graph_alloc(int V) {
    Graph g;
    g.fitness = 0.0;
    g.adj = malloc(V * sizeof(int *));
    for (int i = 0; i < V; i++)
        g.adj[i] = calloc(V, sizeof(int));
    return g;
}

/**
 * @brief Zwalnia pamięć zajmowaną przez graf
 *
 * @param g wskaźnik na graf
 * @param V liczba wierzchołków
 */
void graph_free(Graph *g, int V) {
    for (int i = 0; i < V; i++)
        free(g->adj[i]);
    free(g->adj);
}

/**
 * @brief Przeszukiwanie DFS
 *
 * @param v aktualny wierzchołek
 * @param V liczba wierzchołków
 * @param vis tablica odwiedzin
 * @param A macierz sąsiedztwa
 */
void dfs(int v, int V, int *vis, int **A) {
    vis[v] = 1;
    for (int i = 0; i < V; i++)
        if (A[v][i] && !vis[i])
            dfs(i, V, vis, A);
}

/**
 * @brief Sprawdza czy graf jest spójny
 *
 * @param g wskaźnik na graf
 * @param V liczba wierzchołków
 * @return 1 jeśli spójny, 0 w przeciwnym razie
 */
int is_connected(Graph *g, int V) {
    int *vis = calloc(V, sizeof(int));
    dfs(0, V, vis, g->adj);

    for (int i = 0; i < V; i++) {
        if (!vis[i]) {
            free(vis);
            return 0;
        }
    }

    free(vis);
    return 1;
}

/**
 * @brief Zlicza liczbę krawędzi grafu
 *
 * @param g wskaźnik na graf
 * @param V liczba wierzchołków
 * @return liczba krawędzi
 */
int edge_count(Graph *g, int V) {
    int c = 0;
    for (int i = 0; i < V; i++)
        for (int j = i + 1; j < V; j++)
            c += g->adj[i][j];
    return c;
}

/**
 * @brief Funkcja porównująca grafy wg fitness (do qsort)
 */
int cmp_graph(const void *a, const void *b) {
    const Graph *ga = a;
    const Graph *gb = b;
    if (ga->fitness < gb->fitness) return -1;
    if (ga->fitness > gb->fitness) return 1;
    return 0;
}

/**
 * @brief Oblicza wartości własne macierzy sąsiedztwa
 *
 * @param A macierz sąsiedztwa
 * @param v liczba wierzchołków
 * @param eigenvalues tablica wynikowa
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
 * @brief Sprawdza czy graf jest całkowity
 *
 * @param g wskaźnik na graf
 * @param V liczba wierzchołków
 * @return 1 jeśli graf całkowity, 0 w przeciwnym razie
 */
int is_integral(Graph *g, int V) {
    double *eigenvalues = malloc(V * sizeof(double));
    get_eigenvalues(g->adj, V, eigenvalues);

    for (int i = 0; i < V; i++) {
        if (fabs(eigenvalues[i] - round(eigenvalues[i])) > EPS) {
            free(eigenvalues);
            return 0;
        }
    }

    free(eigenvalues);
    return 1;
}

/**
 * @brief Kara za niecałkowite wartości własne
 *
 * @param A macierz sąsiedztwa
 * @param v liczba wierzchołków
 * @param score akumulowana kara
 */
void integral_grade(int **A, int v, double *score) {
    double *eigenvalues = malloc(v * sizeof(double));
    get_eigenvalues(A, v, eigenvalues);

    for (int i = 0; i < v; i++)
    {
        if (fabs(eigenvalues[i] - round(eigenvalues[i])) > EPS)
            *score += 100.0;
        *score += fabs(eigenvalues[i] - round(eigenvalues[i])) * 10.0;
    }
    free(eigenvalues);
}

/**
 * @brief Funkcja przystosowania (minimalizowana)
 *
 * @param g wskaźnik na graf
 * @param V liczba wierzchołków
 * @param E docelowa liczba krawędzi
 * @return wartość fitness
 */
double fitness(Graph *g, int V, int E) {
    double f = 0.0;
    if (!is_connected(g, V)) f += 50000;
    f += abs((double)(edge_count(g, V) - E)) * 5000;
    integral_grade(g->adj, V, &f);
    return f;
}

/**
 * @brief Selekcja turniejowa
 *
 * @param pop populacja
 * @param POP liczba osobników
 * @return wybrany osobnik
 */
Graph tournament(Graph *pop, int POP) {
    int a = rand() % POP;
    int b = rand() % POP;
    return pop[a].fitness < pop[b].fitness ? pop[a] : pop[b];
}


/**
 * @brief Mutacja grafu
 *
 * @param g wskaźnik na graf
 * @param V liczba wierzchołków
 * @param MUT_RATE parametr mutacji
 */
void mutate(Graph *g, int V, int MUT_RATE) {
    for(int i=0; i<MUT_RATE; i++) {
        int i = rand() % V;
        int j = rand() % V;
        if (i != j)
            g->adj[i][j] = g->adj[j][i] ^= 1;
    }
}

/**
 * @brief Konwersja macierzy sąsiedztwa do graph6
 */
void adjacency_to_graph6(int **A, int n, char *output) {
    int bit_count = 0, value = 0, out_index = 0;
    output[out_index++] = (char)(n + 63);

    for (int i = 0; i < n; i++)
        for (int j = 0; j < i; j++) {
            value = (value << 1) | A[i][j];
            if (++bit_count == 6) {
                output[out_index++] = (char)(value + 63);
                bit_count = value = 0;
            }
        }

    if (bit_count)
        output[out_index++] = (char)((value << (6 - bit_count)) + 63);

    output[out_index] = '\0';
}

/**
 * @brief Konwersja graph6 do macierzy sąsiedztwa
 */
void graph6_to_adjacency(const char *input, int **A, int *n) {
    *n = input[0] - 63;
    int bit_index = 0, byte_index = 1;

    for (int i = 0; i < *n; i++)
        for (int j = 0; j < i; j++) {
            if (bit_index == 6) bit_index = 0, byte_index++;
            int bit = (input[byte_index] - 63) >> (5 - bit_index++) & 1;
            A[i][j] = A[j][i] = bit;
        }
}

/**
 * @brief Inicjalizuje populację z pliku graph6
 */
void initial_population(Graph *pop, int POP, int V, int E, char *in_file) {
    FILE *f = fopen(in_file, "r");
    char line[1000];
    int count = 0;

    while (fgets(line, sizeof(line), f) && count < POP) {
        pop[count] = graph_alloc(V);
        graph6_to_adjacency(line, pop[count].adj, &V);
        pop[count].fitness = fitness(&pop[count], V, E);
        count++;
    }

    fclose(f);
}

/**
 * @brief Główna funkcja programu
 */
int main(int argc, char **argv) {
    if (argc != 8) {
        fprintf(stderr,
            "Użycie: %s V E POP GEN MUT_RATE ELITES in_file\n", argv[0]);
        return 1;
    }

    int V = atoi(argv[1]);
    int E = atoi(argv[2]);
    int POP = atoi(argv[3]);
    int GEN = atoi(argv[4]);
    int MUT_RATE = atoi(argv[5]);
    int ELITES = atoi(argv[6]);
    char *in_file = argv[7];

    srand(time(NULL));

    Graph *pop = malloc(POP * sizeof(Graph));
    Graph *newpop = malloc(POP * sizeof(Graph));

    for (int i = 0; i < POP; i++) {
        pop[i] = graph_alloc(V);
        newpop[i] = graph_alloc(V);
    }

    initial_population(pop, POP, V, E, in_file);
    for (int g = 0; g < GEN; g++) {
        qsort(pop, POP, sizeof(Graph), cmp_graph);

        for (int e = 0; e < ELITES; e++) {
            for (int i = 0; i < V; i++)
                memcpy(newpop[e].adj[i], pop[e].adj[i], V * sizeof(int));
            newpop[e].fitness = pop[e].fitness;
        }

        for (int i = ELITES; i < POP; i++) {
            Graph g = tournament(pop, POP);
            for (int k = 0; k < V; k++)
                memcpy(newpop[i].adj[k], g.adj[k], sizeof(int) * V);
            newpop[i].fitness = g.fitness;
            mutate(&newpop[i], V, MUT_RATE);
            newpop[i].fitness = fitness(&newpop[i], V, E);
        }

        Graph *tmp = pop; 
        pop = newpop; 
        newpop = tmp;

        if (g % 200 == 0)
            fprintf(stderr, "Gen %d | best fitness = %.7f\n", g, pop[0].fitness);

        if (is_integral(&pop[0], V)) {
            char out[1000];
            adjacency_to_graph6(pop[0].adj, V, out);
            printf("%s\n", out);
            fprintf(stderr, "Gen %d | best fitness = %f\n", g, pop[0].fitness);
            break;
        }
    }

    free(pop);
    free(newpop);
    return 0;
}
