# Optymalizacja Kombinatoryczna

## Krótki wstęp

Projekt służy do poszukiwania grafów całkowitych (integral graphs).
Zaimplementowano dwie strategie konstrukcji: algorytm zachłanny ulosowiony (greedy_random)
i genetyczny (genetic).

### Wymagania

- GCC
- GNU Scientific Library (GSL)
- (opcjonalnie) Doxygen do generowania dokumentacji w `docs`

## Kompilacja

W katalogu głównym uruchom:

```bash
gcc genetic.c -o genetic -lgsl -lgslcblas -lm
gcc greedy_random.c -o greedy_random -lgsl -lgslcblas -lm
```

## Uruchomienie

./greed_random v e greed_number number_of_trials outfile

Argumenty linii poleceń:
 * - <b>v</b> – liczba wierzchołków grafu
 * - <b>e</b> – docelowa liczba krawędzi
 * - <b>greed_number</b> – liczba potencjalnych krawędzi rozważanych na każdym kroku
 * - <b>number_of_trials</b> – liczba wygenerowanych grafów
 * - <b>outfile</b> – plik wyjściowy do zapisu wygenerowanych grafów w formacie graph6
  
  ./genetic v e pop gen mut_rate elites in_file

  Argumenty linii poleceń:
 * - <b>v</b> – liczba wierzchołków grafu
 * - <b>e</b> – docelowa liczba krawędzi
 * - <b>pop</b> – rozmiar populacji
 * - <b>gen</b> – maksymalna liczba generacji
 * - <b>mut_rate</b> – parametr mutacji (mutacja z prawd. 1/MUT_RATE)
 * - <b>elites</b> – liczba najlepszych osobników kopiowanych bez zmian
 * - <b>in_file</b> – plik wejściowy z grafami w formacie graph6
