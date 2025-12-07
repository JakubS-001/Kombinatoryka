#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

/* LAPACK prototype */
extern void dsyev_(char *jobz, char *uplo, int *n,
                   double *a, int *lda, double *w,
                   double *work, int *lwork, int *info);

double tol = 1e-6;

/* ---------- read spectrum ---------- */
double* read_spectrum(const char *path, int *n_out) {
    FILE *f = fopen(path, "r");
    if (!f) { printf("Cannot open spectrum file\n"); exit(1); }

    double *arr = malloc(256 * sizeof(double));
    int cap = 256, n = 0;

    while (fscanf(f, "%lf", &arr[n]) == 1) {
        n++;
        if (n == cap) {
            cap *= 2;
            arr = realloc(arr, cap * sizeof(double));
        }
    }
    fclose(f);
    *n_out = n;
    return arr;
}

/* ---------- check spectrum symmetry ---------- */
int is_symmetric_spectrum(double *spec, int n) {
    for (int i = 0; i < n/2; i++)
        if (fabs(spec[i] + spec[n-1-i]) > 1e-9)
            return 0;
    return 1;
}

/* ---------- sort array ---------- */
int cmp_double(const void *a, const void *b) {
    double x = *(double*)a, y = *(double*)b;
    if (x < y) return -1;
    if (x > y) return 1;
    return 0;
}

/* ---------- eigenvalues using LAPACK ---------- */
void compute_eigenvalues(double *A, int n, double *eig) {
    int info, lwork = -1;
    double wkopt;
    double *work;

    /* query optimal workspace */
    dsyev_("N", "U", &n, A, &n, eig, &wkopt, &lwork, &info);
    lwork = (int)wkopt;
    work = (double*)malloc(lwork * sizeof(double));

    /* compute eigenvalues */
    dsyev_("N", "U", &n, A, &n, eig, work, &lwork, &info);

    free(work);
}

/* ---------- adjacency matrix ---------- */
void adjacency_from_edges(int n, int edgeCount, int edges[][2], double *A) {
    for (int i = 0; i < n*n; i++)
        A[i] = 0;

    for (int i = 0; i < edgeCount; i++) {
        int u = edges[i][0];
        int v = edges[i][1];
        A[u*n + v] = 1;
        A[v*n + u] = 1;
    }
}

/* ---------- all bipartitions ---------- */
int next_subset(int *subset, int n) {
    for (int i = 0; i < n; i++) {
        if (subset[i] == 0) {
            subset[i] = 1;
            for (int j = 0; j < i; j++) subset[j] = 0;
            return 1;
        }
    }
    return 0;
}

/* ---------- main search ---------- */
void search_graph(double *spec, int n, const char *output_path) {
    qsort(spec, n, sizeof(double), cmp_double);

    if (!is_symmetric_spectrum(spec, n)) {
        FILE *out = fopen(output_path, "w");
        fprintf(out, "{ \"result\": null }\n");
        fclose(out);
        return;
    }

    int *subset = calloc(n, sizeof(int));
    int *Aset = malloc(n * sizeof(int));
    int *Bset = malloc(n * sizeof(int));

    double *matrix = malloc(n*n * sizeof(double));
    double *eig = malloc(n * sizeof(double));

    /* iterate partitions */
    do {
        int Ac = 0, Bc = 0;

        for (int i = 0; i < n; i++) {
            if (subset[i]) Aset[Ac++] = i;
            else Bset[Bc++] = i;
        }
        if (Ac == 0 || Bc == 0) continue;

        /* possible edges */
        int ABcount = Ac * Bc;
        int (*AB)[2] = malloc(ABcount * sizeof *AB);

        int idx = 0;
        for (int i = 0; i < Ac; i++)
            for (int j = 0; j < Bc; j++)
                AB[idx][0] = Aset[i],
                AB[idx][1] = Bset[j],
                idx++;

        /* iterate subgraphs */
        long long totalMasks = 1LL << ABcount;

        for (long long mask = 0; mask < totalMasks; mask++) {
            int edges[256][2];
            int ec = 0;

            for (int i = 0; i < ABcount; i++)
                if (mask & (1LL << i))
                    edges[ec][0] = AB[i][0],
                    edges[ec][1] = AB[i][1],
                    ec++;

            adjacency_from_edges(n, ec, edges, matrix);

            double *temp = malloc(n*n*sizeof(double));
            memcpy(temp, matrix, n*n*sizeof(double));

            compute_eigenvalues(temp, n, eig);
            free(temp);

            qsort(eig, n, sizeof(double), cmp_double);

            int ok = 1;
            for (int i = 0; i < n; i++)
                if (fabs(eig[i] - spec[i]) > tol) { ok = 0; break; }

            if (ok) {
                FILE *out = fopen(output_path, "w");
                fprintf(out, "{\n");
                fprintf(out, "  \"result\": {\n");
                fprintf(out, "    \"n\": %d,\n", n);

                fprintf(out, "    \"parts\": [ [");
                for (int i = 0; i < Ac; i++)
                    fprintf(out, "%d%s", Aset[i], i+1<Ac?", ":"");
                fprintf(out, "], [");
                for (int i = 0; i < Bc; i++)
                    fprintf(out, "%d%s", Bset[i], i+1<Bc?", ":"");
                fprintf(out, "] ],\n");

                fprintf(out, "    \"edges\": [");
                for (int i = 0; i < ec; i++)
                    fprintf(out, "[%d,%d]%s",
                            edges[i][0], edges[i][1],
                            i+1<ec?", ":"");
                fprintf(out, "],\n");

                fprintf(out, "    \"adjacency_matrix\": [\n");
                for (int i = 0; i < n; i++) {
                    fprintf(out, "      [");
                    for (int j = 0; j < n; j++)
                        fprintf(out, "%.0f%s", matrix[i*n+j], j+1<n?", ":"");
                    fprintf(out, "]%s\n", i+1<n?",":"");
                }
                fprintf(out, "    ]\n");

                fprintf(out, "  }\n");
                fprintf(out, "}\n");
                fclose(out);

                exit(0);
            }
        }

        free(AB);
    } while (next_subset(subset, n));

    /* nothing found */
    FILE *out = fopen(output_path, "w");
    fprintf(out, "{ \"result\": null }\n");
    fclose(out);
}

int main(int argc, char **argv) {
    if (argc != 3) {
        printf("Usage:\n  %s spectrum.txt output.json\n", argv[0]);
        return 1;
    }

    int n;
    double *spec = read_spectrum(argv[1], &n);

    search_graph(spec, n, argv[2]);
    return 0;
}
