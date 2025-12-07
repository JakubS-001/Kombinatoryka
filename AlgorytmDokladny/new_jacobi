#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <limits.h>
#include <stdint.h>

#define EPS_SYM 1e-9
#define EPS_EIG 1e-8
/* maksymalna liczba iteracji Jacobiego: ~n^3 * faktor */
#define MAX_ITERS_FACTOR 5

/* ---------- helper: read spectrum ---------- */
double* read_spectrum(const char *path, int *n_out) {
    FILE *f = fopen(path, "r");
    if (!f) {
        fprintf(stderr, "Nie można otworzyć pliku: %s\n", path);
        exit(1);
    }

    int cap = 256;
    double *arr = malloc(cap * sizeof(double));
    if (!arr) { fprintf(stderr, "Brak pamięci\n"); exit(1); }

    int n = 0;
    while (fscanf(f, "%lf", &arr[n]) == 1) {
        n++;
        if (n == cap) {
            cap *= 2;
            arr = realloc(arr, cap * sizeof(double));
            if (!arr) { fprintf(stderr, "Brak pamięci\n"); exit(1); }
        }
    }
    fclose(f);
    *n_out = n;
    return arr;
}

/* ---------- comparator for qsort ---------- */
int cmp_double(const void *a, const void *b) {
    double x = *(const double*)a;
    double y = *(const double*)b;
    if (x < y) return -1;
    if (x > y) return 1;
    return 0;
}

/* ---------- check symmetric spectrum (necessary for bipartite) ---------- */
int is_symmetric_spectrum(double *spec, int n) {
    double *tmp = malloc(n * sizeof(double));
    if (!tmp) return 0;
    memcpy(tmp, spec, n * sizeof(double));
    qsort(tmp, n, sizeof(double), cmp_double);
    int ok = 1;
    for (int i = 0; i < n/2; ++i) {
        if (fabs(tmp[i] + tmp[n-1-i]) > EPS_SYM) { ok = 0; break; }
    }
    free(tmp);
    return ok;
}

/* ---------- adjacency matrix builder ---------- */
void adjacency_from_edges(int n, int edgeCount, int edges[][2], double *A) {
    for (int i = 0; i < n*n; ++i) A[i] = 0.0;
    for (int e = 0; e < edgeCount; ++e) {
        int u = edges[e][0];
        int v = edges[e][1];
        A[u*n + v] = 1.0;
        A[v*n + u] = 1.0;
    }
}

/* ---------- Jacobi eigenvalue solver (returns eigenvalues in eig array) ---------- */
/* uses max-element pivot selection */
void jacobi_eigen_decomposition(double *A, int n, double *eig) {
    if (n == 0) return;
    if (n == 1) { eig[0] = A[0]; return; }

    /* work on a copy to avoid modifying original */
    double *M = malloc(n * n * sizeof(double));
    if (!M) { fprintf(stderr, "Brak pamięci\n"); exit(1); }
    memcpy(M, A, n*n*sizeof(double));

    /* optional: matrix of eigenvectors (not needed for eigenvalues only) */
    /* double *V = malloc(n*n*sizeof(double));
       for (int i=0;i<n;i++) for (int j=0;j<n;j++) V[i*n+j] = (i==j)?1.0:0.0; */

    int max_iters = MAX_ITERS_FACTOR * n * n * n;
    for (int iter = 0; iter < max_iters; ++iter) {
        /* find largest off-diagonal element */
        int p = 0, q = 1;
        double maxval = 0.0;
        for (int i = 0; i < n; ++i) {
            for (int j = i+1; j < n; ++j) {
                double v = fabs(M[i*n + j]);
                if (v > maxval) { maxval = v; p = i; q = j; }
            }
        }

        if (maxval < EPS_EIG) break; /* enough diagonalized */

        double App = M[p*n + p];
        double Aqq = M[q*n + q];
        double Apq = M[p*n + q];

        /* compute rotation angle */
        double phi = 0.5 * atan2(2.0*Apq, Aqq - App);
        double c = cos(phi);
        double s = sin(phi);

        /* rotate: update matrix M */
        /* compute new A_pp, A_qq, set A_pq = 0 */
        double new_App = c*c*App - 2.0*c*s*Apq + s*s*Aqq;
        double new_Aqq = s*s*App + 2.0*c*s*Apq + c*c*Aqq;
        M[p*n + p] = new_App;
        M[q*n + q] = new_Aqq;
        M[p*n + q] = 0.0;
        M[q*n + p] = 0.0;

        /* update other entries */
        for (int r = 0; r < n; ++r) {
            if (r == p || r == q) continue;
            double Arp = M[r*n + p];
            double Arq = M[r*n + q];
            double new_Arp = c*Arp - s*Arq;
            double new_Arq = s*Arp + c*Arq;
            M[r*n + p] = new_Arp;
            M[p*n + r] = new_Arp;
            M[r*n + q] = new_Arq;
            M[q*n + r] = new_Arq;
        }

        /* optionally update eigenvector matrix V if needed
           for (int r=0;r<n;++r) {
               double Vrp = V[r*n + p];
               double Vrq = V[r*n + q];
               V[r*n + p] = c*Vrp - s*Vrq;
               V[r*n + q] = s*Vrp + c*Vrq;
           }
        */
    }

    /* collect diagonal as eigenvalues */
    for (int i = 0; i < n; ++i) eig[i] = M[i*n + i];

    free(M);
    /* free(V); */
}

/* ---------- main search: iterate bipartitions and subgraphs ---------- */
void search_graph(double *spec_in, int n, const char *output_path) {
    /* sort target spectrum copy */
    double *spec = malloc(n * sizeof(double));
    memcpy(spec, spec_in, n * sizeof(double));
    qsort(spec, n, sizeof(double), cmp_double);

    /* symmetry check */
    if (!is_symmetric_spectrum(spec, n)) {
        FILE *out = fopen(output_path, "w");
        if (!out) { fprintf(stderr, "Cannot open output\n"); exit(1); }
        fprintf(out, "{ \"result\": null }\n");
        fclose(out);
        free(spec);
        return;
    }

    /* we'll iterate over nontrivial partitions via bitmask 1..(1<<n)-2 */
    int maxPartMask = (n >= (int)(8*sizeof(unsigned long long))) ? -1 : ((1<<n)-1);

    /* buffer allocations */
    double *matrix = malloc(n*n * sizeof(double));
    double *eig = malloc(n * sizeof(double));
    if (!matrix || !eig) { fprintf(stderr, "Brak pamięci\n"); exit(1); }

    for (int part = 1; part < (1<<n)-1; ++part) {
        /* build A and B lists */
        int Ac = 0, Bc = 0;
        int *Aset = malloc(n * sizeof(int));
        int *Bset = malloc(n * sizeof(int));
        for (int v = 0; v < n; ++v) {
            if ( (part >> v) & 1 ) Aset[Ac++] = v;
            else Bset[Bc++] = v;
        }
        if (Ac == 0 || Bc == 0) { free(Aset); free(Bset); continue; }

        int ABcount = Ac * Bc;
        if (ABcount > 60) {
            /* safety: too many edges to iterate (would overflow 1<<ABcount), skip */
            free(Aset); free(Bset);
            continue;
        }

        /* list all AB edges */
        int (*AB)[2] = malloc(ABcount * sizeof *AB);
        int idx = 0;
        for (int i = 0; i < Ac; ++i) for (int j = 0; j < Bc; ++j)
            { AB[idx][0] = Aset[i]; AB[idx][1] = Bset[j]; idx++; }

        long long totalMasks = 1LL << ABcount;
        for (long long mask = 0; mask < totalMasks; ++mask) {
            /* build edge list */
            int edges_count = 0;
            /* small fixed-size edges array; if ABcount large, edges_count will be <= ABcount */
            int (*edges)[2] = malloc(ABcount * sizeof *edges);
            for (int i = 0; i < ABcount; ++i) {
                if (mask & (1LL << i)) {
                    edges[edges_count][0] = AB[i][0];
                    edges[edges_count][1] = AB[i][1];
                    edges_count++;
                }
            }

            adjacency_from_edges(n, edges_count, edges, matrix);

            /* compute eigenvalues via Jacobi */
            jacobi_eigen_decomposition(matrix, n, eig);

            /* sort and compare with target */
            qsort(eig, n, sizeof(double), cmp_double);

            int match = 1;
            for (int i = 0; i < n; ++i) {
                if (fabs(eig[i] - spec[i]) > 1e-5) { match = 0; break; }
            }

            free(edges);

            if (match) {
                /* write JSON output */
                FILE *out = fopen(output_path, "w");
                if (!out) { fprintf(stderr, "Cannot open output\n"); exit(1); }

                fprintf(out, "{\n");
                fprintf(out, "  \"result\": {\n");
                fprintf(out, "    \"n\": %d,\n", n);

                /* write parts */
                fprintf(out, "    \"parts\": [ [");
                for (int i = 0; i < Ac; ++i) {
                    fprintf(out, "%d%s", Aset[i], (i+1<Ac)?", ":"");
                }
                fprintf(out, "], [");
                for (int i = 0; i < Bc; ++i) {
                    fprintf(out, "%d%s", Bset[i], (i+1<Bc)?", ":"");
                }
                fprintf(out, "] ],\n");

                /* write edges */
                fprintf(out, "    \"edges\": [");
                int first = 1;
                for (int i = 0; i < ABcount; ++i) {
                    if ( (mask >> i) & 1 ) {
                        if (!first) fprintf(out, ", ");
                        fprintf(out, "[%d, %d]", AB[i][0], AB[i][1]);
                        first = 0;
                    }
                }
                fprintf(out, "],\n");

                /* adjacency matrix */
                fprintf(out, "    \"adjacency_matrix\": [\n");
                for (int i = 0; i < n; ++i) {
                    fprintf(out, "      [");
                    for (int j = 0; j < n; ++j) {
                        /* print as integer 0/1 */
                        int val = (int) round(matrix[i*n + j]);
                        fprintf(out, "%d%s", val, (j+1<n)?", ":"");
                    }
                    fprintf(out, "]%s\n", (i+1<n)?",":"");
                }
                fprintf(out, "    ]\n");

                fprintf(out, "  }\n");
                fprintf(out, "}\n");
                fclose(out);

                free(Aset); free(Bset); free(AB);
                free(matrix); free(eig); free(spec);
                return;
            }
        }

        free(AB);
        free(Aset); free(Bset);
    }

    /* jeśli nic nie znaleziono */
    FILE *out = fopen(output_path, "w");
    if (!out) { fprintf(stderr, "Cannot open output\n"); exit(1); }
    fprintf(out, "{ \"result\": null }\n");
    fclose(out);

    free(matrix); free(eig); free(spec);
}

/* ---------- main ---------- */
int main(int argc, char **argv) {
    if (argc != 3) {
        fprintf(stderr, "Użycie: %s spectrum.txt output.json\n", argv[0]);
        return 1;
    }

    int n;
    double *spec = read_spectrum(argv[1], &n);
    if (n <= 0) {
        fprintf(stderr, "Brak liczb we spektrum.\n");
        return 1;
    }

    search_graph(spec, n, argv[2]);
    return 0;
}
