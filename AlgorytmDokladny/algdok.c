#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "cJSON.h"

// Algorytm dokładny (bruteforce + backtracking) dla grafu dwudzielnego
// Buduje macierz B (n1 x n2) z wartościami 0/1, a następnie A = [0 B; B^T 0]
// Oblicza wartości własne macierzy A (algorytm Jacobiego) i porównuje z docelowym spektrum

// Prosty parser JSON ograniczony do formatu wygenerowanego przez Generator/graph_to_json.c

double* read_target_spectrum(const char* path, int* out_n) {
    FILE* f = fopen(path, "r");
    if (!f) { perror("fopen target_spectrum"); return NULL; }
    // Najpierw policz linie/liczby
    int capacity = 16;
    double* arr = malloc(sizeof(double)*capacity);
    int n = 0;
    double x;
    while (fscanf(f, "%lf", &x) == 1) {
        if (n >= capacity) { capacity *= 2; arr = realloc(arr, sizeof(double)*capacity); }
        arr[n++] = x;
    }
    fclose(f);
    *out_n = n;
    return arr;
}

// Minimalny JSON "parser" dla pól n1, n2. Nie jest pełny — oczekuje prostego formatu.
// Read graph JSON using cJSON. Also optionally read "edges" array and fill fixed_mask (allocated here).
// Returns 0 on success. On success *fixed_mask_out is allocated (size n1*n2) and filled with 0/1 for fixed edges.
int read_json_graph(const char* path, int* n1, int* n2, unsigned char** fixed_mask_out) {
    FILE* f = fopen(path, "r");
    if (!f) { perror("fopen json"); return 1; }
    fseek(f, 0, SEEK_END);
    long sz = ftell(f);
    fseek(f, 0, SEEK_SET);
    char *buf = malloc(sz+1);
    if (!buf) { fclose(f); return 2; }
    size_t r = fread(buf,1,sz,f);
    buf[r]=0;
    fclose(f);

    cJSON *root = cJSON_Parse(buf);
    free(buf);
    if (!root) return 3;

    cJSON *jn1 = cJSON_GetObjectItemCaseSensitive(root, "n1");
    cJSON *jn2 = cJSON_GetObjectItemCaseSensitive(root, "n2");
    if (!cJSON_IsNumber(jn1) || !cJSON_IsNumber(jn2)) { cJSON_Delete(root); return 4; }
    *n1 = jn1->valueint;
    *n2 = jn2->valueint;

    int total = (*n1) * (*n2);
    unsigned char *fixed = calloc(total, 1);
    if (!fixed) { cJSON_Delete(root); return 5; }

    cJSON *edges = cJSON_GetObjectItemCaseSensitive(root, "edges");
    if (edges && cJSON_IsArray(edges)) {
        cJSON *elem = NULL;
        cJSON_ArrayForEach(elem, edges) {
            if (cJSON_IsArray(elem) && cJSON_GetArraySize(elem) >= 2) {
                cJSON *ja = cJSON_GetArrayItem(elem, 0);
                cJSON *jb = cJSON_GetArrayItem(elem, 1);
                if (cJSON_IsNumber(ja) && cJSON_IsNumber(jb)) {
                    int u = ja->valueint;
                    int v = jb->valueint;
                    if (u >= 0 && u < *n1 && v >= 0 && v < *n2) {
                        fixed[u * (*n2) + v] = 1;
                    }
                }
            }
        }
    }

    *fixed_mask_out = fixed;
    cJSON_Delete(root);
    return 0;
}

// Jacobi method for symmetric matrix eigenvalues (returns array of size n)
// Simple implementation, OK for small n (<= 60)
void jacobi_eigenvalues(double* a, int n, double* eig) {
    // copy a to working matrix
    double *A = malloc(sizeof(double)*n*n);
    memcpy(A, a, sizeof(double)*n*n);
    double *V = calloc(n*n, sizeof(double));
    for (int i=0;i<n;i++) V[i*n+i]=1.0;

    const int MAX_IT = 100*n*n;
    for (int it=0; it<MAX_IT; ++it) {
        // find largest off-diagonal
        int p=0,q=1;
        double max = fabs(A[p*n+q]);
        for (int i=0;i<n;i++) for (int j=i+1;j<n;j++) {
            double val = fabs(A[i*n+j]);
            if (val > max) { max = val; p=i; q=j; }
        }
        if (max < 1e-12) break;
        double app = A[p*n+p];
        double aqq = A[q*n+q];
        double apq = A[p*n+q];
        double phi = 0.5 * atan2(2*apq, aqq - app);
        double c = cos(phi), s = sin(phi);
        // rotate
        for (int i=0;i<n;i++) {
            if (i==p || i==q) continue;
            double aip = A[i*n+p];
            double aiq = A[i*n+q];
            A[i*n+p] = A[p*n+i] = c*aip - s*aiq;
            A[i*n+q] = A[q*n+i] = s*aip + c*aiq;
        }
        double new_pp = c*c*app - 2*s*c*apq + s*s*aqq;
        double new_qq = s*s*app + 2*s*c*apq + c*c*aqq;
        A[p*n+p] = new_pp; A[q*n+q] = new_qq; A[p*n+q] = A[q*n+p] = 0.0;
        // update V
        for (int i=0;i<n;i++) {
            double vip = V[i*n+p];
            double viq = V[i*n+q];
            V[i*n+p] = c*vip - s*viq;
            V[i*n+q] = s*vip + c*viq;
        }
    }
    for (int i=0;i<n;i++) eig[i] = A[i*n+i];
    free(A); free(V);
}

int compare_doubles(const void* a, const void* b) {
    double da = *(const double*)a;
    double db = *(const double*)b;
    if (da < db) return -1; if (da > db) return 1; return 0;
}

// Build full adjacency matrix A from B (n1 x n2) into A (n x n)
void build_full_A(int n1, int n2, const unsigned char* B, double* A) {
    int n = n1 + n2;
    for (int i=0;i<n;i++) for (int j=0;j<n;j++) A[i*n+j]=0.0;
    for (int i=0;i<n1;i++) for (int j=0;j<n2;j++) if (B[i*n2 + j]) {
        int u = i;
        int v = n1 + j;
        A[u*n + v] = 1.0;
        A[v*n + u] = 1.0;
    }
}

// Backtracking search
int found = 0;
// fixed_mask: if fixed_mask[idx]==1 then edge is fixed (must be present)
void search_bt(int n1, int n2, unsigned char* B, const unsigned char* fixed_mask,
               int idx, int total_edges_set, int remaining, double target_sum_sq,
               double* target_eigs, int target_n, double eps) {
    if (found) return;
    if (idx == n1*n2) {
        // complete
        int n = n1 + n2;
        double* A = malloc(sizeof(double)*n*n);
        build_full_A(n1,n2,B,A);
        double* eig = malloc(sizeof(double)*n);
        jacobi_eigenvalues(A,n,eig);
        // copy and sort locally to avoid modifying original target
        double* target_copy = malloc(sizeof(double)*n);
        for (int i=0;i<n;i++) target_copy[i] = target_eigs[i];
        qsort(eig,n,sizeof(double),compare_doubles);
        qsort(target_copy,n,sizeof(double),compare_doubles);
        int ok = 1;
        for (int i=0;i<n;i++) {
            if (fabs(eig[i]-target_copy[i]) > eps) { ok = 0; break; }
        }
        if (ok) {
            printf("FOUND graph with %d x %d bipartition. Edges:\n", n1, n2);
            for (int i=0;i<n1;i++) for (int j=0;j<n2;j++) if (B[i*n2+j]) printf("%d %d\n", i, j);
            found = 1;
        }
        free(A); free(eig);
        free(target_copy);
        return;
    }
    // pruning by sum of squares invariant: sum lambda^2 = trace(A^2) = 2*m
    // target_sum_sq is sum target_eigs^2
    // current edges = total_edges_set, remaining edges possible = remaining
    double curr_min_m = total_edges_set;
    double curr_max_m = total_edges_set + remaining;
    double curr_min_sum_sq = 2.0*curr_min_m;
    double curr_max_sum_sq = 2.0*curr_max_m;
    if (target_sum_sq + 1e-9 < curr_min_sum_sq || target_sum_sq - 1e-9 > curr_max_sum_sq) {
        return;
    }

    int rem_after = remaining - 1;
    if (fixed_mask && fixed_mask[idx]) {
        // forced to 1
        B[idx] = 1;
        search_bt(n1,n2,B,fixed_mask,idx+1,total_edges_set+1,rem_after,target_sum_sq,target_eigs,target_n,eps);
        // do not unset fixed (keep as is)
        return;
    }
    // try 0
    B[idx] = 0;
    search_bt(n1,n2,B,fixed_mask,idx+1,total_edges_set,rem_after,target_sum_sq,target_eigs,target_n,eps);
    if (found) return;
    // try 1
    B[idx] = 1;
    search_bt(n1,n2,B,fixed_mask,idx+1,total_edges_set+1,rem_after,target_sum_sq,target_eigs,target_n,eps);
    B[idx] = 0;
}

int main(int argc, char** argv) {
    if (argc < 3) {
        fprintf(stderr, "Usage: %s graph.json target_spectrum.txt [eps]\n", argv[0]);
        return 1;
    }
    const char* json_path = argv[1];
    const char* spectrum_path = argv[2];
    double eps = 1e-6;
    if (argc >= 4) eps = atof(argv[3]);

    int n1=0,n2=0;
    unsigned char *fixed_mask = NULL;
    if (read_json_graph(json_path, &n1, &n2, &fixed_mask) != 0) {
        fprintf(stderr, "Failed to read graph JSON %s\n", json_path);
        return 2;
    }
    int target_n;
    double* target = read_target_spectrum(spectrum_path, &target_n);
    if (!target) return 3;
    int n = n1 + n2;
    if (target_n != n) {
        fprintf(stderr, "Target spectrum length (%d) does not match n1+n2 (%d)\n", target_n, n);
        free(target); return 4;
    }

    double target_sum_sq = 0.0;
    for (int i=0;i<n;i++) target_sum_sq += target[i]*target[i];

    unsigned char* B = calloc(n1*n2,1);
    found = 0;
    // initialize B for fixed edges (they will be forced in search)
    int fixed_count = 0;
    if (fixed_mask) {
        for (int i=0;i<n1*n2;i++) if (fixed_mask[i]) { B[i]=1; fixed_count++; }
    }
    search_bt(n1,n2,B,fixed_mask,0,fixed_count,n1*n2 - 0,target_sum_sq,target,target_n,eps);
    if (!found) {
        printf("No matching bipartite graph found (for given sizes and spectrum)\n");
    }
    free(B); free(target); if (fixed_mask) free(fixed_mask);
    return 0;
}
