#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "cJSON.h"

// read JSON similar to algdok's reader but returns list of edges
int read_json_edges(const char* path, int* n1, int* n2, unsigned char** B_out) {
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
    unsigned char *B = calloc(total, 1);
    if (!B) { cJSON_Delete(root); return 5; }

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
                        B[u * (*n2) + v] = 1;
                    }
                }
            }
        }
    }

    *B_out = B;
    cJSON_Delete(root);
    return 0;
}

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

int compare_doubles(const void* a, const void* b) {
    double da = *(const double*)a;
    double db = *(const double*)b;
    if (da < db) return -1; if (da > db) return 1; return 0;
}

// Jacobi copied from algdok
void jacobi_eigenvalues(double* a, int n, double* eig) {
    double *A = malloc(sizeof(double)*n*n);
    memcpy(A, a, sizeof(double)*n*n);
    double *V = calloc(n*n, sizeof(double));
    for (int i=0;i<n;i++) V[i*n+i]=1.0;
    const int MAX_IT = 100*n*n;
    for (int it=0; it<MAX_IT; ++it) {
        int p=0,q=1; double max = fabs(A[p*n+q]);
        for (int i=0;i<n;i++) for (int j=i+1;j<n;j++) {
            double val = fabs(A[i*n+j]); if (val > max) { max = val; p=i; q=j; }
        }
        if (max < 1e-12) break;
        double app = A[p*n+p]; double aqq = A[q*n+q]; double apq = A[p*n+q];
        double phi = 0.5 * atan2(2*apq, aqq - app);
        double c = cos(phi), s = sin(phi);
        for (int i=0;i<n;i++) {
            if (i==p || i==q) continue;
            double aip = A[i*n+p]; double aiq = A[i*n+q];
            A[i*n+p] = A[p*n+i] = c*aip - s*aiq;
            A[i*n+q] = A[q*n+i] = s*aip + c*aiq;
        }
        double new_pp = c*c*app - 2*s*c*apq + s*s*aqq;
        double new_qq = s*s*app + 2*s*c*apq + c*c*aqq;
        A[p*n+p] = new_pp; A[q*n+q] = new_qq; A[p*n+q] = A[q*n+p] = 0.0;
        for (int i=0;i<n;i++) {
            double vip = V[i*n+p]; double viq = V[i*n+q];
            V[i*n+p] = c*vip - s*viq; V[i*n+q] = s*vip + c*viq;
        }
    }
    for (int i=0;i<n;i++) eig[i] = A[i*n+i];
    free(A); free(V);
}

int main(int argc, char** argv) {
    if (argc < 2) { fprintf(stderr, "Usage: %s graph.json\n", argv[0]); return 1; }
    const char* json_path = argv[1];
    int n1,n2;
    unsigned char *B = NULL;
    if (read_json_edges(json_path,&n1,&n2,&B)!=0) { fprintf(stderr,"Failed to read JSON\n"); return 2; }
    int n = n1 + n2;
    double *A = malloc(sizeof(double)*n*n);
    build_full_A(n1,n2,B,A);
    double *eig = malloc(sizeof(double)*n);
    jacobi_eigenvalues(A,n,eig);
    qsort(eig,n,sizeof(double),compare_doubles);
    for (int i=0;i<n;i++) printf("%1.12g\n", eig[i]);
    free(A); free(eig); free(B);
    return 0;
}
