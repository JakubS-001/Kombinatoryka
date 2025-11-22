// Now using cJSON for JSON output. Compile with: -lcjson
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "cJSON.h"

// Prosty konwerter: wczytuje format wejściowy (edge list):
// Pierwsza linia: n1 n2 m
// Następne m linii: u v   (u in [0,n1-1], v in [0,n2-1])
// Wypisuje JSON na stdout z polami: n1, n2, edges: [[u,v],...]

int main(int argc, char** argv) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s input_edge_list.txt\n", argv[0]);
        return 1;
    }
    const char* infile = argv[1];
    FILE* f = fopen(infile, "r");
    if (!f) { perror("fopen"); return 2; }

    int n1, n2, m;
    if (fscanf(f, "%d %d %d", &n1, &n2, &m) != 3) {
        fprintf(stderr, "Invalid input format. Expected: n1 n2 m\n");
        fclose(f);
        return 3;
    }

    cJSON *root = cJSON_CreateObject();
    if (!root) { fprintf(stderr, "cJSON_CreateObject failed\n"); fclose(f); return 4; }
    cJSON_AddNumberToObject(root, "n1", n1);
    cJSON_AddNumberToObject(root, "n2", n2);
    cJSON *edges = cJSON_AddArrayToObject(root, "edges");

    int u, v;
    for (int i = 0; i < m; ++i) {
        if (fscanf(f, "%d %d", &u, &v) != 2) break;
        cJSON *pair = cJSON_CreateIntArray((int[]){u, v}, 2);
        if (!pair) continue;
        cJSON_AddItemToArray(edges, pair);
    }

    char *out = cJSON_Print(root);
    if (out) {
        printf("%s\n", out);
        free(out);
    }
    cJSON_Delete(root);
    fclose(f);
    return 0;
}
