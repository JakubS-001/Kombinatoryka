import json
import sys
import itertools
import numpy as np

def is_symmetric_spectrum(spec, tol=1e-9):
    spec_sorted = sorted(spec)
    return all(abs(spec_sorted[i] + spec_sorted[-1-i]) < tol
               for i in range(len(spec)//2))

def bipartitions(n):
    """Zwraca wszystkie podziały {0..n-1} na dwie części."""
    V = list(range(n))
    for r in range(1, n):
        for part in itertools.combinations(V, r):
            A = list(part)
            B = list(set(V) - set(A))
            if len(A) > 0 and len(B) > 0:
                yield A, B

def adjacency_from_edges(n, edges):
    M = np.zeros((n, n), dtype=float)
    for u, v in edges:
        M[u, v] = 1
        M[v, u] = 1
    return M

def find_bipartite_graph_with_spectrum(target_spec, tol=1e-6):
    target_spec = list(target_spec)
    target_spec_sorted = sorted(target_spec)

    if not is_symmetric_spectrum(target_spec):
        return None  # spektrum niedwudzielne

    n = len(target_spec)

    for A, B in bipartitions(n):
        AB = list(itertools.product(A, B))
        m = len(AB)

        for mask in range(1 << m):
            edges = []
            for i in range(m):
                if (mask >> i) & 1:
                    edges.append(list(AB[i]))  # JSON-friendly list

            M = adjacency_from_edges(n, edges)

            eigenvalues = np.linalg.eigvalsh(M)

            if np.allclose(sorted(eigenvalues), target_spec_sorted, atol=tol):
                return {
                    "n": n,
                    "parts": [A, B],
                    "edges": edges,
                    "adjacency_matrix": M.tolist()
                }

    return None


def load_spectrum(path):
    """Ładuje spektrum z pliku .json lub .txt (lista liczb)."""
    if path.endswith(".json"):
        with open(path) as f:
            return json.load(f)

    # .txt → jedna lista liczb, rozdzielona spacjami lub enterami
    with open(path) as f:
        data = f.read().strip().replace(",", " ").split()
        return [float(x) for x in data]


def main():
    if len(sys.argv) != 3:
        print("Użycie: python find_bipartite_from_spectrum.py input_spectrum.json output.json")
        sys.exit(1)

    input_path = sys.argv[1]
    output_path = sys.argv[2]

    spectrum = load_spectrum(input_path)
    result = find_bipartite_graph_with_spectrum(spectrum)

    with open(output_path, "w") as f:
        json.dump({"result": result}, f, indent=2)

    print(f"Zapisano wynik do {output_path}")


if __name__ == "__main__":
    main()
