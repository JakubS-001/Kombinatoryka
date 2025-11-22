# Generowanie grafów dwudzielnych o danym spektrum 

Wylosować graf, obliczyć spektrum, wynik potraktować jak daną wejściową (mamy gwarancję, że co najmniej jeden taki graf istnieje)
Zobacz: Lab 2
Generowanie wszystkich możliwych wartości (spektra całkowite)
Losowanie z zbioru potencjalnych wartości spektrum
jezyk_minizinc


# KTZ: Uwagi
- Pierwszym etapem realizacji projektu jest przegląd literatury i dostępnych w internecie implementacji problemu w wersji klasycznej.
- Można spróbować uruchomić kody dla przykładowych danych
- Sformułować model matematyczny problemu
- Wskazać literaturę ← Zobacz: problemy
- projektowanie_algorytmow_grafowych
- Przygotować słownik pojęć (tabelka)
- Drugi etap to przygotowanie procedur weryfikujących poprawność i jakość wygenerowanego rozwiązania (np. czy dane kolorowanie jest legalne i ile kolorów wykorzystuje)
- skrypty
- Zobacz: przetwarzanie_wsadowe
- Przygotować mały przykład ilustrujący problem i jego rozwiązanie (rysunek).
- Kody najlepiej przygotować w języku C (ułatwi to późniejsze przekształcenie projektu na wersję równoległą, zobacz:openmp)
- Prawidłowa obsługa: int64_t
- Jeżeli program jest w języku Python stosować if _ _name_ _ == „_ _main_ _”:
- Format JSON - jako podstawowy sposób przekazywania danych do/z programu
- Uwaga!. Wstawianie plików JSON do Dokuwiki: <file javascript dane.json> … </file>
- Przygotować skrypty testujące czasy wykonania programu
- kepler_project
- Program powinien reagować na parametry wywołania int main(int argc, char* argv[]) { … }
- Unikać wstawiania do kodu nazw plików bez możliwości ich zmiany przez parametry wywołania programu
- Rozwiązanie dokładne może być generowane przez systemy typu: https://www.minizinc.org/software.html
- Niekiedy wymaga to przygotowania dodatkowych programów pośredniczących (zobacz: minizinc_python_json)
- Kody powinny zawierać komentarze zgodne z : doxygen
- Zobacz: (ang. Literate programming)
- Zobacz: (ang. Self-documenting code)
- Fragmenty kodu wygenerowane za pomocą AI należy prawidłowo oznaczyć (nazwa systemu oraz inne informacje pomocne w późniejszej analizie)
- (ang. Text watermarking)
- Algorytm dokładny powinien mieć możliwość wprowadzenia na wejściu rozwiązania bazowego wygenerowanego np. przez algorytm zachłanny
- Zobacz: przetwarzanie_potokowe ← to może komplikować wymianę danych przez format JSON
- Jeżeli w programie korzysta się z generatora liczb pseudolosowych, to inicjalizacja ziarna powinna być wykonana tylko raz i przekazana przez parametr do metody (dodatkowo/opcjonalnie w programie powinna być możliwość ustawienia wartości ziarna przez argument wywołania - co może być pomocne w czasie testów)
- Zobacz: generowanie_liczb_pseudolosowych
- W algorytmach można dodać mechanizm ograniczający maksymalny czas jego działania (raportowany jest wynik cząstkowy i/lub komunikat o błędzie)
- Zobacz: signal.h
- Zobacz: timeout

status KTZ: nie ok