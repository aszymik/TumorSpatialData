# Analiza TLS w próbkach biopsji nowotworów

## Opis projektu

Ten projekt skupia się na analizie i wizualizacji TLS – trzeciorzędowych struktur limfatycznych (ang. tertiary lymphoid structures), które są obecne w próbkach biopsji od pacjentów z nowotworami. 

Wyniki przedstawione są w aplikacji Streamlit, w której użytkownik może wybrać konkretnego pacjenta i zobaczyć wizualizację fragmentu jego tkanki, z wyszczególnionymi typami komórek, oraz zidentyfikowanych u niego kandydatów na TLS. 

W aplikacji przedstawione są dodatkowo różne statystyki dotyczące zarówno TLS danego pacjenta, jak i zbiorczych danych zebranych od wszystkich badanych pacjentów. Dzięki temu możliwe jest porównywanie wyników i lepsze zrozumienie kontekstu badań oraz natury i różnorodności struktur TLS.

![app screen](graphics/app_screen.png)
*Wygląd aplikacji*

## Struktura projektu

Projekt składa się z trzech głównych skryptów: 

* `main.py`: Zbiera i analizuje dane dotyczące liczby i składu TLS u wszystkich pacjentów i zapisuje je do pliku, aby umożliwić łatwy dostęp do danych.
* `patient_statistics.py`: Zbiera statystyki i wizualizuje dane od poszczególnych pacjentów.
* `app.py`: Podsumowuje wyniki całego projektu, prezentując je w aplikacji Streamlit.

## Wymagania

Wszystkie wymagane pakiety zapisane są w pliku ```requirements.txt``` i mogą być zainstalowane komendą ```pip install -r requirements.txt```.

## Uruchamianie

Po przejściu do katalogu źródłowego projektu, można uruchomić aplikację następującą komendą:

```streamlit run app.py```


Projekt stanowi zadanie zaliczeniowe z Biologii Systemów 2023/2024.
