# Analiza TLS w próbkach biopsji nowotworów

## Opis projektu

Ten projekt skupia się na analizie i wizualizacji TLS (Tertiary Lymphoid Structures) w próbkach biopsji pacjentów z nowotworami.

![app screen](graphics/app_screen.png)

## Struktura projektu

Projekt składa się z trzech głównych skryptów: 

* `main.py`: Zbiera i analizuje dane dotyczące liczby i składu TLS u wszystkich pacjentów i zapisuje je do pliku, aby umożliwić łatwy dostęp do danych.
* `patient_statistics.py`: Zajmuje się wizualizacją danych od poszczególnych pacjentów.
* `app.py`: Zbiera wyniki całego projektu i prezentuje je w aplikacji Streamlit.

## Wymagania

Wszystkie wymagane pakiety zapisane są w pliku ```requirements.txt``` (możesz zainstalować je komendą ```pip install requirements.txt```).

## Jak uruchomić

Projekt może zostać uruchomiony następującą komendą:
```streamlit run app.py```