Zawarte pliki .txt to są macierze użyte w sprawozdaniu.
aby uruchomić te programy należy:
-dla kodu C:
  wpisać do cmd z pobranym gcc następującą komende aby skompilować program "gcc Program.c -O3 -march=native -ffast-math -o Nazwa".
  Następnie uruchomić program komendą "Nazwa.exe Liczba RodzajLiczba.txt".
-dla kodu Python:
  Należy, posiadając zainstalowany kompilator pythona, typu anaconda, uruchomić komendę "python Program.py Liczba RodzajLiczba.txt".
-dla kodu MatLab:
W wbudowanym w program MatLab terminalu wpisać komendę "eigen(load("RodzajLiczba.txt"))",

Wszystkie powyższe komendy muszą być uruchamiane z tego samego folderu do którego zostaną przekopiowane programy.
Wyjaśnienie pojęć z komend:
-Nazwa: nazwa programu C, można podać dowolną która zawiera zwykłe znaki, np. "eigen","cholesky".
-Liczba: tu należy wpisać rozmiar macierzy, którą chce się podać do programu,np. 3 lub 100.
-Rodzaj: tu należy wpisać jaką macierz zamierza się użyć w programie, np. "rand", "hilb".

Nazwe pliku można podmienić na dowolną, w której będzie zawierała się macierz wpisana w następujacym formacie:
x,x,x
x,x,x
x,x,x
gdzie x oznacza dowolną cyfrę.
