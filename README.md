
<p>Zawarte pliki .txt to są macierze użyte w sprawozdaniu.<br>
aby uruchomić te programy należy:<br>
<ul>
<li>dla kodu C:</li>
<ul>
<li>wpisać do cmd z pobranym gcc następującą komende aby skompilować program "gcc Program.c -O3 -march=native -ffast-math -o Nazwa".</li>
<li>Następnie uruchomić program komendą "Nazwa.exe Liczba RodzajLiczba.txt".</li>
</ul>
<li>dla kodu Python:</li>
<ul>
<li>Należy, posiadając zainstalowany kompilator pythona, typu anaconda, uruchomić komendę "python Program.py Liczba RodzajLiczba.txt".</li>
</ul>
<li>dla kodu MatLab:</li>
<ul>
<li>W wbudowanym w program MatLab terminalu wpisać komendę "eigen(load("RodzajLiczba.txt"))".</li>
</ul>
</ul>
</p>
<p></p>Wszystkie powyższe komendy muszą być uruchamiane z tego samego folderu do którego zostaną przekopiowane programy.<br>
Wyjaśnienie pojęć z komend:<br>
<ul>
<li>Program: nazwa programu, jakiego chce się użyć. dostępne są "eigen" i "cholesky".</li>
<li>Nazwa: nazwa programu C, można podać dowolną która zawiera zwykłe znaki, np. "eigen","cholesky".</li>
<li>Liczba: tu należy wpisać rozmiar macierzy, którą chce się podać do programu,np. 3 lub 100.</li>
<li>Rodzaj: tu należy wpisać jaką macierz zamierza się użyć w programie, np. "rand", "hilb".</li>
</ul>
</p>
<p></p>Nazwe pliku można podmienić na dowolną, w której będzie zawierała się macierz wpisana w następujacym formacie:<br>
x,x,x<br>
x,x,x<br>
x,x,x<br>
gdzie x oznacza dowolną cyfrę.<br>
</p>
