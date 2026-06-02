import numpy as np
import time
import sys
import os

def cholesky_decomposition(A, n):
    L=np.zeros((n,n))
    for i in range(n):
        for j in range(i+1):
            sum_val=L[i,:j]@L[j,:j]
            if i==j:
                diff=A[i,i]-sum_val
                if diff<=0.0:
                    return None
                L[i,j]=np.sqrt(diff)
            else:
                L[i,j]=(A[i,j]-sum_val)/L[j,j]
    return L

def save_matrix(filename, matrix):
    np.savetxt(filename,matrix,fmt="%9.5f",delimiter=",")
    print(f"Zapisano macierz do pliku: {filename}")

def main():
    if len(sys.argv)>=3:
        n=int(sys.argv[1])
        filename=sys.argv[2]
    else:
        n=int(input("Podaj rozmiar macierzy n: "))
        filename = input("Podaj nazwe pliku z danymi: ")
        while not os.path.exists(filename):
            filename = input("Blad otwarcia pliku. Podaj poprawna nazwe: ")
    try:
        with open(filename,'r') as f:
            content=f.read().replace(',',' ')
        A=np.fromstring(content,sep=' ').reshape((n, n))
    except Exception as e:
        print(f"Blad podczas wczytywania danych: {e}")
        return
    start_time=time.perf_counter_ns()
    L=cholesky_decomposition(A, n)
    if L is not None:
        LT=L.T
        save_matrix("L_py.txt", L)
        save_matrix("LT_py.txt", LT)
    else:
        print("Blad: Macierz nie jest symetryczna dodatnio okreslona! Dekompozycja niemozliwa.")
    end_time=time.perf_counter_ns()
    duration=(end_time-start_time)/1e9
    print(f"czas wykonania algorytmu: {duration:.9f} sekund")

if __name__ == "__main__":
    main()
