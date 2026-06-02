import numpy as np
import time
import sys
import os

def qr_decomp(A):
    n=A.shape[0]
    Q=A.copy().astype(float)
    R=np.zeros((n,n))
    for j in range(n):
        norm=np.linalg.norm(Q[:,j])
        if norm>1e-12:
            R[j,j]=norm
            Q[:,j]/=norm
        else:
            R[j,j]=0.0
            for k in range(n):
                Q[:,j]=0.0
                Q[k,j]=1.0
                for p in range(j):
                    dot=np.dot(Q[:,p],Q[:,j])
                    Q[:,j]-=dot*Q[:,p]
                f_norm=np.linalg.norm(Q[:,j])
                if f_norm>1e-12:
                    Q[:,j]/=f_norm
                    break
        for k in range(j+1,n):
            r_jk=np.dot(Q[:,j],Q[:,k])
            R[j,k]=r_jk
            Q[:,k]-=r_jk*Q[:,j]
    return Q,R

def qr_method(E,EV):
    n=E.shape[0]
    for k in range(50):
        Q,R=qr_decomp(E)
        E=R@Q
        EV=EV@Q
        if np.max(np.abs(np.tril(E,-1)))<1e-9:
            print(f"Zbieznosc osiagnieta po {k+1} iteracjach.")
            break
    return E,EV

def schur_to_eigenvectors(E,EV):
    n=E.shape[0]
    Y=np.zeros((n,n))
    for j in range(n):
        lambda_val=E[j,j]
        Y[j,j]=1.0
        for i in range(j-1,-1,-1):
            s=0.0
            for k in range(i+1,j+1):
                s+=E[i,k]*Y[k,j]
            diff=lambda_val-E[i,i]
            if abs(diff)<1e-12:
                if diff>=0:
                    diff=1e-12
                else:
                    diff=-1e-12
            Y[i,j]=s/diff
        norm=np.linalg.norm(Y[:j+1,j])
        if norm>1e-12:
            Y[:j+1,j]/=norm
    return EV@Y

def save_matrix(filename, matrix):
    np.savetxt(filename,matrix,fmt="%9.5f",delimiter=",")
    print(f"Zapisano macierz do pliku: {filename}")

def main():
    if len(sys.argv)>=3:
        n=int(sys.argv[1])
        filename=sys.argv[2]
    else:
        n=int(input("Podaj rozmiar macierzy n:"))
        filename=input("Podaj nazwe pliku z danymi:")
        while not os.path.exists(filename):
            filename = input("Blad otwarcia pliku. Podaj poprawna nazwe:")
    try:
        with open(filename, 'r') as f:
            content = f.read().replace(',',' ')
        A=np.fromstring(content,sep=' ').reshape((n,n))
    except Exception as e:
        print(f"Blad: Za malo danych lub niepoprawny format pliku! Szczegóły:{e}")
        return
    E=A.copy()
    EV=np.eye(n)
    start_time=time.perf_counter_ns()
    E,EV=qr_method(E,EV)
    EV=schur_to_eigenvectors(E,EV)
    try:
        EVR=np.linalg.inv(EV)
    except np.linalg.LinAlgError:
        print("Błąd: Macierz jest osobliwa (nieodwracalna)!")
        return
    D=np.diag(np.diag(E))
    end_time=time.perf_counter_ns()
    duration=(end_time-start_time)/1e9
    save_matrix("EVR_py.txt",EVR)
    save_matrix("EV_py.txt",EV)
    save_matrix("D_py.txt",D)
    print(f"czas wykonania algorytmu: {duration:4.9f} sekund")

if __name__ == "__main__":
    main()
