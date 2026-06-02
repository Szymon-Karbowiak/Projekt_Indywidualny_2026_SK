#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <windows.h>

long long get_nanoseconds() {
    LARGE_INTEGER freq, counter;
    QueryPerformanceFrequency(&freq);
    QueryPerformanceCounter(&counter);
    return (long long)(counter.QuadPart * 1000000000LL / freq.QuadPart);
}

double* mem_allocate(int n){
    double *A=(double *) malloc(n*n*sizeof(double));
    if(A==NULL){
        exit(1);
    }
    return A;
}

void write_to_file(const char *filename, double *A, int n) {
    FILE *file = fopen(filename, "w");
    if (file == NULL) {
        printf("Blad: Nie udalo sie utworzyc pliku %s\n", filename);
        return;
    }

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            fprintf(file, "%9.5f", A[i*n+j]);
            if (j < n - 1) {
                fprintf(file, ",");
            }
        }
        fprintf(file, "\n");
    }
    fclose(file);
    printf("Zapisano macierz do pliku: %s\n", filename);
}

void multiply_matrices(double *A, double *B, double *C, int n){
    //makes an equation of C=A*B, !!! WILL OVERRITE VALUES IN C !!!
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){
            C[i*n+j]=0;
        }
    }
    for(int i=0;i<n;i++){
        for(int k=0;k<n;k++){
            double temp=A[i*n+k];
            for(int j=0;j<n;j++){
             C[i*n+j]+=temp*B[k*n+j];
            }
        }
    }
    return;
}

void assign_I(double *A, int n){
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){
            if(i==j){
                A[i*n+j]=1;
            }else{
                A[i*n+j]=0;
            }
        }
    }
return;
}

int is_converged(double *A, int n) {
    for (int i=1;i<n;i++) {
        for (int j=0;j<i;j++) {
            if (fabs(A[i*n+j])>1e-9) {
                return 0;
            }
        }
    }
    return 1;
}

void invert_matrix(double *EV, double *EVR, int n) {
    double *Temp=mem_allocate(n);
    for (int i=0;i<n;i++) {
        for (int j=0;j<n;j++) {
            Temp[i*n+j]=EV[i*n+j];
        }
    }
    assign_I(EVR,n);
    for (int i=0;i<n;i++) {
        double pivot=Temp[i*n+i];
        if (fabs(pivot)<1e-18) {
            fprintf(stderr,"Błąd: Macierz jest osobliwa (nieodwracalna)!\n");
            return;
        }
        pivot=1/pivot;
        for (int j=0;j<n;j++) {
            Temp[i*n+j]*=pivot;
            EVR[i*n+j]*=pivot;
        }

        for (int k=0;k<n;k++) {
            if (k!=i) {
                double factor=Temp[k*n+i];
                for (int j=0;j<n;j++) {
                    Temp[k*n+j]-=factor*Temp[i*n+j];
                    EVR[k*n+j]-=factor*EVR[i*n+j];
                }
            }
        }
    }
    free(Temp);
}

void QR_decomp(double *A, double *Q, double *R, int n){
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){
            R[i*n+j]=0.0;
            Q[i*n+j]=A[i*n+j];
        }
    }

    double norm;
    for(int j=0;j<n;j++){
        norm=0.0;
        for (int i=0;i<n;i++) {
            norm+=Q[i*n+j]*Q[i*n+j];
        }
        norm=sqrt(norm);
        if (norm>1e-12) {
            R[j*n+j]=norm;
            norm=1/norm;
            for (int i=0;i<n;i++) {
                Q[i*n+j]*=norm;
            }
        } else {
            R[j*n+j]=0.0;
            for (int k=0;k<n;k++) {
                for (int i=0;i<n;i++) {
                    if(i==k){
                        Q[i*n+j]=1.0;
                    }else{
                        Q[i*n+j]=0.0;
                    }
                }
                for (int p=0;p<j;p++) {
                    double dot=0.0;
                    for (int i=0;i<n;i++) {
                        dot+=Q[i*n+p]*Q[i*n+j];
                    }
                    for (int i=0;i<n;i++) {
                        Q[i*n+j]-=dot*Q[i*n+p];
                    }
                }
                double f_norm=0.0;
                for (int i=0;i<n;i++) {
                    f_norm+=Q[i*n+j]*Q[i*n+j];
                }
                f_norm=sqrt(f_norm);
                if (f_norm>1e-12) {
                    f_norm=1/f_norm;
                    for (int i=0;i<n;i++) {
                        Q[i*n+j]*=f_norm;
                    }
                    break;
                }
            }
        }
        for (int k=j+1;k<n;k++) {
            double r_jk=0.0;
            for (int i=0;i<n;i++) {
                r_jk+=Q[i*n+j]*Q[i*n+k];
            }
            R[j*n+k]=r_jk;

            for (int i=0;i<n;i++) {
                Q[i*n+k]-=R[j*n+k]*Q[i*n+j];
            }
        }
    }
}

void QR_Method(double *E, double *EV, int n)
{
    // E is input array, and at the same time, resulting array with eigenvalues, because it becomes diagonal version of matrix Q,
    //and since you compute eigenvalues by det(A-Lamda*I), diagonal matrix values are its eigenvalues.
    double *Q, *R, *Temp;
    Q=mem_allocate(n);
    R=mem_allocate(n);
    Temp=mem_allocate(n);

    for(int k=0;k<50;k++){
        QR_decomp(E,Q,R,n);
        multiply_matrices(R,Q,E,n);
        multiply_matrices(EV,Q,Temp,n);
        for(int i=0;i<n;i++){
            for(int j=0;j<n;j++){
               EV[i*n+j]=Temp[i*n+j];
            }
        }
        if (is_converged(E,n)) {
            printf("Zbieznosc osiagnieta po %d iteracjach.\n", k + 1);
            break;
        }
    }
    free(Q);
    free(R);
    free(Temp);
}

void Schur_to_Eigenvectors(double *E, double *EV, int n) {
    double *Y = mem_allocate(n);
    double *Temp = mem_allocate(n);

    for (int j=0;j<n;j++) {
        double lambda=E[j*n+j];

        for (int i=0;i<n;i++) {
            Y[i*n+j]=0.0;
        }
        Y[j*n+j]=1.0;

        for (int i=j-1;i>=0;i--) {
            double sum=0.0;
            for (int k=i+1;k<=j;k++) {
                sum+=E[i*n+k]*Y[k*n+j];
            }
            double diff=lambda-E[i*n+i];
            if (fabs(diff)<1e-12) {
                if(diff>=0){
                    diff=1e-12;
                }else{
                    diff=-1e-12;
                }
            }
            Y[i*n+j]=sum/diff;
        }

        double norm=0.0;
        for (int i=0;i<=j;i++) {
            norm+=Y[i*n+j]*Y[i*n+j];
        }
        norm=sqrt(norm);
        if(norm>1e-12) {
            norm=1/norm;
            for(int i=0;i<=j;i++) {
                Y[i*n+j]*=norm;
            }
        }
    }
    multiply_matrices(EV, Y, Temp, n);

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            EV[i*n+j] = Temp[i*n+j];
        }
    }
    free(Y);
    free(Temp);
}

int main(int argc, char **argv)
{
double *A,*E,*EV,*EVR,*D;
int n;
FILE *file = NULL;
if (argc >= 3) {
        n = atoi(argv[1]);
        file = fopen(argv[2], "r");
        if (!file) {
            printf("Nie udało się otworzyć pliku: %s\n", argv[2]);
            return 1;
        }
    } else {
        // Opcja B: Wpisywanie ręczne
        printf("Podaj rozmiar macierzy n: ");
        scanf("%d", &n);

        char filename[128];
        printf("Podaj nazwe pliku z danymi: ");
        scanf("%s", filename);

        file = fopen(filename, "r");
        while (!file) {
            printf("Blad otwarcia pliku. Podaj poprawna nazwe: ");
            scanf("%s", filename);
            file = fopen(filename, "r");
        }
    }
A=mem_allocate(n);
E=mem_allocate(n);
EV=mem_allocate(n);
EVR=mem_allocate(n);
D=mem_allocate(n);
for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
        if (fscanf(file, "%lf", &A[i*n+j]) != 1) {
            printf("Blad: Za malo danych w pliku dla macierzy %dx%d!\n", n, n);
            fclose(file);
            return 1;
        }
        if (j<n-1) {
                char writeme;
                if (fscanf(file, " %c", &writeme) == 1 && writeme != ',') {
                    printf("Ostrzeżenie: Oczekiwano przecinka, znaleziono: %c\n", writeme);
                }
            }
        E[i*n+j] = A[i*n+j];
    }
}
fclose(file);
long long start,end;
double duration;
start=get_nanoseconds();
assign_I(EV,n);
QR_Method(E,EV,n);
Schur_to_Eigenvectors(E,EV,n);
invert_matrix(EV,EVR,n);
assign_I(D,n);
for(int i=0;i<n;i++){
    D[i*n+i]=E[i*n+i];
}//zamiast liczyć osobno macierz D wzorem D=EV*A*EVR, zauważyłem że macierz D jest macierzą diagonalną składającą się z wartości własnych macierzy A.
write_to_file("EVR_c.txt",EVR,n);
write_to_file("EV_c.txt",EV,n);
write_to_file("D_c.txt",D,n);
//koniec programu, obliczone wszystkie składniki dekompozycji, czyli A=EV(P)*D*EVR(P^-1).
end=get_nanoseconds();
duration=(end-start)/1e9;
printf("czas wykonania algorytmu: %.9lf sekund",duration);
free(A);
free(E);
free(EV);
free(EVR);
free(D);
return 0;
}
