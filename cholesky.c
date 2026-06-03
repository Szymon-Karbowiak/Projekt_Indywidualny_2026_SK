#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <windows.h>

long long get_nanoseconds(){
    LARGE_INTEGER freq,counter;
    QueryPerformanceFrequency(&freq);
    QueryPerformanceCounter(&counter);
    return (long long)(counter.QuadPart*1000000000LL/freq.QuadPart);
}

double* mem_allocate(int n){
    double *A=(double *) malloc(n*n*sizeof(double));
    if(A==NULL){
        exit(1);
    }
    return A;
}

void write_to_file(const char *filename, double *A, int n) {
    FILE *file=fopen(filename,"w");
    if(file==NULL){
        printf("Blad: Nie udalo sie utworzyc pliku %s\n",filename);
        return;
    }
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){
            fprintf(file,"%3.15f",A[i*n+j]);
            if(j<n-1){
                fprintf(file,",");
            }
        }
        fprintf(file,"\n");
    }
    fclose(file);
    printf("Zapisano macierz do pliku: %s\n",filename);
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

void transpose_matrix(double *src, double *dst, int n) {
    for (int i=0;i<n;i++){
        for (int j=0;j<n;j++){
            dst[j*n+i]=src[i*n+j];
        }
    }
}

int cholesky_decomposition(double *A, double *L, int n) {
    for (int i=0;i<n*n;i++){
        L[i]=0.0;
    }

    for (int i=0;i<n;i++){
        for (int j=0;j<=i;j++){
            double sum=0.0;
            for (int k=0;k<j;k++){
                sum+=L[i*n+k]*L[j*n+k];
            }
            if(i==j){
                double diff=A[i*n+i]-sum;
                if (diff<=0.0) {
                    return 0;
                }
                L[i*n+j]=sqrt(diff);
            }else{
                L[i*n+j]=(A[i*n+j]-sum)/L[j*n+j];
            }
        }
    }
    return 1;
}

int main(int argc, char **argv){
    double *A,*L,*LT;
    int n;
    FILE*file=NULL;
    if(argc>=3) {
        n=atoi(argv[1]);
        file=fopen(argv[2],"r");
        if (!file){
            printf("Nie udało się otworzyć pliku: %s\n",argv[2]);
            return 1;
        }
    } else {
        printf("Podaj rozmiar macierzy n: ");
        if(scanf("%d",&n)!=1){
            return 1;
        }
        char filename[128];
        printf("Podaj nazwe pliku z danymi: ");
        if(scanf("%s",filename)!=1){
            return 1;
        }
        file=fopen(filename,"r");
        while (!file){
            printf("Blad otwarcia pliku. Podaj poprawna nazwe: ");
            if (scanf("%s",filename)!=1){
                return 1;
            }
            file=fopen(filename,"r");
        }
    }
    A = mem_allocate(n);
    L = mem_allocate(n);
    LT = mem_allocate(n);
    for(int i=0;i<n;i++){
        for (int j=0;j<n;j++){
            if (fscanf(file,"%lf",&A[i*n+j])!=1){
                printf("Blad: Za malo danych w pliku dla macierzy %dx%d!\n",n,n);
                fclose(file);
                return 1;
            }
            if(j<n-1) {
                char writeme;
                if (fscanf(file,"%c",&writeme)==1 && writeme!=','){
                    printf("Ostrzeżenie: Oczekiwano przecinka, znaleziono: %c\n",writeme);
                }
            }
        }
    }
    fclose(file);

    long long start,end;
    double duration;
    start = get_nanoseconds();
    if(cholesky_decomposition(A,L,n)) {
        transpose_matrix(L,LT,n);
        write_to_file("L_c.txt",L,n);
        write_to_file("LT_c.txt",LT,n);
    }else{
        printf("Blad: Macierz nie jest symetryczna dodatnio okreslona! Dekompozycja niemozliwa.\n");
    }
    end=get_nanoseconds();
    duration=(end-start)/1e9;
    printf("czas wykonania algorytmu: %.9lf sekund\n",duration);
    free(A);
    free(L);
    free(LT);
    return 0;
}
