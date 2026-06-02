function t = cholesky(A)
tic;
n=length(A);
L = zeros(n, n);
for i = 1:n
    for j = 1:i
        sum_val = L(i, 1:j-1) * L(j, 1:j-1)';
        if i == j
            diff = A(i, i) - sum_val;
            if diff <= 0
                L = [];
                fprintf('Blad: Macierz nie jest symetryczna dodatnio okreslona! Dekompozycja niemozliwa.\n');
                return
            end
            L(i, j) = sqrt(diff);
        else
            L(i, j) = (A(i, j) - sum_val) / L(j, j);
        end
    end
end
    LT = L';
    writematrix(L, 'L_mat_opt.txt');
    writematrix(LT, 'LT_mat_opt.txt');
t=toc;
end