function t = eigen(A)
tic
d=1e-9;
n=length(A);
V=eye(n);
%metoda QR:
for i=1:50
    sigma=A(n,n);
    A=A-sigma*eye(n);
    Q=A;
    R=zeros(n);
    for j=1:n
        R(j,j)=norm(Q(:,j));
        Q(:,j)=Q(:,j)/R(j,j);
        R(j,j+1:n)=Q(:,j)'*Q(:,j+1:n);
        Q(:,j+1:n)=Q(:,j+1:n)-Q(:,j)*R(j,j+1:n);
    end
    A=R*Q+sigma*eye(n);
    V=V*Q;
    if (max(max(abs(tril(A,-1))))<d)
        fprintf("osiagnieto zbieznosc po %d iteracjach.",i)
        break;
    end
end
%Schur na wektory własne
tolerance=1e-12;
Y=eye(n);
for j=1:n
    lambda=A(j,j);
    for i=(j-1):-1:1
        S=A(i,i+1:j)*Y(i+1:j,j);
        diff=lambda-A(i,i);
        if abs(diff)<tolerance
            if diff>=0
                diff=tolerance;
            else
                diff=-tolerance;
            end
        end
        Y(i,j)=S/diff;
    end
    norm_Y=norm(Y(1:j,j));
    if norm_Y>tolerance
        Y(1:j,j)=Y(1:j,j)/norm_Y;
    end
end
V=V*Y;
VR=V';
writematrix(VR,"EVR_m.txt")
writematrix(diag(A),"D_m.txt")
writematrix(V,"EV_m.txt")
t=toc;
