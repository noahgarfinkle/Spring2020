function a= poly_multiindex (p,d)
% 0-th order polynomial
a(1 ,1:p)= zeros (1,p);
% First order polynomials
a(2:p+1 ,1:p)= eye(p,p);
P=p+1;
pmat =[];
pmat (1:p ,1)=1;
disp(pmat)
for k=2:d
    L=P;
    for i=1:p
        pmat (i,k)= sum( pmat (i:p,k -1));
    end
    for j=1:p
        for m=L- pmat (j,k )+1: L
            disp(m)
            P=P+1;
            a(P ,1:p)=a(m ,1:p);
            a(P,j)=a(P,j)+1;
        end
    end
end
disp(pmat)