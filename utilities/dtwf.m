function d = dtwf(a,b,p)

n = size(a,1);
w = p{1};
D = p{2};

for i=1:n
    Na = find(isnan(a(i,:)),1)-1;
    Nb = find(isnan(b(i,:)),1)-1;
    A = reshape(a(i,1:Na),[Na/D,D]);
    B = reshape(b(i,1:Nb),[Nb/D,D]);
    d(i,1) = dtw(A,B,w);
end

