function D = lcss(X,Y,w)
%%%Calculates the longest common substring between to strings.
%%%Code written by David Cumin
%%%email: d.cumin@auckland.ac.nz
%%%INPUT
%%%X, Y - both are strings e.g. 'test' or 'stingtocompare'
%%%OUTPUT
%%%D is the substring over the length of the shortest string
%%%dist is the length of the substring
%%%aLongestString is a sting of length dist (only one of potentially many)

%%%For example
%%% X = 'abcabc';
%%% Y = 'adcbac';
%%% [D dist str] = LCS(X,Y);
%%% results in:
%%% D = 0.6667 
%%% dist = 4
%%% str = acbc
%%% this is seen for X: 'a-c-bc' and Y: 'a-cb-c'

%%%Make matrix
n =size(X,1);
m =size(Y,1);
L=zeros(n+1,m+1);
L(1,:)=0;
L(:,1)=0;
b = zeros(n+1,m+1);
b(:,1)=1;%%%Up
b(1,:)=2;%%%Left

for i = 2:n+1
    for j = 2:m+1
        if (norm(X(i-1,:)-Y(j-1,:)) < 10)
            L(i,j) = L(i-1,j-1) + 1;
            b(i,j) = 3;%%%Up and left
        else
            L(i,j) = L(i-1,j-1);
        end
        if(L(i-1,j) >= L(i,j))
            L(i,j) = L(i-1,j);
            b(i,j) = 1;%Up
        end
        if(L(i,j-1) >= L(i,j))
            L(i,j) = L(i,j-1);
            b(i,j) = 2;%Left
        end
    end
end
L(:,1) = [];
L(1,:) = [];
b(:,1) = [];
b(1,:) = [];
dist = L(n,m);

D = 1-(dist / min(m,n));

D = exp(-D^2/2/w^2);