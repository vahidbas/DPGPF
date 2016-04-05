function p = mvnpdf2(x,m,S)

k = length(x);

p = exp(-.5*(x-m)/S*(x-m)')/sqrt((2*pi)^k*det(S));