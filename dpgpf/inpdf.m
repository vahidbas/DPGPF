function p = inpdf(data, x,mu)
m = zeros(size(data(1).m0));
S = zeros(size(data(1).S0));
for i=1:length(data)
    m = m+data(i).m0*mu(i);
    S = S+data(i).S0*mu(i);
end
m = m./sum(mu);
S = S./sum(mu);

p = mvnpdf(x,m,S);