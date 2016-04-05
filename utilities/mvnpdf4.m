function p = mvnpdf4(x,m,S1,S2)
X = repmat(x,size(m,1),1);
p = prod(exp(-(X-m).^2./[S1 S2]./2)./sqrt(2*pi*[S1 S2]),2);
