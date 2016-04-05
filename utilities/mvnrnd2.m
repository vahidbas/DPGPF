function  f = mvnrnd2(m, S)
if ~(S(1,1) >0 && S(2,2) >0)
    error('bad covariance matrix')
end
f(1) = randn*sqrt(S(1,1))+m(1);
f(2) = randn*sqrt(S(2,2))+m(2);