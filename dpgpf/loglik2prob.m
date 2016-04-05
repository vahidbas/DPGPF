function p = loglik2prob(ll)

for i=1:length(ll)
    if isinf(ll(i))
        p(i,1) = 0;
    else
         p(i,1) = 1/(sum(exp(ll-ll(i))));
    end
end