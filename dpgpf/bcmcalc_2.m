function lambda = bcmcalc_2(gs,mu)

Q = size(mu,2);


rho_k = mu.^2./ repmat(sum(mu,2),[1,Q]);
rho_k(:,Q+1) = 1-sum(rho_k,2);
rho_k(end+1,Q+1) = 1;

lambda = rho_k*gs;
