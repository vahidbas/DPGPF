function gamma = gamma_correction(gamma,t,tau,alpha,mu,F)
F = F./(t-1);
%a=1/2;
%b=.1;
%C = 1./(1+exp((F-a)/b))+1./(1+exp(a/b));
%C = (F>0.15)*0.001+(F<0.15);
C=.5*erf(sqrt(pi/2)*-(F-tau)*5)+.5;
if isempty(mu)
    N=[];
else
    N = sum(mu,2);
end

%gamma = gamma +[N; t*tau+alpha];
%gamma = gamma/t;
gamma = gamma + log([C(1:end-1)';1])*t+log([N; alpha]);