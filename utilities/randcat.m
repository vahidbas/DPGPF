function r = randcat(p,n)
if abs(sum(p)- 1) > 1e-10
    error('p should some to one');
end
u = rand(n,1);
c = zeros(n,1);
r = zeros(n,1);
for i=1:length(p)
    c = c+p(i);
    
    r = (u <= c).*(r.*(r ~=0)+i.*(r == 0));
%     for k=1:n
%         if u(k)<=c(k) && r(k) == 0
%         r(k) = i;
%         end
%     end
end

