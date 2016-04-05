function [fx, fy, sig2x, sig2y] = bcm23(data, inf, covfunc, likfunc, x)

P = length(data);

DS = 10;
for i=1:length(data)
     L = length(data(i).vx);
     DS = floor(L/20);
     [~, ~, fx(:,i), sig2x(:,i) ] = gp(data(i).hypex, inf, [], covfunc,...
             likfunc, data(i).x(1:DS:end,:), data(i).vx(1:DS:end,:), x);
     [~, ~, fy(:,i), sig2y(:,i) ] = gp(data(i).hypey, inf, [], covfunc,...
             likfunc, data(i).x(1:DS:end,:), data(i).vy(1:DS:end,:), x);

end

