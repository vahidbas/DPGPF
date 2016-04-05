function [varargout] = wgp(hyp, w, inf, mean, cov, lik, x, y, xs, ys)
% x : cell array of inputs categories
% y : cell array of output categories

% input check copied from gp
if nargin<8 || nargin>10
  disp('Usage: [nlZ dnlZ          ] = wgp(hyp, w, inf, mean, cov, lik, x, y);')
  disp('   or: [ymu ys2 fmu fs2   ] = wgp(hyp, w, inf, mean, cov, lik, x, y, xs);')
  disp('   or: [ymu ys2 fmu fs2 lp] = wgp(hyp, w,, inf, mean, cov, lik, x, y, xs, ys);')
  return
end


K = length(w); % number of categories
w=w/sum(w);

for i = 1:K
    if nargin == 8
        [nlZi{i}, dnlZi{i}] = gp(hyp, inf, mean, cov, lik, x{i}, y{i});
    elseif nargin == 9
        [ymui{i}, ys2i{i}, fmui{i}, fs2i{i}   ] = gp(hyp, inf, mean, cov, lik, x{i}, y{i}, xs);
    elseif nargin == 10
        [ymui{i} ys2i{i} fmui{i} fs2i{i} lpi{i}] = gp(hyp, inf, mean, cov, lik, x{i}, y{i}, xs, ys);
    end
end

if nargin == 8
    fields = fieldnames(dnlZi{i});
    for f=1:length(fields)
        if ~isempty(dnlZi{1,i}.(fields{f}))
            dnlZ.(fields{f})=0;
        else
            dnlZ.(fields{f})=[];
        end
    end
    for i = 1:K
        g(i) = exp(-nlZi{i});
    end
    Z = sum(g.*w);
    for i = 1:K        
        for f=1:length(fields)
            if ~isempty(dnlZi{1,i}.(fields{f}))
                dnlZ.(fields{f}) = (-Z*dnlZ.(fields{f})-g(i).*dnlZi{1,i}.(fields{f})*w(i))/-Z;
            end
        end
    end
    nlZ = -log(Z);
    varargout =  {nlZ, dnlZ};
end

if nargin > 8
    ymu = zeros(size(ymui{1}));
    ys2 = zeros(size(ys2i{1}));
    fmu = zeros(size(fmui{1}));
    fs2 = zeros(size(fs2i{1}));    
    for i = 1:K
        ymu = ymu+ymui{i}*w(i);
        ys2 = ys2+ys2i{i}*w(i);
        fmu = fmu+fmui{i}*w(i);
        fs2 = fs2+fs2i{i}*w(i);
    end
    if nargin == 10
        lp = zeros(size(lpi{1}));
        for i = 1:K
            lp = lp+exp(lpi{i})*w(i);
        end
        lp = log(lp);
        
        varargout = {ymu, ys2, fmu, fs2, lp};
    else
        varargout = {ymu, ys2, fmu, fs2};
    end
end