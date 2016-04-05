function i=itcdg(seqs,sig,l,tau,alpha,testid)
mkdir(['./' testid '/'])


%% Loading required Data
load('./data/pdtv_kalman_data.mat'); % results of traking and observations of objects
load('C:\Users\vahid.bastani\Documents\MyDrive\Datasets\PDTV\Matlab\M.mat'); % detection times
%% Adding path of requred directories
% path to GPML library
addpath(genpath('C:\Users\vahid.bastani\Documents\MyDrive\Projects\gpml-matlab-v3.5-2014-12-08'))

%% code control related parameters
monitor = 0; % 1: plot intermediate algorithm data 0: no

%% Algotithm related parameters and variables

%1) RBPF
par_num = 500; % number ofm particles
par_min = 500; % minimum number of effective particles
par_nin = 10000; % number of initioal paricles
resampler = makedist('Multinomial'); % resampler distribution

%x_p
%gamma
catdist = makedist('Multinomial'); % categorical distribution object

%2) GP
%sig = 15;
%l = 1;
likfunc = @likGauss; % Gaussian liklihood function
covfunc = @covSEard; % Squerd Exponential covariance function with ARD
prior.lik = {'priorDelta'}; % prior distribution of liklihood
prior.cov  = {{@priorGauss,l,1};... % prior distribution of length scale 1st coordinate
    {@priorGauss,l,1};... % prior distribution of length scale 2nd coordinate
    {@priorGauss,l,1};
    {'priorDelta'}}; % prior distribution of scale
inf = {@infPrior,@infExact,prior}; % inference method of GP
inf2 = {@infExact};
%3) BCM
cat(1).x = [];
cat(1).svx = [];
cat(1).vx = [];
cat(1).svy = [];
cat(1).vy = []; % initial empty DP subset data

cat(1).hypex.cov = [];
cat(1).hypex.lik = [];
cat(1).hypey.cov = [];
cat(1).hypey.lik = []; % initial empty DP subset hyperparameters

%4) DP
%alpha = 2; % concentration parameter

K = 0; % initial number of clusters
%tau = .15;
%5) Clustering
cls_mat = zeros(0,0); %initial empty clustering matrix

%% Execution variables (only used for introduction)
t = []; % time index
q = []; % sequence index
seq = []; % current sequence being processed
i = []; % general lop variable
j = []; % general lop variable
dt = []; % detection time of obervations
cnt1 = []; % general counter

%% Global loop


for q=1:length(seqs) % for all sequences
    tic
    seq = seqs(q)
    %CCC
    F = zeros(1,K+1);
    %CCC
    % find detection times
    dt = M(seq,find(M(seq,:) == 1):find(M(seq,:) == 3))>0;
    
    % TRACKER -------------------------
    % tracker initialization
    x_i = rand(2,par_nin).*repmat([640 480]',[1 par_nin]); % random initial particles
    for j =1:par_nin
        d(j) = (x_i(:,j)-Data(seq).obs_mean(:,1))'...
            /Data(seq).obs_cov{1}*...
            (x_i(:,j)-Data(seq).obs_mean(:,1));
    end % for j
    [val, idx] = sort(d);
    for i = 1:par_num
        xp{i}(:,1) = x_i(:,idx(i));
        w(i) = exp(-val(i));
        for k = 1:K
            %gamma{i}(k,a) = log(1e-5+sum(cls_mat(:,k,a)).*inpdf(cat, xp{i}(:,1),cls_mat(:,k,a)'));
            %gamma{i}(k,a) = log(sum(cls_mat(:,k,a)));
            gamma{i}(k,1) = 0;
        end
        %gamma{i}(K+1,a) = log(alpha(a));
        gamma{i}(K+1,1) = 0;
    end % for i
    w = w/sum(w);
    clear idx val d x_i
    
    cnt1 = 1;
    tflo = 0; % time from last observation;
    for t=2:1:length(dt) % for all times >= 2
        
        %tic
        if dt(t)
            cnt1 = cnt1+1;
        end
        
        % resampling
        k=randcat(w,par_num);
        x_old = xp;
        w_old = w;
        gamma_old = gamma;
        for i = 1:par_num
            if t==1
                xp{i} = x_old{k(i)};
                w(i) = w_old(k(i));
                gamma{i} = gamma_old{k(i)};
            else
                xp{i}(:,t-1) = x_old{k(i)}(:,t-1);
                gamma{i} = gamma_old{k(i)};
                w(i) = w_old(k(i));
            end
        end
        w = w/sum(w);
        
        for i = 1:par_num
            XI(i,:) = [xp{i}(:,t-1)' (t-1)/21];
        end
        
        idx = 1;
        %tic
        if q>1
            [fx(:,idx:idx+q-2), fy(:,idx:idx+q-2), sig2x(:,idx:idx+q-2), sig2y(:,idx:idx+q-2)] =...
                bcm23(cat, inf2, covfunc, likfunc,XI);
            
        end
        %et(t) = toc;
        idx = q;
        fx(1:par_num,idx) = 0; fy(1:par_num,idx) = 0; sig2x(1:par_num,idx) = sig; sig2y(1:par_num,idx) = sig;
        
        
        % propgate particles
        for i =1:par_num
            
            % sampling
            if q>1
                temp = gamma_correction(gamma{i},t,tau,alpha,cls_mat(:,:)',F);
                gamma_n = loglik2prob(temp);
                rho = sampling_mix_weights(gamma_n,cls_mat(:,:)');
                r = randcat(rho,1);
            else
                r =1;
                gamma_n = 1;
                rho =1;
            end
            vp{i}(:,t-1) = mvnrnd2([fx(i,r) fy(i,r)], diag([sig2x(i,r),sig2y(i,r)]))';
            
            xp{i}(:,t) = xp{i}(:,t-1)+vp{i}(:,t-1);
            
            if dt(t) % if observation is avaliable at this time
                nu = 1e-6+mvnpdf2(xp{i}(:,t)',Data(seq).obs_mean(:,cnt1)',Data(seq).obs_cov{cnt1});
            else % if dt(t)
                nu = 1;
            end % if dt(t)
            
            idx = 1:q;
            g_s = 1e-4+mvnpdf4(vp{i}(:,t-1)',[fx(i,idx)' fy(i,idx)'], sig2x(i,idx)',sig2y(i,idx)');
            % \propto p(x(i),k|x(i)_, y_)
            
           
            lambda = bcmcalc_2(g_s,cls_mat(:,:)');
            gamma{i} = gamma{i}+log(lambda) + log(nu); % nu eliminated
            
            wi(i) = nu*(gamma_n'*lambda)/(rho*g_s) * w(i);
            
            %CCC
            fkk(:,i) = lambda-lambda(K+1);
            %CCC
        end %for i
        w = wi./sum(wi);
        %CCC
        F = F+(w*fkk' < 0.001);
        %CCC
        if dt(t)
            % sttimtion speed
            
            for i=1:par_num
                x1x(i) = xp{i}(1,t-tflo-1);
                x2x(i,1) = xp{i}(1,t);
                x1y(i) = xp{i}(2,t-tflo-1);
                x2y(i,1) = xp{i}(2,t);
            end
            v1t = (repmat(x2x,[1 par_num]) - repmat(x1x,[par_num 1]))./(tflo+1);
            v2t = (repmat(x2y,[1 par_num]) - repmat(x1y,[par_num 1]))./(tflo+1);
            w_vt = w'*w;
            
            vt = [sum(sum(v1t.*w_vt));sum(sum(v2t.*w_vt))];
            sv = sqrt([sum(sum(v1t.^2.*w_vt))-vt(1); sum(sum(v2t.^2.*w_vt))-vt(2)]);
            
            V(:,t-tflo-1:t-1) = repmat(vt,[1 tflo+1]);
            %     plot( xp{i}(1,1:t), xp{i}(2,1:t))
            %     hold on
            SV(:,t-tflo-1:t-1) = repmat(sv,[1 tflo+1]);
            tflo = 0;
        else
            tflo = tflo+1;
        end


        G  =  gamma{i} * 0;
        for i=1:par_num
            G = G+ gamma{i}.*w(i);
        end
        G = gamma_correction(G,t,tau,alpha,cls_mat(:,:)',F);
        Gtl(:,t) = G;
        Gt(:,t) = loglik2prob(G);
    end % for t
    
    X  =  xp{i} * 0;
    for i=1:par_num
        X = X+ xp{i}.*w(i);
    end
    X = X(:,1:end-1);
    X(3,:) =  ((1:length(X))-1)/21;
    % DP
    cat(q).x =  X';
    cat(q).vx = V(1,:)';
    cat(q).vy = V(2,:)';
    cat(q).svx = SV(1,:)';
    cat(q).svy = SV(2,:)';
    
    cat(q).hypex.cov = [1; 1; 1; log(sqrt(sig))];
    cat(q).hypey.cov = [1; 1; 1; log(sqrt(sig))];
    
    if mean(cat(q).svx) > 0
        cat(q).hypex.lik = log(mean(cat(q).svx));
    else
        cat(q).hypex.lik = log(1.5);
    end
    
    
    if mean(cat(q).svy) > 0
        cat(q).hypey.lik = log(mean(cat(q).svy));
    else
        cat(q).hypey.lik = log(1.5);
    end
    
    cat(q).hypex = minimize(cat(q).hypex, @gp, -200, inf, [],...
        covfunc, likfunc, cat(q).x, cat(q).vx);
    cat(q).hypey = minimize(cat(q).hypey, @gp, -200, inf, [],...
        covfunc, likfunc, cat(q).x, cat(q).vy);
    
    
    N0 = 3;
    cat(q).m0 = mean(X(1:2,1:N0),2);
    cat(q).S0 = cov((X(1:2,1:N0)-repmat(cat(q).m0,[1,N0]))')+diag([10 10]);
    
    if q == 1
        cls_mat(1,1,:) = 1;
        K = K+1;
    else
        
        G  =  gamma{i} * 0;
        for i=1:par_num
            G = G+ gamma{i}.*w(i);
        end
        G = gamma_correction(G,t,tau,alpha,cls_mat(:,:)',F);
        G = loglik2prob(G);
        if sum(G(K+1,:)) > .1
            for k=1:K+1
                cls_mat(q,k) = G(k)/sum(G);
            end
            K= K+1;
        else
            G= G(1:K,:);
            G = G./sum(G(:));
            for k=1:K
                cls_mat(q,k) = G(k)/sum(G(:));
            end
        end
        
    end
    ET(q) = toc
    F
    cls_mat
    save(['./' testid '/' num2str(seq) '_track.mat'],'xp','w','Gtl','Gt')
    clear xp v_1 v_2 w t x_old w_old wi w_v v1t v2t w_vt x1x x2x x1y x2y  X V vt sv et g_s SV fkk F Gt B E Gtl
    
end %for q

for q = 1:length(seqs)
    seq = seqs(q);
    [i,C(q)] = max(cls_mat(q,:));
    GC(q) = Data(seq).gtc;
end

mi = ami(GC,C)
I = zeros(max(C),max(GC));
for i=1:length(C)
    I(C(i),GC(i)) = 1+I(C(i),GC(i));
end

acc = mean(max(I')./sum(I'))

save(['./' testid '/final.mat'])