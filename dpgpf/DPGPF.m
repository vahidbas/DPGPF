classdef DPGPF < handle
    properties (Access=private)
        % dataset parameters
        image_width; % video frame width (pixels)
        image_height; % video frame hight (pixels)
        Data; % the dataset
        seqs; % list of sequence to be proccessed
        
        % Particle filter parameters
        num_paticle; % number of particles
        min_paticle; % minimum number of effective particles
        nin_particle; % number of initioal paricles
        det_time;
        
        % Dirichlet process parameters
        correction_param; % tau
        dp_concentration; % alpha
        cluster_ind; % cluster indicator
        num_cluster % number of clusters
        
        % Gaussian process parameters
        gp_lenghscale; % l
        gp_likfunc; % Gaussian process liklihood function
        gp_covfunc; % Gaussian process covariance function
        gp_prior; % prior distribution of gaussian process
        gp_inf_learning; % inference method for learning
        gp_inf_prediction; % inference method for prediction
        gp_var; % prior variance
        
        % other parameters
        correction_m % the measure for correction term
    end
    
    properties (Access=private)
        memory;
        in_parser;
    end
    
    methods
        % Constructor
        function obj = DPGPF(im_size, num_par, eff_par_percent,...
                            prior_var, length_scale, tau, alpha)
            
            obj.num_paticle = num_par;
            obj.min_paticle = eff_par_percent/100*num_par;
            obj.nin_particle = obj.num_paticle*20;
            
            obj.image_width = im_size(2);
            obj.image_height = im_size(1);
            
            obj.num_cluster = 0; % no initial cluster
            obj.cluster_ind = zeros(0,0);
            
            obj.gp_var = prior_var;            
            obj.correction_param = tau;
            obj.dp_concentration = alpha;
            obj.gp_lenghscale = length_scale;
            
            obj.gp_likfunc = @likGauss; % Gaussian liklihood function
            obj.gp_covfunc = @covSEard; % Squerd Exponential covariance function with ARD
            obj.gp_prior.lik = {'priorDelta'}; % prior distribution of liklihood
            obj.gp_prior.cov  = {{@priorGauss,obj.gp_lenghscale,1};... % prior distribution of length scale 1st coordinate
                {@priorGauss,obj.gp_lenghscale,1};... % prior distribution of length scale 2nd coordinate
                {@priorGauss,obj.gp_lenghscale,1};
                {'priorDelta'}}; % prior distribution of scale
            obj.gp_inf_learning = {@infPrior,@infExact,obj.gp_prior}; % inference method of GP
            obj.gp_inf_prediction = {@infExact};
            
            % configure input parser
            obj.in_parser = inputParser;
            
            default_plot = 'off';
            valid_plot = {'on','off'};
            checkPlot = @(x) any(validatestring(x,valid_plot));
            
            default_WriteOutput = 'off';
            valid_WriteOutput = {'on','off'};
            check_WriteOutput = @(x) any(validatestring(x,valid_WriteOutput));
            
            default_OutputFile = 'result.mat';
            default_TrackOutputDirectory = './results/';
            
            addRequired(obj.in_parser,'Sequences',@isnumeric)
            addOptional(obj.in_parser,'Plot',default_plot,checkPlot);
            addOptional(obj.in_parser,'WriteOutput',default_WriteOutput,check_WriteOutput);
            addOptional(obj.in_parser,'OutputFile',default_OutputFile,@ischar); 
            
            
        end
        
        function processDataSet(obj, Data, varargin)
            obj.Data = Data;
            
            % parse input
            parse(obj.in_parser,varargin{:})
            
            %make folder fo tracking results
            if isequal(obj.in_parser.Results.WriteOutput, 'on')
                [fpathstr,fname,~] = fileparts(obj.in_parser.Results.OutputFile);
                fnfull = [fpathstr '/' fname];
                if exist(fnfull, 'dir')
                    warning([obj.in_parser.Results.OutputFile '/ will be removed'])
                    rmdir(fnfull,'s');                
                end
                mkdir(fnfull)
            end
                
            
            
            obj.seqs = obj.in_parser.Results.Sequences;
            active_plot = 0;
            if isequal(obj.in_parser.Results.Plot, 'on')
                active_plot = 1;
            end
            
            for q=1:length(obj.seqs)
                processSequnce(obj,q);
                
                if active_plot == 1
                    col = lines(50);
                    for i=1:q
                        seq_i = obj.seqs(i);
                        [~,C(i)] = max(obj.cluster_ind(i,:));
                        plot(Data(seq_i).obs_mean(1,:),...
                             Data(seq_i).obs_mean(2,:),...
                             '.',...
                             'Color', col(C(i),:));
                        hold on
                    end
                    hold off
                    xlim([0 obj.image_width])
                    ylim([0 obj.image_width])
                    pause(0.1)
                end
            end
            
            for q = 1:length(obj.seqs)
                seq = obj.seqs(q);
                [~,C(q)] = max(obj.cluster_ind(q,:));
                GC(q) = obj.Data(seq).groundtruth_class;
            end
            
            mi = ami(GC,C)
            I = zeros(max(C),max(GC));
            for i=1:length(C)
                I(C(i),GC(i)) = 1+I(C(i),GC(i));
            end

            acc = mean(max(I')./sum(I'))
            K = max(C);
            if isequal(obj.in_parser.Results.WriteOutput, 'on')
                save(['./' obj.in_parser.Results.OutputFile],...
                     'C','GC','I','acc','mi','K');
            end
        end
    end
    
    methods (Access=private)
        
        function tracker = initializeParticles(obj, init_obs_mean,...
                init_obs_cov,...
                num_clasters)
            
            Ni = obj.nin_particle;
            N = obj.num_paticle;
            imW = obj.image_width;
            imH = obj.image_height;
            K = num_clasters;
            
            % random initial particles
            x_i = rand(2,Ni).*repmat([imW imH]',[1 Ni]);
            
            for j =1:Ni
                d(j) = (x_i(:,j)-init_obs_mean)'...
                    /init_obs_cov*...
                    (x_i(:,j)-init_obs_mean);
            end % for j
            [val, idx] = sort(d);
            for i = 1:N
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
            tracker.w = w/sum(w);
            tracker.xp = xp;
            tracker.gamma = gamma;
        end % initializeParticles
        
        function [tracker,tflo,cnt1] = filterStep(obj,tracker,q,t,tflo,cnt1)
            seq = obj.seqs(q);
            Np = obj.num_paticle;
            
            
            if obj.det_time(t)
                cnt1 = cnt1+1;
            end
            tic % statit ticking
            % resampling
            k=randcat(tracker.w,Np);
            x_old = tracker.xp;
            w_old = tracker.w;
            gamma_old = tracker.gamma;
            for i = 1:Np
                if t==1
                    tracker.xp{i} = x_old{k(i)};
                    w(i) = w_old(k(i));
                    tracker.gamma{i} = gamma_old{k(i)};
                else
                    tracker.xp{i}(:,t-1) = x_old{k(i)}(:,t-1);
                    tracker.gamma{i} = gamma_old{k(i)};
                    w(i) = w_old(k(i));
                end
            end
            tracker.w = w/sum(w);
            
            % put time in the state vector
            for i = 1:Np
                XI(i,:) = [tracker.xp{i}(:,t-1)' (t-1)/21];
            end
            
            idx = 1;
            %tic
            if q>1
                [fx(:,idx:idx+q-2), fy(:,idx:idx+q-2), sig2x(:,idx:idx+q-2), sig2y(:,idx:idx+q-2)] =...
                    bcm23(obj.memory, obj.gp_inf_prediction, obj.gp_covfunc, obj.gp_likfunc,XI);
                
            end
            
            idx = q;
            fx(1:Np,idx) = 0; fy(1:Np,idx) = 0;
            sig2x(1:Np,idx) = obj.gp_var; sig2y(1:Np,idx) = obj.gp_var;
            
            
            % propgate particles
            for i =1:Np
                
                % sampling
                if q>1
                    temp = gamma_correction(tracker.gamma{i},t,...
                                obj.correction_param,...
                                obj.dp_concentration,...
                                obj.cluster_ind(:,:)',...
                                obj.correction_m);
                    gamma_n = loglik2prob(temp);
                    rho = sampling_mix_weights(gamma_n,obj.cluster_ind(:,:)');
                    r = randcat(rho,1);
                else
                    r =1;
                    gamma_n = 1;
                    rho =1;
                end
                vp{i}(:,t-1) = mvnrnd2([fx(i,r) fy(i,r)], diag([sig2x(i,r),sig2y(i,r)]))';
                
                tracker.xp{i}(:,t) = tracker.xp{i}(:,t-1)+vp{i}(:,t-1);
                
                if obj.det_time(t) % if observation is avaliable at this time
                    nu = 1e-6+mvnpdf2(tracker.xp{i}(:,t)', ...
                        obj.Data(seq).obs_mean(:,cnt1)',...
                        obj.Data(seq).obs_cov{cnt1});
                else % if dt(t)
                    nu = 1;
                end % if dt(t)
                
                idx = 1:q;
                g_s = 1e-4+mvnpdf4(vp{i}(:,t-1)',[fx(i,idx)' fy(i,idx)'], sig2x(i,idx)',sig2y(i,idx)');
                % \propto p(x(i),k|x(i)_, y_)
                
                
                lambda = bcmcalc_2(g_s,obj.cluster_ind(:,:)');
                tracker.gamma{i} = tracker.gamma{i}+log(lambda) + log(nu); % nu eliminated
                
                wi(i) = nu*(gamma_n'*lambda)/(rho*g_s) * w(i);
                
                %CCC
                fkk(:,i) = lambda-lambda(obj.num_cluster+1);
                %CCC
            end %for i
            tracker.w = wi./sum(wi);
            %CCC
            obj.correction_m = obj.correction_m+(w*fkk' < 0.001);
            %CCC
            
            if obj.det_time(t)
                % sttimtion speed
                
                for i=1:Np
                    x1x(i) = tracker.xp{i}(1,t-tflo-1);
                    x2x(i,1) = tracker.xp{i}(1,t);
                    x1y(i) = tracker.xp{i}(2,t-tflo-1);
                    x2y(i,1) = tracker.xp{i}(2,t);
                end
                v1t = (repmat(x2x,[1 Np]) - repmat(x1x,[Np 1]))./(tflo+1);
                v2t = (repmat(x2y,[1 Np]) - repmat(x1y,[Np 1]))./(tflo+1);
                w_vt = tracker.w'*tracker.w;
                
                vt = [sum(sum(v1t.*w_vt));sum(sum(v2t.*w_vt))];
                sv = sqrt([sum(sum(v1t.^2.*w_vt))-vt(1); sum(sum(v2t.^2.*w_vt))-vt(2)]);
                
                tracker.V(:,t-tflo-1:t-1) = repmat(vt,[1 tflo+1]);
                 
                tracker.SV(:,t-tflo-1:t-1) = repmat(sv,[1 tflo+1]);
                tflo = 0;
            else
                tflo = tflo+1;
            end
%             figure(1)
%             for i =1:Np
%             plot( tracker.xp{i}(1,1:t), tracker.xp{i}(2,1:t))
%             hold on
%             end
            
            G  =  tracker.gamma{i} * 0;
            for i=1:Np
                G = G+ tracker.gamma{i}.*w(i);
            end
            G = gamma_correction(G,t,...
                obj.correction_param,...
                obj.dp_concentration,...
                obj.cluster_ind(:,:)',...
                obj.correction_m);
            tracker.Gtl(:,t) = G;
            tracker.Gt(:,t) = loglik2prob(G);
            
            tracker.frame_processing_time(t) = toc;% stop ticking
        end % filterStep
        
        function processSequnce(obj,q)
            seq = obj.seqs(q)
            M = obj.Data(seq).detection_time;
            obj.det_time = M(find(M(:) == 1):find(M(:) == 3))>0;
            obj.correction_m = zeros(1,obj.num_cluster+1);
            
            tracker = initializeParticles(obj, ...
                obj.Data(seq).obs_mean(:,1), obj.Data(seq).obs_cov{1},obj.num_cluster);
            
            cnt1 = 1;
            tflo = 0; % time from last observation;
            for t=2:1:length(obj.det_time) % for all times >= 2               
                [tracker, tflo,cnt1] = filterStep(obj,tracker,q,t,tflo,cnt1);
%                 figure(1)
%                 hold off
                
            end
            
            
            X  =  tracker.xp{1} * 0;
            for i=1:obj.num_paticle
                X = X+ tracker.xp{i}.*tracker.w(i);
            end
            X = X(:,1:end-1);
            X(3,:) =  ((1:length(X))-1)/21;
            % DP
            obj.memory(q).x =  X';
            obj.memory(q).vx = tracker.V(1,:)';
            obj.memory(q).vy = tracker.V(2,:)';
            obj.memory(q).svx = tracker.SV(1,:)';
            obj.memory(q).svy = tracker.SV(2,:)';
            
            obj.memory(q).hypex.cov = [1; 1; 1; log(sqrt(obj.gp_var))];
            obj.memory(q).hypey.cov = [1; 1; 1; log(sqrt(obj.gp_var))];
            
            if mean(obj.memory(q).svx) > 0
                obj.memory(q).hypex.lik = log(mean(obj.memory(q).svx));
            else
                obj.memory(q).hypex.lik = log(1.5);
            end
            
            
            if mean(obj.memory(q).svy) > 0
                obj.memory(q).hypey.lik = log(mean(obj.memory(q).svy));
            else
                obj.memory(q).hypey.lik = log(1.5);
            end
            
            obj.memory(q).hypex = minimize(obj.memory(q).hypex,...
                @gp, -200, obj.gp_inf_learning, [],...
                obj.gp_covfunc, obj.gp_likfunc, ...
                obj.memory(q).x, obj.memory(q).vx);
            obj.memory(q).hypey = minimize(obj.memory(q).hypey,...
                @gp, -200, obj.gp_inf_learning, [],...
                obj.gp_covfunc, obj.gp_likfunc,...
                obj.memory(q).x, obj.memory(q).vy);
            
            
            N0 = 3;
            obj.memory(q).m0 = mean(X(1:2,1:N0),2);
            obj.memory(q).S0 = cov((X(1:2,1:N0)-repmat(obj.memory(q).m0,[1,N0]))')+diag([10 10]);
            
            if q == 1
                obj.cluster_ind(1,1,:) = 1;
                obj.num_cluster = obj.num_cluster+1;
            else
                
                G  =  tracker.gamma{i} * 0;
                for i=1:obj.num_paticle
                    G = G+ tracker.gamma{i}.*tracker.w(i);
                end
                G = gamma_correction(G,t,...
                obj.correction_param,...
                obj.dp_concentration,...
                obj.cluster_ind(:,:)',...
                obj.correction_m);
                G = loglik2prob(G);
                if sum(G(obj.num_cluster+1,:)) > .1
                    for k=1:obj.num_cluster+1
                        obj.cluster_ind(q,k) = G(k)/sum(G);
                    end
                    obj.num_cluster= obj.num_cluster+1;
                else
                    G= G(1:obj.num_cluster,:);
                    G = G./sum(G(:));
                    for k=1:obj.num_cluster
                        obj.cluster_ind(q,k) = G(k)/sum(G(:));
                    end
                end
                
            end
            
            if isequal(obj.in_parser.Results.WriteOutput, 'on')
                [fpathstr,fname,~] = fileparts(obj.in_parser.Results.OutputFile);
                fnfull = [fpathstr '/' fname];
                file_name = [fnfull '/' num2str(seq) '_track.mat'];
                xp = tracker.xp;
                w = tracker.w;
                Gt = tracker.Gt;
                Gtl = tracker.Gtl;
                fpt = tracker.frame_processing_time;
                save(file_name,'xp','w','Gtl','Gt','fpt')
            end
                        
        end % processSequnce
    end
end