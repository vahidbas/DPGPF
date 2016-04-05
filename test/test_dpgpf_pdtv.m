
sig = [15 30];
l = [1 2];
tau = [0.5 .15 .25];
alpha = [.2 .5 2 5];

data_set = 'pdtv';
load([data_set '.mat']);
seqs = 1:length(Data);
out_put_dir = ['output/result_' data_set];

num_particle = 500;
eff_particle_percent = 80;

for i1=1:length(sig)
  for i2=1:length(l)
      for i3=1:length(tau)
          for i4=1:length(alpha)
              output_file = [out_put_dir '_' num2str(sig(i1)) '_' ...
                             num2str(l(i2)) '_' num2str(tau(i3)) '_' ...
                             num2str(alpha(i4)) '.mat'];
              alg = DPGPF(image_size,num_particle,...
                          eff_particle_percent,...
                          sig(i1),l(i2),tau(i3),alpha(i4));
              alg.processDataSet(Data,seqs,...
                                'WriteOutput','on',...
                                'Plot','on', ...
                                 'OutputFile',output_file);
          end
      end
  end
end
