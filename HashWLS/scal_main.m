clear all
clc
sizes=[5000,10000,15000,20000];
iterations = 5;
hashDims = 50:50:300;

method = 'hashwls';
times = zeros(length(sizes),iterations);
accs = zeros(length(sizes),iterations);
for i = 1:length(sizes)
    
   dataname =  'STH_DHFR';
   load(['data/', dataname, '/',dataname, '.mat']);
   
   for iteration = 1:iterations
        rng(1)
       [accs(i,iteration),times(i,iteration)] = scal(graphs, labels,iteration,sizes(i),300);

   end  
   
   if ~exist(['results/', dataname, '/'], 'dir')
        mkdir(['results/', dataname, '/']);
   end
   %accs_mean = mean(accs, 3);
   %cpus_mean = mean(cpus, 3);
   save(['results/','scal.mat'],'times','accs')
end