clear all
clc

datanames = {'BZR', 'COX2', 'DHFR', 'PROTEINS_full', 'AIDS'};

iterations = 5;
hashDims = 50:50:300;
method = 'hashwls';

accs = zeros(iterations, length(hashDims));
cpus = zeros(iterations, length(hashDims));

for idataname = 1:length(datanames)
    
   dataname =  datanames{idataname};
   load(['data/', dataname, '/',dataname, '.mat']);
   
   for r = 1:iterations
       for iHash = 1:length(hashDims)
                rng(2);
               [accs(r, iHash), cpus(r, iHash)] = hashwls(graphs, labels, r, hashDims(iHash));

           end
   end
   
   if ~exist(['results/', dataname, '/'], 'dir')
        mkdir(['results/', dataname, '/']);
   end
   % save(['results/','classfication_results/' dataname, '/', dataname, '_', method, '_sum_results.mat'], 'accs', 'cpus')
end