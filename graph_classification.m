clear all
clc

datanames = {'BZR', 'COX2','DHFR','PROTEINS_full', 'AIDS','STH_DHFR'};
iterations_list = [1,2,3];
hashDims = [100,200,300];
method = 'hashwls_adaptive';
c = 10.0; 

opt_config.use_std_agg = 1; 
opt_config.test_c_values = 0;
% Aggregation mode combinations for pooling
agg_modes = {
    [1, 0, 0], 'sum';
    [0, 1, 0], 'mean';
    [0, 0, 1], 'max';
    [1, 1, 0], 'sum_mean';
    [1, 0, 1], 'sum_max';
    [0, 1, 1], 'mean_max';
    [1, 1, 1], 'sum_mean_max';
};

% Print optimization settings
fprintf('============================================================\n');
fprintf('ADAPTIVE OPTIMIZATION (Fast Path for sum-only)\n');
fprintf('============================================================\n');
fprintf('STD Aggregation:          %d\n', opt_config.use_std_agg);
fprintf('Test C Values:            %d\n', opt_config.test_c_values);
fprintf('============================================================\n\n');

for idataname = 1:length(datanames)
    
    dataname = datanames{idataname};
    fprintf('============================================================\n');
    fprintf('Processing dataset: %s\n', dataname);
    fprintf('============================================================\n');
    
    % Load dataset
    data_path = ['data/', dataname, '/', dataname, '.mat'];
    if ~exist(data_path, 'file')
        fprintf('Warning: Data file not found: %s\n', data_path);
        fprintf('Skipping dataset %s\n\n', dataname);
        continue;
    end
    
    load(data_path);
    
    % Initialize result matrices
    num_agg = size(agg_modes, 1);
    accs = zeros(length(iterations_list), length(hashDims), num_agg);
    cpus = zeros(length(iterations_list), length(hashDims), num_agg);

    for iAgg = 1:num_agg
        agg_mode = agg_modes{iAgg, 1};
        agg_name = agg_modes{iAgg, 2};
        
        fprintf('\n--- Aggregation: %s ---\n', agg_name);
        
        for iIter = 1:length(iterations_list)
            iteration = iterations_list(iIter);
            
            for iHash = 1:length(hashDims)
                hashDim = hashDims(iHash);
                
                rng(42);
                
                [accs(iIter, iHash, iAgg), cpus(iIter, iHash, iAgg)] = ...
                    hashwls_adaptive(graphs, labels, iteration, hashDim, c, agg_mode, opt_config);
                
                fprintf('  Iter: %d, Dim: %d, Acc: %.2f%%, Time: %.2fs\n', ...
                    iteration, hashDim, accs(iIter, iHash, iAgg), cpus(iIter, iHash, iAgg));
            end
        end
    end
    

    accs_mean = mean(accs, 3);
    cpus_mean = mean(cpus, 3);

    results_dir = ['results/', dataname, '/'];
    if ~exist(results_dir, 'dir')
        mkdir(results_dir);
    end

    results.accs = accs;
    results.cpus = cpus;
    results.accs_mean = accs_mean;
    results.cpus_mean = cpus_mean;
    results.iterations = iterations_list;
    results.hashDims = hashDims;
    results.agg_modes = agg_modes;
    results.method = method;
    results.c = c;
    results.dataname = dataname;
    results.opt_config = opt_config;

    results_path = [results_dir, dataname, '_', method, '_results.mat'];
    save(results_path, '-struct', 'results');

    fprintf('\nResults saved to: %s\n', results_path);
    [best_acc, best_idx] = max(accs(:));
    [best_iter_idx, best_hash_idx, best_agg_idx] = ind2sub(size(accs), best_idx);

    txt_path = [results_dir, dataname, '_', method, '_results.txt'];
    fid = fopen(txt_path, 'w');

    fprintf(fid, "Dataset: %s\n", dataname);
    fprintf(fid, "Method: %s (Hamming Linear Kernel + Adaptive Optimization)\n", method);
    fprintf(fid, "C: %d\n\n", c);
    
    % Save optimization settings
    fprintf(fid, "OPTIMIZATION SETTINGS:\n");
    fprintf(fid, "STD Aggregation: %d\n\n", opt_config.use_std_agg);

    for iAgg = 1:num_agg
        agg_name = agg_modes{iAgg, 2};
        fprintf(fid, "================ Aggregation: %s ================\n", agg_name);

        for iIter = 1:length(iterations_list)
            iteration = iterations_list(iIter);
            fprintf(fid, "Iteration = %d\n", iteration);

            for iHash = 1:length(hashDims)
                hashDim = hashDims(iHash);
                fprintf(fid, "  hashDim=%3d  Acc=%.4f  Time=%.4f\n", ...
                    hashDim, accs(iIter, iHash, iAgg), cpus(iIter, iHash, iAgg));
            end
            fprintf(fid, "\n");
        end

        fprintf(fid, "--------------------------------------------------\n\n");
    end

    fprintf(fid, "**************** SUMMARY ****************\n");
    fprintf(fid, "Best accuracy = %.4f\n", best_acc);
    fprintf(fid, "Best iteration = %d\n", iterations_list(best_iter_idx));
    fprintf(fid, "Best hashDim = %d\n", hashDims(best_hash_idx));
    fprintf(fid, "Best aggregation = %s\n\n", agg_modes{best_agg_idx, 2});
    fprintf(fid, "Mean accuracy = %.4f\n", mean(accs(:)));
    fprintf(fid, "Total time    = %.4f\n", sum(cpus(:)));
    fclose(fid);
    fprintf('TXT results saved to: %s\n', txt_path);
    fprintf('============================================================\n\n');

end

fprintf('All classification experiments completed!\n');



function [acc, cpu] = hashwls_adaptive(graphs, labels, iterations, hashDims, c, agg_mode, opt_config)
if nargin < 5
    c = 1;
end

if nargin < 6
    agg_mode = [1, 1, 1];  % Default: sum+mean+max
end

if nargin < 7
    opt_config.use_std_agg = 0;
    opt_config.test_c_values = 0;
end

tic;

% Check if we can use fast path: sum only + no std
use_fast_path = isequal(agg_mode, [1, 0, 0]) && ~opt_config.use_std_agg;

if use_fast_path
    % Fast path: optimized for sum-only
    fingerprints_by_layer = generate_fingerprints_fast(graphs, iterations, hashDims);
else
    % General path: flexible but slightly slower
    fingerprints_by_layer = generate_fingerprints_general(graphs, iterations, hashDims, agg_mode, opt_config);
end

graphNum = length(graphs);

% Compute Hamming distance kernel
K = zeros(graphNum);
for r = 1:iterations
    layer_fingerprints = fingerprints_by_layer{r};
    K = K + (1 - squareform(pdist(layer_fingerprints, 'hamming')));
end

% SVM with precomputed kernel
acc = svmtrain(labels, [(1:size(K,1))' K], ['-t 4 -v 10 -q -c ', num2str(c)]);
cpu = toc;

% Test different C values if enabled
if opt_config.test_c_values
    test_c_values(labels, K, c, acc);
end

end


function test_c_values(labels, K, default_c, default_acc)
    fprintf('    [C Test] Default(%.1f): %.2f%% | ', default_c, default_acc);
    c_list = [1, 10, 100];
    
    for c_val = c_list
        acc_test = svmtrain(labels, [(1:size(K,1))' K], ['-t 4 -v 10 -q -c ', num2str(c_val)]);
        fprintf('C=%.1f: %.2f%% | ', c_val, acc_test);
    end
    fprintf('\n');
    drawnow;
end


function layer_fingerprints = generate_fingerprints_fast(graphs, iterations, hash_dims)
    % Optimized fast path for sum aggregation only
    % No conditional checks, direct execution
    
    graph_num = length(graphs);
    feature_size = size(graphs{1}.fv, 2);
    
    % Pre-generate hyperplanes
    hyperplanes = cell(1, iterations);
    for r = 1:iterations
        if r == 1
            input_dim = feature_size;
        else
            input_dim = hash_dims * 2;  % [self, sum]
        end
        hyperplanes{r} = randn(input_dim, hash_dims);
    end
    
    % Initialize storage
    layer_fingerprints = cell(1, iterations);
    for r = 1:iterations
        layer_fingerprints{r} = zeros(graph_num, hash_dims * 2);
    end
    
    % Process each graph
    for i_graph = 1:graph_num
        node_num = size(graphs{i_graph}.am, 1);
        
        % First iteration
        features = double(graphs{i_graph}.fv);
        dot_results = features * hyperplanes{1};
        transformed_fingerprints = double(dot_results >= 0);
        
        concatenated_fingerprints = zeros(node_num, hash_dims * 2);
        for i_node = 1:node_num
            neighbors = graphs{i_graph}.al{i_node};
            % MATLAB handles empty arrays automatically
            concatenated_fingerprints(i_node, :) = [...
                transformed_fingerprints(i_node, :), ...
                sum(transformed_fingerprints(neighbors, :), 1)];
        end
        
        layer_fingerprints{1}(i_graph, :) = sum(concatenated_fingerprints, 1);
        
        % Subsequent iterations
        features = concatenated_fingerprints;
        for r = 2:iterations
            dot_results = features * hyperplanes{r};
            transformed_fingerprints = double(dot_results >= 0);
            
            concatenated_fingerprints = zeros(node_num, hash_dims * 2);
            for i_node = 1:node_num
                neighbors = graphs{i_graph}.al{i_node};
                concatenated_fingerprints(i_node, :) = [...
                    transformed_fingerprints(i_node, :), ...
                    sum(transformed_fingerprints(neighbors, :), 1)];
            end
            
            layer_fingerprints{r}(i_graph, :) = sum(concatenated_fingerprints, 1);
            features = concatenated_fingerprints;
        end
    end
end



function layer_fingerprints = generate_fingerprints_general(graphs, iterations, hash_dims, agg_mode, opt_config)
    % General path with flexible aggregation options
    % Still optimized with batch projection
    
    graph_num = length(graphs);
    feature_size = size(graphs{1}.fv, 2);
    
    % Calculate feature dimension
    num_agg = sum(agg_mode);
    if opt_config.use_std_agg
        num_agg = num_agg + 1;
    end
    
    % Pre-generate hyperplanes
    hyperplanes = cell(1, iterations);
    for r = 1:iterations
        if r == 1
            input_dim = feature_size;
        else
            input_dim = hash_dims * (1 + num_agg);
        end
        hyperplanes{r} = randn(input_dim, hash_dims);
    end
    
    % Initialize storage
    layer_fingerprints = cell(1, iterations);
    for r = 1:iterations
        layer_fingerprints{r} = zeros(graph_num, hash_dims * (1 + num_agg));
    end
    
    % Process each graph
    for i_graph = 1:graph_num
        node_num = size(graphs{i_graph}.am, 1);
        
        % First iteration
        features = double(graphs{i_graph}.fv);
        dot_results = features * hyperplanes{1};
        transformed_fingerprints = double(dot_results >= 0);
        
        concatenated_fingerprints = zeros(node_num, hash_dims * (1 + num_agg));
        for i_node = 1:node_num
            neighbors = graphs{i_graph}.al{i_node};
            self_fp = transformed_fingerprints(i_node, :);
            
            % Aggregate neighbor features
            agg_features = aggregate_neighbors(transformed_fingerprints, neighbors, ...
                                              hash_dims, agg_mode, opt_config);
            
            concatenated_fingerprints(i_node, :) = [self_fp, agg_features];
        end
        
        layer_fingerprints{1}(i_graph, :) = sum(concatenated_fingerprints, 1);
        
        % Subsequent iterations
        features = concatenated_fingerprints;
        for r = 2:iterations
            dot_results = features * hyperplanes{r};
            transformed_fingerprints = double(dot_results >= 0);
            
            concatenated_fingerprints = zeros(node_num, hash_dims * (1 + num_agg));
            for i_node = 1:node_num
                neighbors = graphs{i_graph}.al{i_node};
                self_fp = transformed_fingerprints(i_node, :);
                
                agg_features = aggregate_neighbors(transformed_fingerprints, neighbors, ...
                                                  hash_dims, agg_mode, opt_config);
                
                concatenated_fingerprints(i_node, :) = [self_fp, agg_features];
            end
            
            layer_fingerprints{r}(i_graph, :) = sum(concatenated_fingerprints, 1);
            features = concatenated_fingerprints;
        end
    end
end


function agg_features = aggregate_neighbors(transformed_fingerprints, neighbors, hash_dims, agg_mode, opt_config)
    if isempty(neighbors)
        % Empty neighbors: return zeros
        num_agg = sum(agg_mode);
        if opt_config.use_std_agg
            num_agg = num_agg + 1;
        end
        agg_features = zeros(1, hash_dims * num_agg);
        return;
    end
    
    neighbor_features = transformed_fingerprints(neighbors, :);
    agg_features = [];
    
    % Apply selected aggregations
    if agg_mode(1)  % sum
        agg_features = [agg_features, sum(neighbor_features, 1)];
    end
    if agg_mode(2)  % mean
        agg_features = [agg_features, mean(neighbor_features, 1)];
    end
    if agg_mode(3)  % max
        agg_features = [agg_features, max(neighbor_features, [], 1)];
    end
    if opt_config.use_std_agg  % std
        agg_features = [agg_features, std(neighbor_features, 0, 1)];
    end
end
