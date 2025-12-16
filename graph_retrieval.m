function graph_retrieval_experiments()
    datanames = {'BZR', 'COX2','DHFR','PROTEINS_full', 'AIDS','STH_DHFR'};
    iterations_list = [1, 2, 3];
    hashDims = [100, 200, 300];
    k = 50; % Top-k for MAP computation
    opt_config.use_std_agg = 1;
    
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
    fprintf('ADAPTIVE OPTIMIZATION (Fast Path for sum-only) - Retrieval\n');
    fprintf('============================================================\n');
    fprintf('STD Aggregation:          %d\n', opt_config.use_std_agg);
    fprintf('============================================================\n\n');
    
    for d = 1:length(datanames)
        dataname = datanames{d};
        fprintf('============================================================\n');
        fprintf('Processing dataset: %s\n', dataname);
        fprintf('============================================================\n');
        
        % Load query and database data
        query_path = sprintf('data/%s_query/%s_query.mat', dataname, dataname);
        database_path = sprintf('data/%s_database/%s_database.mat', dataname, dataname);
        
        if ~exist(query_path, 'file') || ~exist(database_path, 'file')
            fprintf('Warning: Data files for %s not found, skipping...\n', dataname);
            continue;
        end
        
        try
            [query_graphs, query_labels] = load_matlab_graphs(query_path);
            [database_graphs, database_labels] = load_matlab_graphs(database_path);
        catch ME
            fprintf('Error loading %s: %s\n', dataname, ME.message);
            continue;
        end
        
        num_queries = length(query_graphs);
        num_database = length(database_graphs);
        
        fprintf('  Queries: %d, Database: %d\n', num_queries, num_database);
        
        % Initialize result matrices
        num_agg = size(agg_modes, 1);
        maps = zeros(length(iterations_list), length(hashDims), num_agg);
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
                    start_time = tic;
                    
                    % Check if we can use fast path: sum only + no std
                    use_fast_path = isequal(agg_mode, [1, 0, 0]) && ~opt_config.use_std_agg;
                    
                    if use_fast_path
                        % Fast path: optimized for sum-only
                        [query_fingerprints, hyperplanes] = generate_fingerprints_fast(...
                            query_graphs, iteration, hashDim);
                        database_fingerprints = generate_fingerprints_fast_with_hyperplanes(...
                            database_graphs, iteration, hashDim, hyperplanes);
                    else
                        % General path: flexible aggregations
                        [query_fingerprints, hyperplanes] = generate_fingerprints_general(...
                            query_graphs, iteration, hashDim, agg_mode, opt_config);
                        database_fingerprints = generate_fingerprints_general_with_hyperplanes(...
                            database_graphs, iteration, hashDim, agg_mode, opt_config, hyperplanes);
                    end
                    
                    % Compute similarity matrix using Hamming distance
                    gram_matrix = compute_hamming_similarity(...
                        query_fingerprints, database_fingerprints, iteration);
                    
                    % Compute MAP
                    [map_score, ~] = compute_map(gram_matrix, query_labels, database_labels, k);
                    
                    cpu_time = toc(start_time);
                    
                    maps(iIter, iHash, iAgg) = map_score;
                    cpus(iIter, iHash, iAgg) = cpu_time;
                    
                    fprintf('  Iter: %d, Dim: %d, MAP@%d: %.4f, Time: %.2fs\n', ...
                        iteration, hashDim, k, map_score, cpu_time);
                end
            end
        end
        

        results_dir = sprintf('results/retrieval/%s/', dataname);
        if ~exist(results_dir, 'dir')
            mkdir(results_dir);
        end
        
        % Save MAT
        results.maps = maps;
        results.cpus = cpus;
        results.iterations = iterations_list;
        results.hashDims = hashDims;
        results.agg_modes = agg_modes;
        results.opt_config = opt_config;
        results.k = k;
        results.dataname = dataname;
        
        results_path = [results_dir, dataname, '_hashwls_retrieval_results.mat'];
        save(results_path, '-struct', 'results');
        fprintf('\nResults saved to: %s\n', results_path);
        
        % ============================================================
        % SUMMARY
        % ============================================================
        [best_map, best_idx] = max(maps(:));
        [best_iter_idx, best_hash_idx, best_agg_idx] = ind2sub(size(maps), best_idx);
        
        % ============================================================
        % SAVE TXT
        % ============================================================
        txt_path = [results_dir, dataname, '_hashwls_retrieval_results.txt'];
        fid = fopen(txt_path, 'w');
        
        fprintf(fid, "Dataset: %s\n", dataname);
        fprintf(fid, "Method: hashwls_retrieval (Adaptive Optimization)\n");
        fprintf(fid, "k: %d\n\n", k);
        
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
                    fprintf(fid, "  hashDim=%3d  MAP=%.4f  Time=%.4f\n", ...
                        hashDim, maps(iIter, iHash, iAgg), cpus(iIter, iHash, iAgg));
                end
                fprintf(fid, "\n");
            end
            
            fprintf(fid, "--------------------------------------------------\n\n");
        end
        
        fprintf(fid, "**************** SUMMARY ****************\n");
        fprintf(fid, "Best MAP@%d = %.4f\n", k, best_map);
        fprintf(fid, "Best iteration = %d\n", iterations_list(best_iter_idx));
        fprintf(fid, "Best hashDim = %d\n", hashDims(best_hash_idx));
        fprintf(fid, "Best aggregation = %s\n\n", agg_modes{best_agg_idx, 2});
        fprintf(fid, "Mean MAP = %.4f\n", mean(maps(:)));
        fprintf(fid, "Total time = %.4f\n", sum(cpus(:)));
        
        fclose(fid);
        fprintf('TXT results saved to: %s\n', txt_path);
        
        fprintf('============================================================\n\n');
    end
    
    fprintf('All retrieval experiments completed!\n');
end

function [layer_fingerprints, hyperplanes] = generate_fingerprints_fast(graphs, iterations, hash_dims)
    % Optimized fast path for sum aggregation only
    % Returns both fingerprints and hyperplanes for database use
    
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


function layer_fingerprints = generate_fingerprints_fast_with_hyperplanes(graphs, iterations, hash_dims, hyperplanes)
    % Fast path for database graphs using pre-computed hyperplanes
    
    graph_num = length(graphs);
    
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


function [layer_fingerprints, hyperplanes] = generate_fingerprints_general(graphs, iterations, hash_dims, agg_mode, opt_config)
    % General path with flexible aggregation options
    % Returns both fingerprints and hyperplanes for database use
    
    graph_num = length(graphs);
    feature_size = size(graphs{1}.fv, 2);
    
    % Calculate feature dimension
    num_agg = sum(agg_mode);
    if opt_config.use_std_agg
        num_agg = num_agg + 1;
    end
    concat_dim = hash_dims * (1 + num_agg);
    
    % Pre-generate hyperplanes
    hyperplanes = cell(1, iterations);
    for r = 1:iterations
        if r == 1
            input_dim = feature_size;
        else
            input_dim = concat_dim;
        end
        hyperplanes{r} = randn(input_dim, hash_dims);
    end
    
    % Initialize storage
    layer_fingerprints = cell(1, iterations);
    for r = 1:iterations
        layer_fingerprints{r} = zeros(graph_num, concat_dim);
    end
    
    % Process each graph
    for i_graph = 1:graph_num
        node_num = size(graphs{i_graph}.am, 1);
        
        % First iteration
        features = double(graphs{i_graph}.fv);
        dot_results = features * hyperplanes{1};
        transformed_fingerprints = double(dot_results >= 0);
        
        concatenated_fingerprints = zeros(node_num, concat_dim);
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
            
            concatenated_fingerprints = zeros(node_num, concat_dim);
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


function layer_fingerprints = generate_fingerprints_general_with_hyperplanes(graphs, iterations, hash_dims, agg_mode, opt_config, hyperplanes)
    % General path for database graphs using pre-computed hyperplanes
    
    graph_num = length(graphs);
    
    % Calculate feature dimension
    num_agg = sum(agg_mode);
    if opt_config.use_std_agg
        num_agg = num_agg + 1;
    end
    concat_dim = hash_dims * (1 + num_agg);
    
    % Initialize storage
    layer_fingerprints = cell(1, iterations);
    for r = 1:iterations
        layer_fingerprints{r} = zeros(graph_num, concat_dim);
    end
    
    % Process each graph
    for i_graph = 1:graph_num
        node_num = size(graphs{i_graph}.am, 1);
        
        % First iteration
        features = double(graphs{i_graph}.fv);
        dot_results = features * hyperplanes{1};
        transformed_fingerprints = double(dot_results >= 0);
        
        concatenated_fingerprints = zeros(node_num, concat_dim);
        for i_node = 1:node_num
            neighbors = graphs{i_graph}.al{i_node};
            self_fp = transformed_fingerprints(i_node, :);
            
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
            
            concatenated_fingerprints = zeros(node_num, concat_dim);
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
    % Aggregate neighbor features based on configuration
    % Optimized: pre-allocate and build in one pass
    
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


function gram_matrix = compute_hamming_similarity(query_fingerprints, database_fingerprints, iterations)
    % Compute similarity matrix using Hamming distance
    % Optimized with vectorized operations
    
    num_queries = size(query_fingerprints{1}, 1);
    num_database = size(database_fingerprints{1}, 1);
    
    gram_matrix = zeros(num_queries, num_database);
    
    for r = 1:iterations
        query_fp = query_fingerprints{r};
        db_fp = database_fingerprints{r};
        dim = size(query_fp, 2);
        
        % Vectorized Hamming distance computation
        % For each query, compute distance to all database entries
        for iQ = 1:num_queries
            % Broadcast query fingerprint and compute XOR with all database
            diff = bsxfun(@ne, query_fp(iQ, :), db_fp);
            hamming_dist = sum(diff, 2)' / dim;
            gram_matrix(iQ, :) = gram_matrix(iQ, :) + (1 - hamming_dist);
        end
    end
end


function [map_score, map_list] = compute_map(gram_matrix, query_labels, database_labels, k)
    % Compute Mean Average Precision (MAP) for retrieval
    
    if nargin < 4
        k = 50;
    end
    
    num_queries = size(gram_matrix, 1);
    num_database = size(gram_matrix, 2);
    
    % Adjust k if database is smaller
    k = min(k, num_database);
    
    map_list = zeros(num_queries, 1);
    
    for iQ = 1:num_queries
        % Sort database graphs by similarity (descending)
        [~, sorted_indices] = sort(gram_matrix(iQ, :), 'descend');
        topk_indices = sorted_indices(1:k);
        
        % Compute Average Precision for this query
        curr = 0;
        precision = 0;
        
        for idx = 1:k
            if query_labels(iQ) == database_labels(topk_indices(idx))
                curr = curr + 1;
                precision = precision + curr / idx;
            end
        end
        
        if curr > 0
            precision = precision / curr;
        end
        
        map_list(iQ) = precision;
    end
    
    map_score = mean(map_list);
end


function [graphs, labels] = load_matlab_graphs(filepath)
    % Load graph data from MATLAB .mat file
    
    data = load(filepath);
    
    % Try different possible key names for graphs and labels
    graph_keys = {'graphs', 'query_graphs', 'database_graphs'};
    label_keys = {'labels', 'query_labels', 'database_labels'};
    
    graphs_mat = [];
    labels = [];
    
    for k = 1:length(graph_keys)
        if isfield(data, graph_keys{k})
            graphs_mat = data.(graph_keys{k});
            break;
        end
    end
    
    for k = 1:length(label_keys)
        if isfield(data, label_keys{k})
            labels = data.(label_keys{k});
            labels = labels(:); % Ensure column vector
            break;
        end
    end
    
    if isempty(graphs_mat)
        error('Could not find graph data in %s', filepath);
    end
    
    graphs = cell(size(graphs_mat, 1), 1);
    for i = 1:size(graphs_mat, 1)
        graphs{i}.fv = graphs_mat{i}.fv; % feature vectors
        graphs{i}.am = graphs_mat{i}.am; % adjacency matrix
        graphs{i}.al = cell(size(graphs_mat{i}.al, 1), 1); % adjacency list
        
        % Convert adjacency list (already in MATLAB 1-based indexing)
        for j = 1:size(graphs_mat{i}.al, 1)
            if ~isempty(graphs_mat{i}.al{j})
                neighbors = graphs_mat{i}.al{j}(:)';
                graphs{i}.al{j} = neighbors;
            else
                graphs{i}.al{j} = [];
            end
        end
    end
end
