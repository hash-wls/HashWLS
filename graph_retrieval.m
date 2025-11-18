function graph_retrieval_experiments()
    % Main function to run graph retrieval experiments
    
    % datanames = {'BZR'}; % {'BZR', 'COX2', 'DHFR', 'PROTEINS_full', 'AIDS', 'STH_DHFR'}
    datanames = {'BZR', 'COX2', 'DHFR', 'PROTEINS_full', 'AIDS', 'STH_DHFR'};
    iterations = [1, 2, 3, 4, 5]; % Number of WL iterations
    hash_dims = [50, 100, 150, 200, 250, 300]; % Hash dimensions
    k = 50; % Top-k for MAP computation
    
    for d = 1:length(datanames)
        dataname = datanames{d};
        fprintf('Processing dataset: %s\n', dataname);
        
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
        fprintf('  Iterations: [%s], Hash dims: [%s]\n', num2str(iterations), num2str(hash_dims));
        
        for h = 1:length(hash_dims)
            hash_dim = hash_dims(h);
            for it = 1:length(iterations)
                iteration = iterations(it);
                start_time = tic;
                
                % Set random seed for reproducibility
                rng(42);
                
                % Generate fingerprints
                query_fingerprints = generate_fingerprints_retrieval(...
                    query_graphs, iteration, hash_dim, 22);
                database_fingerprints = generate_fingerprints_retrieval(...
                    database_graphs, iteration, hash_dim, 22);
                
                % Compute Gram matrix (similarity matrix)
                num_queries_local = size(query_fingerprints, 1);
                num_database_local = size(database_fingerprints, 1);
                gram_matrix = zeros(num_queries_local, num_database_local);
                
                % Use cosine similarity for multi-scale features
                query_norm = query_fingerprints ./ (vecnorm(query_fingerprints, 2, 2) + 1e-10);
                db_norm = database_fingerprints ./ (vecnorm(database_fingerprints, 2, 2) + 1e-10);
                gram_matrix = query_norm * db_norm';
                
                % Compute MAP
                [map_score, map_list] = compute_map(...
                    gram_matrix, query_labels, database_labels, k);
                
                cpu_time = toc(start_time);
                fprintf('%s, Hash dim: %d, Iteration: %d\n', dataname, hash_dim, iteration);
                fprintf('  MAP@%d: %.4f\n', k, map_score);
                fprintf('  Time: %.2fs\n', cpu_time);
            end
        end
        
        % Save results
        % results_dir = sprintf('results/retrieval_results/%s', dataname);
        if ~exist(results_dir, 'dir')
            mkdir(results_dir);
        end
        
        results.map = map_score;
        results.map_list = map_list;
        results.gram_matrix = gram_matrix;
        results.iterations = iterations;
        results.hash_dims = hash_dims;
        results.k = k;
        results.cpu_time = cpu_time;
        
        % Save as .mat file
        % results_path = sprintf('%s/%s.mat', results_dir, dataname);
        % save(results_path, '-struct', 'results');
        
        % fprintf('  Results saved to %s\n\n', results_dir);
    end
    
    fprintf('All experiments completed!\n');
end

function result = matlab_dot_columnwise(A, B)
    % Mimics MATLAB's dot function for matrices.
    % MATLAB's dot(A,B) computes dot product of corresponding columns.
    result = sum(A .* B, 1);
end

function all_layer_fingerprints = generate_fingerprints_retrieval(graphs, iterations, hash_dims, seed)
    % Improved fingerprint generation with multiple enhancements
    
    if nargin < 4
        seed = 22;
    end
    
    graph_num = length(graphs);
    feature_size = size(graphs{1}.fv, 2);
    
    hyperplanes = {};
    all_layer_fingerprints = [];
    
    % Fix: Use different seeds for different iterations
    rng(seed);
    for r = 1:iterations
        if r == 1
            hyperplanes{r} = randn(feature_size, hash_dims);
        else
            % Use progressively larger feature spaces
            hyperplanes{r} = randn(hash_dims * 3, hash_dims);
        end
    end
    
    % Process each graph
    for i_graph = 1:graph_num
        node_num = size(graphs{i_graph}.am, 1);
        layer_features_list = [];
        
        % First iteration with improved aggregation
        features = double(graphs{i_graph}.fv);
        transformed_fingerprints = zeros(node_num, hash_dims);
        
        for i_node = 1:node_num
            feature = features(i_node, :);
            feature_replicated = repmat(feature', 1, hash_dims);
            dot_result = matlab_dot_columnwise(feature_replicated, hyperplanes{1});
            transformed_fingerprints(i_node, :) = double(dot_result >= 0);
        end
        
        % Improved neighbor aggregation: sum, mean, max
        concatenated_fingerprints = zeros(node_num, hash_dims * 3);
        for i_node = 1:node_num
            neighbors = graphs{i_graph}.al{i_node};
            if ~isempty(neighbors)
                neighbor_features = transformed_fingerprints(neighbors, :);
                sum_fingerprint = sum(neighbor_features, 1);
                mean_fingerprint = mean(neighbor_features, 1);
                max_fingerprint = max(neighbor_features, [], 1);
            else
                sum_fingerprint = zeros(1, hash_dims);
                mean_fingerprint = zeros(1, hash_dims);
                max_fingerprint = zeros(1, hash_dims);
            end
            
            concatenated_fingerprints(i_node, :) = [...
                sum_fingerprint, mean_fingerprint, max_fingerprint];
        end
        
        % Store layer 0 features
        layer_features_list = [layer_features_list; sum(concatenated_fingerprints, 1)];
        
        % Subsequent iterations
        features = concatenated_fingerprints;
        for r = 2:iterations
            transformed_fingerprints = zeros(node_num, hash_dims);
            
            for i_node = 1:node_num
                feature = features(i_node, :);
                feature_replicated = repmat(feature', 1, hash_dims);
                dot_result = matlab_dot_columnwise(feature_replicated, hyperplanes{r});
                transformed_fingerprints(i_node, :) = double(dot_result >= 0);
            end
            
            concatenated_fingerprints = zeros(node_num, hash_dims * 3);
            for i_node = 1:node_num
                neighbors = graphs{i_graph}.al{i_node};
                if ~isempty(neighbors)
                    neighbor_features = transformed_fingerprints(neighbors, :);
                    sum_fingerprint = sum(neighbor_features, 1);
                    mean_fingerprint = mean(neighbor_features, 1);
                    max_fingerprint = max(neighbor_features, [], 1);
                else
                    sum_fingerprint = zeros(1, hash_dims);
                    mean_fingerprint = zeros(1, hash_dims);
                    max_fingerprint = zeros(1, hash_dims);
                end
                
                concatenated_fingerprints(i_node, :) = [...
                    sum_fingerprint, mean_fingerprint, max_fingerprint];
            end
            
            features = concatenated_fingerprints;
            layer_features_list = [layer_features_list; sum(concatenated_fingerprints, 1)];
        end
        
        % Multi-scale: concatenate all layers
        all_features = [];
        for l = 1:size(layer_features_list, 1)
            all_features = [all_features, layer_features_list(l, :)];
        end
        all_layer_fingerprints = [all_layer_fingerprints; all_features];
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
    
    % Convert MATLAB cell array to cell array of structures
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
