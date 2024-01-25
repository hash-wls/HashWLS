datanames={'BZR','COX2','DHFR','PROTEINS_full','AIDS','STH_DHFR'};
for i=1:length(datanames)
    dataname=datanames{i};
    load(['data/',dataname,'_query/',dataname,'_query.mat']);
    load(['data/',dataname,'_database/',dataname,'_database.mat']);
    iterations=1;
    query_fingerprints=generate_fingerprints(query_graphs,itration,hashDim);
    database_fingerprints=generate_fingerprints(database_graphs,iteration,hashDim);
    num_queries=length(query_graphs);
    num_database=length(database_graphs);
    hamming_distances = zeros(num_queries, num_database);
    gram_matrix = zeros(num_queries, num_database);
    
    for r = 1:iteration
        for i = 1:num_queries
            for j = 1:num_database
                hamming_distances(i, j) = sum(query_fingerprints{r}(i,:) ~= database_fingerprints{r}(j,:));
            end
        end
        gram_matrix = gram_matrix + (1 - hamming_distances);
    end
    map_list = zeros(1, num_queries);
   
    
    k=50;
    for iQ = 1:num_queries
        [~, topk_indices] = sort(gram_matrix(iQ, :), 'descend');
            topk_indices = topk_indices(1:k);
            
            % Compute MAP
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
                p=curr/k;
    
            end
    
            map_list(iQ) = precision;
            
    
    end
    map = mean(map_list);
end
save(['results/', 'retrieval_results/',dataname, '/', dataname, '.mat'], 'map')