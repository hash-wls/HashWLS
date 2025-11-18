function [acc,cpu] = hashwls(graphs, labels, iterations, hashDims)

graphNum = length(graphs); % number of graphs
featureSize = size(graphs{1,1}.fv, 2); % dimensions of feature vectors

hyperplanes = cell(1, iterations);
fingperints = cell(1, iterations);
for r = 1:iterations
    fingperints{r} = zeros(graphNum, hashDims*2);
    if r==1
        hyperplanes{r} = randn(featureSize, hashDims);
    else
        hyperplanes{r} = randn(hashDims*2, hashDims);
    end
end

tic;
for iGraph = 1:graphNum
    nodeNum = size(graphs{iGraph, 1}.am, 1);
    transformedFingerprints = zeros(nodeNum, hashDims);
    for iNode = 1:nodeNum
        features = graphs{iGraph, 1}.fv;
        feature = features(iNode, :)';
        transformedFingerprints(iNode, :) = dot(repmat(feature, 1, hashDims), hyperplanes{1})>=0;
    end

    concatenatedFingerprints = zeros(nodeNum, hashDims*2);
    for iNode = 1:nodeNum

        neighbors = graphs{iGraph, 1}.al{iNode,1};
        sumFingerprint = sum(transformedFingerprints(neighbors, :), 1);
        concatenatedFingerprints(iNode,:) = [transformedFingerprints(iNode,:), sumFingerprint];
    end
    
    fingperints{1}(iGraph,:) = sum(concatenatedFingerprints, 1);

    features = concatenatedFingerprints;
    for r = 2:iterations
        
        transformedFingerprints = zeros(nodeNum, hashDims);
        for iNode = 1:nodeNum

            feature = features(iNode, :)';
            transformedFingerprints(iNode, :) = dot(repmat(feature, 1, hashDims), hyperplanes{r})>=0;
        end

        concatenatedFingerprints = zeros(nodeNum, hashDims*2);
        for iNode = 1:nodeNum

            neighbors = graphs{iGraph, 1}.al{iNode,1};
            sumFingerprint = sum(transformedFingerprints(neighbors, :), 1);
            concatenatedFingerprints(iNode,:) = [transformedFingerprints(iNode,:), sumFingerprint];
        end
        features = concatenatedFingerprints;
        fingperints{r}(iGraph,:) = sum(concatenatedFingerprints, 1);

    end
end

K = zeros(graphNum);
for r = 1:iterations
    K = K+(1-squareform(pdist(fingperints{r},'hamming')));
end

acc = svmtrain(labels,[(1:size(K,1))' K], ['-t 4 -v 10 -q -c ', num2str(c)]);
cpu = toc;
