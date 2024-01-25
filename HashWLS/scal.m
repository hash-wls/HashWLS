function [acc,time]= scal(graphs, labels, iterations,siz, hashDims)

graphNum = length(graphs); % number of graphs
random_numbers = randperm(20000, siz);
label=zeros(siz,1);
featureSize = size(graphs{1,1}.fv, 2); % dimensions of feature vectors

hyperplanes = cell(1, iterations);
fingperints = cell(1, iterations);
for r = 1:iterations
    fingperints{r} = zeros(siz, hashDims*2);
    if r==1
        hyperplanes{r} = randn(featureSize, hashDims);
    else
        hyperplanes{r} = randn(hashDims*2, hashDims);
    end
end

tic;
runtime=cputime;
for i = 1:siz
    iGraph = random_numbers(i);
    label(i)=labels(iGraph);
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
    
    fingperints{1}(i,:) = sum(concatenatedFingerprints, 1);

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
        fingperints{r}(i,:) = sum(concatenatedFingerprints, 1);

    end
end
K = zeros(siz);
for r = 1:iterations
    K = K+(1-squareform(pdist(fingperints{r},'hamming')));
end

acc = svmtrain(label,[(1:size(K,1))' K],'-t 4 -v 10 -q');
time = toc; 
