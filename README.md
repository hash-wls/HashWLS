The datasets and source code of #WLS are for Fast Hash Sketching for Graphs with Continuous Attributes.

The steps of running the experiments:

Preparation work
```
1. compile svmtrain in libsvm-matlab by referring to libsvm-matlab/matlab/README
2. put the mex file in the same directory as 'main.m'
3. cd HashWLS
```   
Graph Classification

```
run graph_classification.m
```

Graph Retrieval
```
run graph_retrieval.m
```

Parameter Sensitivity
```
run draw_parameters.m
```


Scalability
```
run scal_main.m
```

