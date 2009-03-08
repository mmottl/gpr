function D = dist(x0,x1)

% dist: compute pairwise distance matrix from two column vectors x0 and x1

warning(['To speed up gradient calculation compile mex' ...
	   ' file dist.c'])

n0 = length(x0); n1 = length(x1);
D = repmat(x0,1,n1)-repmat(x1',n0,1);