n=2000000;
est_matrix=zeros(n,50);

estimatedRatiotrn=rand(n);
iter=3




@time est_matrix[:,iter+1]=estimatedRatiotrn[:];
@time write_column!(est_matrix,iter+1,estimatedRatiotrn)
