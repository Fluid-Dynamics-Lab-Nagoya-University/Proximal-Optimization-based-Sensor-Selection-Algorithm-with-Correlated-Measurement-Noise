function prox=F_proximal_operator(V,gamma,lambda,p)
    %prox=F_BlockSoftThreshold(V,gamma,lambda,n,m);
    %prox=F_BlockHardThreshold(V,gamma,lambda,n,m);
     prox=F_Lzero_const_BlockHardThreshold(V,p);

end
