function V=F_Lzero_const_BlockHardThreshold(V,p)
    bnorm = sum(V.*V,1);
    [th,~] = maxk(bnorm,p);
    V = V.*(bnorm>=th(p));
end

%{
    [~,n]=size(V);
    %bnorm=gpuArray(zeros(n,1));
    bnorm=zeros(n,1);
    %sv=gather(V);
    sv=V;
    lp=2;
   % tic
    for j=1:n
        %bnorm(j)=norm(V(:,j),lp);
        bnorm(j)=norm(sv(:,j),lp);
    end
   % toc
    tmp=bnorm;
    for pp=1:p
        [th,k]=max(tmp);
        tmp(k)=0;
    end
    %sv=V;
    V(:,find(bnorm<th))=0;
end

%}
