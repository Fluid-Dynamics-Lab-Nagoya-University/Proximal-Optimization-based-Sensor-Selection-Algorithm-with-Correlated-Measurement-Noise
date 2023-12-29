function [Psi, Ur, Sr, dSdiag, Rdiag, R, Q]=F_pre_SVD_corrnoise(Xorg,r1,r2)
    [n,m] = size(Xorg);
    [Uorg,Sorg,Vorg] = svd(Xorg,'econ');
    [~,r]=size(Sorg);
    Psi = Uorg(:,1:r1);
    
    if isempty(r2) == 0
        Ur=Uorg(:,r1+1:r2);
        Sr=Sorg(r1+1:r2,r1+1:r2);
        Xorglr=Uorg(:,1:r1)*Sorg(1:r1,1:r1)*Vorg(:,1:r1)';
        Rlrdiag=diag(Ur*Sr*Sr*Ur');
        dSdiag=zeros(n,1);
        for i=1:n
            tmp=0;
            for j=1:r
                tmp=tmp+(Xorg(i,j)-Xorglr(i,j))^2;
            end
            dSdiag(i)=tmp-Rlrdiag(i);
        end
        Rdiag=Rlrdiag+dSdiag;
        R=1/r*(Ur*Sr*Sr*Ur'+diag(dSdiag));
    else
        Ur=1;
        Sr=1;
        dSdiag=1;
        Rdiag=1;
        R=1;
    end
    Qdiag=1/r*diag(Sorg(1:r1,1:r1)*Sorg(1:r1,1:r1));
    Q=diag(Qdiag);
end