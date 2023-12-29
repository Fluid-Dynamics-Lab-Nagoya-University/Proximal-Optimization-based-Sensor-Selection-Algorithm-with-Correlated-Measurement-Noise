function [isensors,D,Rinv,itr,elapsetime]=F_sensor_ADMMCN(Xorg,r1,r2,p)
    %params
    thre       = 1e-4;
    gamma_init = 1.0;
    gamma_lim  = 2e-5;
    eta        = 0.99;
    tolX       = 8e-6;
    tol_constraint = 1e-3;
    gamma_intval   = 5000;
    stdout_intval  = 5000;
    lambda=1; %only for BST amd BHT

    %initia
    epsiX = 1e6;
    itr = 0;
    epsi_constraint = 200;
    gamma=gamma_init;

    %pre
    [Psi, Ur, Sr, dSdiag, Rdiag, R, ~]=F_pre_SVD_corrnoise(Xorg,r1,r2);

tic
    %normalize
    Psinorm=diag(Rdiag.^(-0.5))*Psi;
    Psi=Psinorm;
    dSdiag=Rdiag.^(-0.5).*dSdiag.*Rdiag.^(-0.5);
    
    Ur=diag(Rdiag.^(-0.5))*Ur;
    
    %initial solution
    X=pinv(Psi)';
    A=Psi';
    Is=eye(r1,r1);
    Z1 = zeros(size(X));
%    Z1 = X;
    Z2 = Is;
    Y1 = zeros(size(X));
    Y2 = zeros(size(Is));
%    Y1 = Y1 +     X - Z1;
%    Y2 = Y2 + A * X - Z2;

    %pre
    dSIinvdiag=1./(2*dSdiag+1/gamma);
    AdSIiAt=inv(Is+(A*diag(dSIinvdiag))*1/gamma*A');

    dSIinvr1c=repmat(dSIinvdiag',r1,1);
    dSIinvr1r=repmat(dSIinvdiag,1,r1);
    dSIinvr1r2c=repmat(dSIinvdiag',r2-r1,1);
    dSIinvr1r2r=repmat(dSIinvdiag,1,r2-r1);

    %main
    while epsiX > tolX || epsi_constraint > tol_constraint
        Xpre = X;
        tmpZY = 1 / gamma * ( ( Z1 - Y1 )  + A' * ( Z2 - Y2 ) );

        %step 1
        X=(  (dSIinvr1r.*tmpZY  )-( (dSIinvr1r.*(1/gamma).*A')*(AdSIiAt* ((A.*dSIinvr1c)*tmpZY)) )  ) ...  %ARinv  o
        - (  dSIinvr1r2r.*Ur - ((dSIinvr1r.*(1/gamma).*A'))*(AdSIiAt*((A.*dSIinvr1c)*Ur))  ) ...   % ARinv*Ur  o
        *(( ...
        inv( inv(Sr^2) ...
        + (  (Ur'.*dSIinvr1r2c)*Ur - (Ur'*(dSIinvr1r.*(1/gamma).*A'))*(AdSIiAt*((A.*dSIinvr1c)*Ur))  ) ...  %Ur'*ARinv*Ur  o
        ) ...
        * (  Ur'.*dSIinvr1r2c - (Ur'*(dSIinvr1r.*(1/gamma).*A'))*(AdSIiAt*(A.*dSIinvr1c))   )...   %n*n(Ur'*ARinv) o
        )*tmpZY);

        %step 2
        Z1 = F_proximal_operator( (X + Y1)', gamma, lambda, p);   %A(m,n) X(n,m)

        Z1 = Z1';
        Z2 = Is;  %constraint  AX=I

        %step 3
        Y1 = Y1 +     X - Z1;
        Y2 = Y2 + A * X - Z2;
        epsiX=norm(X-Xpre,'fro');

        itr=itr+1;
        if mod(itr,stdout_intval) == 0
            if mod(itr,gamma_intval) == 0 && gamma > gamma_lim
                gamma=gamma*eta;
                dSIinvdiag=1./(2*dSdiag+1/gamma);
                dSIinvr1c=repmat(dSIinvdiag',r1,1);
                dSIinvr1r=repmat(dSIinvdiag,1,r1);
                dSIinvr1r2c=repmat(dSIinvdiag',r2-r1,1);
                dSIinvr1r2r=repmat(dSIinvdiag,1,r2-r1);
                AdSIiAt=inv(Is+(A.*dSIinvr1c)*1/gamma*A');
            end
            epsi_constraint   = norm(A*X-Is,'fro');
            trnorm = trace(X'*X);
            disp(['step: ' num2str(itr,'%8d') ', gamma: ' num2str(gamma) ', ||Xpre-Xnew||_F: ' num2str(epsiX) ', ||AX-I||_F: ' num2str(epsi_constraint) ', tr(X^TX): ' num2str(trnorm)])
        end
        %toc
    end
    bnorm = sum(X.*X,2);
%{
    Xnorm=zeros(n,1);
    for i=1:n
        Xnorm(i)=norm(X(i,:),2);
    end
%}
    isensors = find(bnorm>thre);
    psel=size(isensors,1);
    disp(['selected sensor: ' num2str(psel)])

    D=X';
    Rinv=inv(R(isensors,isensors));
    
    elapsetime=toc;
end


