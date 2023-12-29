function [sensors,K,itr,time]=F_sensor_ADMMWN(Xorg,r1,p)
    %params
    thre       = 1e-5;
    gamma_init = 0.5;
    gamma_lim  = 2e-5;
    eta        = 0.99;
    tolX       = 8e-6;
    tol_constraint = 1e-3;
    gamma_intval   = 5000;
    stdout_intval  =  500;
    lambda=1000; %only for BST amd BHT

    %initia
    epsiX = 1e6;
    itr = 0;
    epsi_constraint = 200;
    gamma=gamma_init;

    %pre
    [Psi, ~, ~, ~, ~, ~, ~]=F_pre_SVD_corrnoise(Xorg,r1,[]);

    tic
    X = pinv(Psi)';
    A = Psi';
    I = eye(r1);

    Z1 = zeros(size(X));
    Z2 = I;
    Y1 = zeros(size(X));
    Y2 = zeros(size(I));
    
    gami  = 1/(2 + 1/gamma);
    AAinv = inv( I + gami/gamma * A * A' );
    
    while epsiX > tolX || epsi_constraint > tol_constraint
        Xpre = X;
        tmpZY = 1 / gamma * ( ( Z1 - Y1 )  + A' * ( Z2 - Y2 ) );
        X =  gami * tmpZY  - gami * gami / gamma * A' * ( AAinv * ( A * tmpZY ) );
    
        Z1 = F_proximal_operator( (X + Y1)', gamma, lambda, p);
        Z1 = Z1';
        Z2 = I;  %constraint AX=I
    
        Y1 = Y1 +     X - Z1;
        Y2 = Y2 + A * X - Z2;
    
        epsiX=norm(X-Xpre,'fro')/norm(X,'fro');
        itr=itr+1;
        if mod(itr,stdout_intval) == 0
            if mod(itr,gamma_intval) == 0 && gamma > gamma_lim
                gamma=gamma*eta;
                gami  = 1/(2 + 1/gamma);
                AAinv = inv( I + gami/gamma * A * A' );
            end
            epsi_constraint   = norm(A*X-I,'fro');
            trnorm = trace(X'*X);
            disp(['step: ' num2str(itr,'%8d') ', gamma: ' num2str(gamma) ', ||Xpre-Xnew||_F: ' num2str(epsiX) ', ||AX-I||_F: ' num2str(epsi_constraint) ', tr(X^TX): ' num2str(trnorm)])
        end
    end
    bnorm = sum(X.*X,2);
    sensors = find(bnorm>thre);
    psel=size(sensors,1);
    disp(['selected sensor: ' num2str(psel)])

    K=X';

    time=toc;
end

