function [isensors, w, Rinv, elapsetime] = F_sensor_SDPCN(Xorg,r1,r2,p)
%% notations
%{
    U   ... \mathbb{R}^{n*r} sensor candidate matrix
    Rnd ... \mathbb{R}^{n*n} nondiagonal component of noise covariance matrix
    a   ... R-Rnd = a*eye(n)
%}

    %pre
    [Psi, ~, ~, ~, ~, R, Q]=F_pre_SVD_corrnoise(Xorg,r1,r2);

    tic
    [n,r] = size(Psi);
  %  Rnd=Rlr-diag(Rlrdiag);
    a=0.002;  %positive constant
    Rnd=R-a*eye(n);
    Rnd=(Rnd+Rnd')/2;
    Q=0;  %for comparison

%%}

    solver_now = cvx_solver;

    cvx_solver('Mosek')
    Rndinv = inv(Rnd);
    HRH = Psi'*Rndinv*Psi;
%    C=inv(Q)+HRH;  %for comparison
    C=HRH;
    B=Rndinv*Psi;
% cvx_solver_settings('MSK_IPAR_NUM_THREADS',1)
    cvx_begin sdp
        variable w(n)
        variable Z(r,r) symmetric
        variable V(r,r) symmetric
        variable W(n,n) symmetric
        minimize( trace(Z) )
        subject to 
            [ (C-V+(C-V)')./2, eye(r1); ...
                 eye(r1)     , Z      ] >=0;  %eq. 19-1
            [   V,          B'; ...
                B,   Rndinv+diag(w)./a] >=0;       %eq. 19-2
            trace(W) <= p;
            diag(W) == w;
            [W,w;w',1] >= 0;
    cvx_end
    cvx_solver(solver_now);
    [~,isensors] = maxk(w,p);
    Rinv=inv(R(isensors,isensors));
 
    elapsetime=toc;
    disp(cvx_cputime-elapsetime)
end