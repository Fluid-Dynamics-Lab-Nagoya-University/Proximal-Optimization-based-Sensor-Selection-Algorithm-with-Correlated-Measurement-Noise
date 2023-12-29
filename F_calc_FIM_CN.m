function [FIM] = F_calc_FIM_CN(U,sensors,R)
    Rinv = inv(R(sensors,sensors));
    C = U(sensors,:);
    [p,r] = size(C);
    if p < r
        FIM = Rinv*C*C';
    else
        FIM = C'*Rinv*C;
    end
end