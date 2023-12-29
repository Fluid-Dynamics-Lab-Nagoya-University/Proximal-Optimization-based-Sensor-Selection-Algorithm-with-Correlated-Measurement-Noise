function [FIM] = F_calc_FIM_WN(U,sensors)
    C = U(sensors,:);
    [p,r] = size(C);
    if p < r
        FIM = C*C';
    else
        FIM = C'*C;
    end
end