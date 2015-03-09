function [c,b] = cvx_foopsi(y,b,sn,b_lb,G,w)

% implementation of constrained foopsi in CVX
% Written by Eftychios Pnevmatikakis

    if isempty(b)
        bas_est = 1;
    else
        bas_est = 0;
    end
    T = length(y);
    cvx_begin quiet
        variable c2(T)
        if bas_est; variable b; end
        minimize(w'*(G*c2))
        subject to
            G*c2>=0;
            norm(y-c2-b)<=sqrt(T)*sn;
            if bas_est; b>=b_lb; end
    cvx_end
    if strcmpi(cvx_status,'Infeasible');
        warning('Problem is infeasible, perhaps noise is underestimated. Projecting original data to the cone induced by the constraints. \n');
        cvx_begin quiet
            variable c2(T)
            if bas_est; variable b; end
            minimize(norm(y-c2-b))
            subject to
                G*c2>=0;
                if bas_est; b>=b_lb; end
        cvx_end
        sn = cvx_optval/sqrt(T);
    end
    c = c2;
end