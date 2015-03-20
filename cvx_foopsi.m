function [c,b,c1] = cvx_foopsi(y,b,c1,sn,b_lb,g,w)

% implementation of constrained foopsi in CVX
% Written by Eftychios Pnevmatikakis

    if isempty(b)
        bas_est = 1;
    else
        bas_est = 0;
    end
    if isempty(c1)
        c1_est = 1;
    else
        c1_est = 0;
    end
    gd = max(roots([1,-g(:)']));
    T = length(y);
    G = spdiags(ones(T,1)*[-g(end:-1:1)',1],-length(g):0,T,T);
    cvx_begin quiet
        variable c2(T)
        if bas_est; variable b; end
        if c1_est; variable c1; end
        minimize(w'*(G*c2))
        subject to
            G*c2>=0;
            norm(y-c2-b-c1*gd.^((0:T-1)'))<=sqrt(T)*sn;
            if bas_est; b>=b_lb; end
            if c1_est; c1>=0; end
    cvx_end
    if strcmpi(cvx_status,'Infeasible');
        warning('Problem is infeasible, perhaps noise is underestimated. Projecting original data to the cone induced by the constraints. \n');
        cvx_begin quiet
            variable c2(T)
            if bas_est; variable b; end
            if c1_est; variable c1; end
            minimize(norm(y-c2-b-c1*gd.^((0:T-1)')))
            subject to
                G*c2>=0;
                if bas_est; b>=b_lb; end
                if c1_est; c1>=0; end
        cvx_end
        sn = cvx_optval/sqrt(T);
    end
    c = c2;
end