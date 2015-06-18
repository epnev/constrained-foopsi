function [c,b,c1] = cvx_mixed_foopsi(y,b,c1,sn,b_lb,g,w,keep)

% implementation of constrained foopsi in CVX with multiple time constants
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
    if iscell(g) % multiple time constants
        nt = length(g);
        T = length(y);
        G = cell(1,nt);
        gd_vec = zeros(T,nt);
        for i = 1:nt
            G{i} = spdiags(ones(T,1)*[-g{i}(end:-1:1)',1],-length(g{i}):0,T,T);
            gd = max(roots([1,-g{i}(:)']));
            gd_vec(:,i) = gd.^((0:T-1)');
        end
        G_mat = cell2mat(G);
        cvx_begin quiet
            variable c2(T,nt)
            if bas_est; variable b; end
            if c1_est; variable c1(nt); end
            C_mat = sparse(1:T*nt,kron((1:nt)',ones(T,1)),c2(:),T*nt,nt);
            minimize(sum(norms(G_mat*C_mat,2,2)));
            subject to
                G_mat*C_mat >= 0;
                norm(y(keep) - sum(c2(keep,:),2) - b - gd_vec(keep,:)*c1)<=sqrt(sum(keep))*sn;
                if bas_est; b>=b_lb; end
                if c1_est; c1>=0; end
        cvx_end
        if strcmpi(cvx_status,'Infeasible');
            disp('Problem is infeasible, adjusting noise value.');
            cvx_begin quiet
                variable c2(T,nt)
                C_mat = sparse(1:T*nt,kron((1:nt)',ones(T,1)),c2(:),T*nt,nt);
                if bas_est; variable b; end
                if c1_est; variable c1(nt); end
                minimize(norm(y(keep)-sum(c2(keep,:),2)-b-gd_vec(keep,:)*c1))
                subject to
                    G_mat*C_mat>=0;
                    if bas_est; b>=b_lb; end
                    if c1_est; c1>=0; end
            cvx_end
            sn = cvx_optval/sqrt(sum(keep));
        end
        c = c2;
    else
        gd = max(roots([1,-g(:)']));
        T = length(y);
        G = spdiags(ones(T,1)*[-g(end:-1:1)',1],-length(g):0,T,T);
        gd_vec = gd.^((0:T-1)');
        cvx_begin quiet
            variable c2(T)
            if bas_est; variable b; end
            if c1_est; variable c1; end
            minimize(w'*(G*c2))
            subject to
                G*c2>=0;
                norm(y(keep)-c2(keep)-b-c1*gd_vec(keep))<=sqrt(sum(keep))*sn;
                if bas_est; b>=b_lb; end
                if c1_est; c1>=0; end
        cvx_end
        if strcmpi(cvx_status,'Infeasible');
            disp('Problem is infeasible, adjusting noise value.');
            cvx_begin quiet
                variable c2(T)
                if bas_est; variable b; end
                if c1_est; variable c1; end
                minimize(norm(y(keep)-c2(keep)-b-c1*gd_vec(keep)))
                subject to
                    G*c2>=0;
                    if bas_est; b>=b_lb; end
                    if c1_est; c1>=0; end
            cvx_end
            sn = cvx_optval/sqrt(sum(keep));
        end
        c = c2;
    end
end