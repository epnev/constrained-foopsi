function [c,b,c1,g,sn,sp] = constrained_foopsi(y,b,c1,g,sn,options)
% spike inference using a constrained foopsi approach:
%      min      sum(sp)
%    c,sp,b,c1
%      subject to: sp >= 0
%                   b >= 0
%                  G*c = sp
%                   c1 >= 0
%           ||y-b-c - c_in|| <= sn*sqrt(T)

%   Variables:
%   y:      raw fluorescence data (vector of length(T))
%   c:      denoised calcium concentration (Tx1 vector)
%   b:      baseline concentration (scalar)
%  c1:      initial concentration (scalar)
%   g:      discrete time constant(s) (scalar or 2x1 vector)
%  sn:      noise standard deviation (scalar)
%  sp:      spike vector (Tx1 vector)

%   USAGE:
%   [c,b,c1,g,sn,sp] = constrained_foopsi(y,b,c1,g,sn,OPTIONS)
%   The parameters b,cin,g,sn can be given or else are estimated from the data

%   OPTIONS: (stuct for specifying options)
%         p: order for AR model, used when g is not given (default 2)
%    method: methods for performing spike inference
%   available methods: 'dual' uses dual ascent
%                       'cvx' uses the cvx package available from cvxr.com (default)
%                      'lars' uses the least regression algorithm 
%                     'spgl1' uses the spgl1 package available from
%                     math.ucdavis.edu/~mpf/spgl1/  (usually fastest)
%   bas_nonneg:   flag for setting the baseline lower bound. if 1, then b >= 0 else b >= min(y)
%   noise_range:  frequency range over which the noise power is estimated. Default [Fs/4,Fs/2]
%   noise_method: method to average the PSD in order to obtain a robust noise level estimate
%   lags:         number of extra autocovariance lags to be considered when estimating the time constants
%   resparse:     number of times that the solution is resparsened (default 0). Currently available only with methods 'cvx', 'spgl'

% Written by Eftychios Pnevmatikakis 

defoptions.p = 2;
defoptions.method = 'cvx';
defoptions.bas_nonneg = 1;              % nonnegativity option for baseline estimation
defoptions.noise_range = [0.25,0.5];    % frequency range over which to estimate the noise
defoptions.noise_method = 'logmexp';    % method for which to estimate the noise level
defoptions.lags = 5;                    % number of extra lags when computing the AR coefficients
defoptions.resparse = 0;                % number of times to re-sparse solution

if nargin < 6
    options = defoptions;
    if nargin < 5
        sn = [];
        if nargin < 4
            g = [];
            if nargin < 3
                c1 = [];
                if nargin < 2
                    b = [];
                end
            end
        end
    end
end
      
if ~isfield(options,'p');  options.p = defoptions.p;  end
if ~isfield(options,'method'); options.method = defoptions.method; end
if ~isfield(options,'bas_nonneg'); options.bas_nonneg = defoptions.bas_nonneg; end
if ~isfield(options,'noise_range'); options.noise_range = defoptions.noise_range; end
if ~isfield(options,'noise_method'); options.noise_method = defoptions.noise_method; end
if ~isfield(options,'lags'); options.lags = defoptions.lags; end
if ~isfield(options,'resparse'); options.resparse = defoptions.resparse; end

method = options.method;    
if isempty(b);
    bas_est = 1;
else
    bas_est = 0;
end
if isempty(c1)
    c1_est = 1;
else
    c1_est = 0;
end
if isempty(sn)
    sn = GetSn(y,options.noise_range,options.noise_method);
end
if isempty(g)
    g = estimate_time_constants(y,options.p,sn,options.lags);
    while max(abs(roots([1,-g(:)']))>1)
        warning('No stable AR(%i) model found. Checking for AR(%i) model \n',options.p,options.p+1);
        options.p = options.p + 1;
        g = estimate_time_constants(y,options.p,sn,options.lags);
    end
    %fprintf('Stable AR(%i) model found \n',options.p);
end
if options.bas_nonneg  % lower bound for baseline
    b_lb = 0;
else
    b_lb = min(y);
end

if strcmpi(method,'dual'); method = 'dual';
elseif strcmpi(method,'cvx'); method = 'cvx';
elseif strcmpi(method,'lars'); method = 'lars';
elseif strcmpi(method,'spgl1'); method = 'spgl1';
else fprintf('Invalid choice of method. Using dual ascent. \n'); method = 'dual';
end
    
pathCell = regexp(path, pathsep, 'split');
g = g(:);
y = y(:);
T = length(y);
G = spdiags(ones(T,1)*[-g(end:-1:1)',1],-length(g):0,T,T);
gd = max(roots([1,-g']));  % decay time constant for initial concentration
gd_vec = gd.^((0:T-1)');

switch method
    case 'dual'
         v = G'*ones(T,1);
        thr = sn*sqrt(T);
        if bas_est; b = 0; end
        if c1_est; c1 = 0; end
        myfun = @(Ald) lagrangian_temporal_gradient(Ald,thr^2,y-b-c1*gd_vec,bas_est,c1_est);
        c = [G\max(G*y,0);zeros(bas_est);zeros(c1_est)];
        options_dual = optimset('GradObj','On','Display','Off','Algorithm','interior-point','TolX',1e-8);
        ld_in = 10;
        [ld,~,flag] = fmincon(myfun,ld_in,[],[],[],[],0,[],[],options_dual);        
        if (flag == -2) || (flag == -3)
            warning('Problem seems unbounded or infeasible. Try a different method.');
        end
        if bas_est; b = c(T+bas_est); end
        if c1_est; c1 = c(end); end
        c = c(1:T);
        sp = G*c;
    case 'cvx'
        onPath = ~isempty(strfind(pathCell, 'cvx'));
        if onPath
            c = zeros(T,1+options.resparse);
            sp = zeros(T,1+options.resparse);
            bas = zeros(1+options.resparse,1);
            cin = zeros(1+options.resparse,1);
            w_ = ones(T,1);
            for rep = 1:options.resparse+1
                [c(:,rep),bas(rep),cin(rep)] = cvx_foopsi(y,b,c1,sn,b_lb,g,w_);
                sp(:,rep) = G*c(:,rep);                
                w_ = 1./(max(sp(:,rep),0) + 1e-8);
            end
            sp(sp<1e-6) = 0;
            c = G\sp;
            b = bas;
            c1 = cin;
        else
            error('CVX does not appear to be on the MATLAB path. It can be downloaded from cvxr.com \n');
        end
    case 'lars'
         Ginv = [full(G\speye(T)),ones(T,bas_est),gd_vec*ones(1,c1_est)];
         if bas_est; b = 0; end
         if c1_est; c1 = 0; end    
         [~, ~, spikes, ~, ~] = lars_regression_noise(y-b_lb*bas_est - b - c1*gd_vec, Ginv, 1, sn^2*T);
         sp = spikes(1:T);
         b = (spikes(T+bas_est)+b_lb)*bas_est + b*(1-bas_est);
         c1 = spikes(end)*c1_est + c1*(1-c1_est);
         c = G\sp;
    case 'spgl1'
        onPath = ~isempty(strfind(pathCell, 'spgl1'));
        if onPath
            Gx = @(x,mode) G_inv_mat(x,mode,T,g,gd_vec,bas_est,c1_est);
            c = zeros(T,1+options.resparse);
            sp = zeros(T,1+options.resparse);
            bas = zeros(1+options.resparse,1);
            cin = zeros(1+options.resparse,1);
            w_ = ones(T,1);
            for rep = 1:options.resparse+1
                if bas_est; b = 0; w_ = [w_;1e-10]; end
                if c1_est; c1 = 0; w_ = [w_;1e-10]; end
                options_spgl = spgSetParms('project',@NormL1NN_project ,'primal_norm', @NormL1NN_primal,'dual_norm',@NormL1NN_dual,'verbosity',0,'weights',w_);
                [spikes,r,~,~] = spg_bpdn( Gx, y-b_lb*bas_est - (1-bas_est)*b-(1-c1_est)*c1*gd_vec, sn*sqrt(T),options_spgl);
                c(:,rep) = Gx([spikes(1:T);0],1);                                  %% calcium signal
                bas(rep) = b*(1-bas_est) + bas_est*spikes(T+bas_est)+b_lb*bas_est;       %% baseline
                cin(rep) = c1*(1-c1_est) + c1_est*spikes(end);
                sp(:,rep) = spikes(1:T);                                           %% spiking signal
                w_ = 1./(spikes(1:T)+1e-8);
            end
            b = bas;
            c1 = cin;
            %sn = norm(r)/sqrt(T);
        else
            error('SPGL1 does not appear to be on the MATLAB path. It can be downloaded from math.ucdavis.edu/~mpf/spgl1 \n');
        end
end

    function sn = GetSn(Y,range_ff,method)
        % estimate noise level with a power spectral density method
        L=length(Y);
        if ~isempty(which('pwelch'));
            [psd_Y,ff]=pwelch(Y,round(L/8),[],1000,1);
        else
            xdft = fft(Y);
            xdft = xdft(:,1:round(N/2)+1);
            psd_Y = (1/L) * abs(xdft).^2;
            ff = 0:1/N:1/2;
            psd_Y(2:end-1) = 2*psd_Y(2:end-1);
        end
        ind=ff>range_ff(1);
        ind(ff>range_ff(2))=0;
        switch method
            case 'mean'
                sn=sqrt(mean(psd_Y(ind)/2));
            case 'median'
                sn=sqrt(median(psd_Y(ind)/2));
            case 'logmexp'
                sn = sqrt(exp(mean(log(psd_Y(ind)/2))));
        end
    end

    
    function g = estimate_time_constants(y,p,sn,lags)
        % estimate time constants from autocorrelation function
        
        lags = lags + p;
        xc = xcov(y,lags,'biased');
        xc = xc(:);
        A = toeplitz(xc(lags+(1:lags)),xc(lags+(1:p))) - sn^2*eye(lags,p);
        g = pinv(A)*xc(lags+2:end);            
    end


%     function Zin = plain_foopsi(H,D,I_est,eps)
% 
%         % solves argmin ||X-H||^2 subject to D*X>=0 with an interior point method
%         % using I_est as the initial value and eps as the initial barrier weight
% 
%         ln = length(H);
%         step_back_frac = 0.5;
%         iter = 0;
%         if nargin == 2
%             I_est = 1e-3*ones(ln,1);
%             eps = 1;
%         end
%         Zin = I_est(:);
% 
%         if nargin == 3
%             eps = 1;
%         end
%         while eps>1e-8
%             n = D*Zin;
%             nnd = 10;
%             E = norm(Zin-H)^2 - eps*sum(log(D*Zin));
%             grad = 2*(Zin-H) - eps*D'*(n.^(-1));
%             Hs = 2*speye(ln) + eps*D'*spdiags(n.^(-2),0,ln,ln)*D;          
%             while nnd/2>1
%                 iter = iter + 1;
%                 Z_dir = -Hs\grad;
%                 hit = -n./(D*Z_dir);
%                 if all(hit<0)
%                     s = 1;
%                 else
%                     s = min(1,.9*min(hit(hit>=0)));
%                 end
%                 E_new = E; s = s/step_back_frac;
%                 x_dg = grad'*Z_dir;
%                 while E_new > E + 0.25*s*x_dg
%                     s=s*step_back_frac; 
%                     Z_new = Zin + s*Z_dir;
%                     n = D*Zin;
%                     E_new = norm(Z_new-H)^2 - eps*sum(log(D*Z_new));
%                 end
%                 %E = E_new;
%                 Zin = Zin + s*Z_dir;
%                 nnd = -x_dg;
%                 E = norm(Zin-H)^2 - eps*sum(log(D*Zin));
%                 n = D*Zin;
%                 grad = 2*(Zin-H) - eps*D'*(n.^(-1));
%                 Hs = 2*speye(ln) + eps*D'*spdiags(n.^(-2),0,ln,ln)*D;
%             end
%             eps = eps/10;
%         end
%     end
    
%     function [f,grad] = lagrangian_temporal_grad(Al,thr,y_raw)
%         %options2 = optimset('Display','Off','Algorithm','interior-point-convex');
%         %c = quadprog(2*sum(Al.*(a.^2))*speye(T),v-2*((la>1)*sum(spdiags(Al.*a,0,la,la)*y,1)' + (la==1)*y'*(Al.*a)),-G,zeros(T,1),[],[],[],[],[],options2);
%         c = plain_foopsi((y_raw*(Al)-v/2)/sum(Al),G);
%         f = v'*c;    
%         grad = (sum((c-y_raw).^2)-thr);
%         f = f + Al(:)'*grad;
%         f = -f;
%         grad = -grad;
%     end

    function [f,grad] = lagrangian_temporal_gradient(Al,thr,y_raw,bas_flag,c1_flag)
        options_qp = optimset('Display','Off','Algorithm','interior-point-convex');
        H = [speye(T),ones(T,bas_flag),gd_vec*ones(1,c1_flag);...
            ones(bas_flag,T),T*ones(bas_flag),(1-gd^T)/(1-gd)*ones(c1_flag,bas_flag);...
            (gd_vec*ones(1,c1_flag))',(1-gd^T)/(1-gd)*ones(c1_flag,bas_flag),(1-gd^(2*T))/(1-gd^2)*ones(c1_flag,c1_flag)];
        Ay = [y_raw;sum(y_raw)*ones(bas_flag);gd_vec'*y_raw*ones(c1_flag)];
        c = quadprog(2*Al(1)*H,[v;zeros(bas_flag+c1_flag,1)]-2*Al(1)*Ay,[-G,sparse(T,bas_flag+c1_flag);sparse(bas_flag+c1_flag,T),-speye(bas_flag+c1_flag)]...
            ,[sparse(T,1);-b_lb*ones(bas_flag);zeros(c1_flag)],[],[],[],[],c,options_qp);
        f = v'*c(1:T);    
        grad = [sum((c(1:T)-y_raw + c(T+bas_flag)*bas_flag + c(end)*gd_vec*c1_flag).^2)-thr];
        f = f + Al(:)'*grad;
        f = -f;
        grad = -grad;
    end

    function b = G_inv_mat(x,mode,NT,gs,gd_vec,bas_flag,c1_flag)
        if mode == 1
            b = filter(1,[1;-gs(:)],x(1:NT)) + bas_flag*x(NT+bas_flag) + c1_flag*gd_vec*x(end);
           %b = G\x(1:NT) + x(NT+bas_flag)*bas_flag + x(end)*c1_flag;
        elseif mode == 2
            b = [flipud(filter(1,[1;-gs(:)],flipud(x)));ones(bas_flag,1)*sum(x);ones(c1_flag,1)*(gd_vec'*x)];
           %b = [G'\x;ones(bas_flag,1)*sum(x);ones(c1_flag,1)*(gd_vec'*x)] ;
        end
    end

end
