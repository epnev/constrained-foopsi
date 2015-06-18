function [tau_new, c_new, accept] = propose_new_time_constant(y,sp,tau_old,tau_min,tau_max,tau_std,sn,c1)

% y:        data signal 
% sp:           spiking signal
% tau_old:      current time constant
% tau_min:      lower bound
% tau_max:      upper bound
% tau_std:      width of proposal distribution
% sn:           noise standard deviation
% c1:           initial value

% tau_new:      new time constant
% c_new:        updated fluorescence signal
% accept:       flag for move acceptance

ITER = 100;
T = length(y);
tau_old = sort(tau_old);
tau_new = tau_old;
p = length(tau_old);
gr = exp(-1./tau_old);
if p == 1
    gr = [0,gr];
    tau_old = [0,tau_old];
end
gd_vec = max(gr).^((0:T-1)');
accept = zeros(p,1);
GSP2 = GSP(gr(2),sp);
GSP1 = GSP(gr(1),sp);
c = (GSP2 - GSP1)/diff(gr);
logC = -norm(y - c1*gd_vec - c)^2;
gr_ = gr;

for iter = 1:ITER
    if p > 1    % update rise time constant
        tau_temp = tau_old(1) + tau_std(1)*randn;
        while tau_temp > tau_old(2) || tau_temp < tau_min
            tau_temp = tau_old(1) + tau_std(1)*randn;
        end
        gr_(1) = exp(-1/tau_temp);
        GSP1_ = GSP(gr_(1),sp);
        c_ = (GSP2 - GSP1_)/diff(gr_);
        logC_ = -norm(y(:)-c_(:)-c1*gd_vec)^2;
        prior_ratio = 1;
        % prior_ratio = gampdf(tau_(2),12,1)/gampdf(tau(2),12,1);
        ratio = exp((logC_-logC)/(2*sn^2))*prior_ratio;
        if rand < ratio %accept
            tau_new(1) = tau_temp;
            GSP1 = GSP1_;
            c = c_;
            accept(1) = accept(1) + 1
            logC = logC_;
            gr = gr_;
        end    
    end

    tau_temp = tau_old(2)+(tau_std(2)*randn);
    while tau_temp > tau_max || tau_temp<tau_new(1)
        tau_temp = tau_old(2)+(tau_std(2)*randn);
    end  
    gr_(2) = exp(-1/tau_temp);
    GSP2_ = GSP(gr_(2),sp);
    c_ = (GSP2_ - GSP1)/diff(gr_);
    gd_vec_ = gr_(2).^((0:T-1)');
    logC_ = -norm(y(:)-c_(:) - c1*gd_vec_)^2;
    prior_ratio = 1;
    ratio = exp((logC_-logC)/(2*sn^2))*prior_ratio;
    if rand < ratio
        tau_new(p) = tau_temp;
        GSP2 = GSP2_;
        c = c_;
        accept(p) = accept(p) + 1
        logC = logC_;
        gr = gr_;
        gd_vec = gd_vec_;
    end
end
    c_new = c + c1*gd_vec;

    function c_ = GSP(gt,sp)
        if gt == 0
            c_ = zeros(length(sp),1);
        else
            c_ = filter(1,[1,-gt],sp);
        end
    end

end