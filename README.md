# constrained-foopsi
Implementation of the constrained deconvolution spike inference algorithm in Matlab.

File constained_foopsi.m
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
%   resparse: number of times that the solution is resparsened (default 0). Currently available only with methods 'cvx', 'spgl'
%   

The noise is estimated with a power spectral density approach and the time constants from the signal autocovariance. 
