function [ model, model_amp, model_pvel, model_vel ] = fun_sacc_model( t, p ) 
% generate saccade model derived from
% main sequence formula I: P.V. = eta * (1 - exp(-A/c))
% p = [eta, c, tau, alpha, beta] are paramters

% If you find the paper and/or code helpful, please cite:
% W.-W. Dai, I. W. Selesnick, J.-R. Rizzo, J. C. Rucker and T. E. Hudson.
% A parametric model for saccadic eye movement.
% IEEE Signal Processing in Medicine and Biology Symposium (SPMB), December 2016.

% eta = 0.6, c = 8, tau = 20, alpha = -10, beta = 0

% Weiwei Dai
% Last edit: 01/07/17

eta = p(1);
c = p(2);
tau = p(3);
alpha = p(4);
beta = p(5);

f = @(t) c/4*exp(-2*eta*t.*sign(t)/c) + eta*t.*(t>=0);
f_vel = @(t) 1000*(eta/2.*exp(2*eta*t/c).*(t<=0) + (eta - eta/2.*exp(-2*eta*t/c)).*(t>0));
fpvel = @(x) 1000*eta*(1-exp(-x/c)); 

model = f(t - alpha) - f(t - alpha - tau) + beta;
model_amp = abs(eta*tau);
model_pvel = fpvel(model_amp);

model_vel = f_vel(t - alpha) - f_vel(t - alpha - tau);

end

