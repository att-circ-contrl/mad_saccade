function [ model, model_amp, model_pvel, model_vel ] = sacc_model_type1( Stype, t, p )
% generate a saccade with different function (Stype)
% Stype can be 'ms', 'sigmoid', 'gumbel', 'gaussian'
% t and p are the corresponding time (msec) and parameters

% If you find the paper and/or code helpful, please cite:
% W.-W. Dai, I. W. Selesnick, J.-R. Rizzo, J. C. Rucker and T. E. Hudson.
% A parametric model for saccadic eye movement.
% IEEE Signal Processing in Medicine and Biology Symposium (SPMB), December 2016.

% Weiwei Dai
% Last edit: 01/07/17

if nargin < 1
    Stype = 'ms';
    t = -100:0.1:100;
    p = [0.6, 8, 20, -10, 0];
end

switch Stype
    case 'ms'
        % our proposed model
        % eta, c, tau, alpha, beta
        [model, model_amp, model_pvel, model_vel] = fun_sacc_model(t, p);
    case 'sigmoid'
        % a, b, alpha, beta
        fmodel = @(t, p) p(1)*(1./(1+exp((-t-p(3))/p(2)))) + p(4);
        fvel = @(t, p) 1000*((p(1)*exp(-(t-p(3))/p(2)))./(p(2)*(1+exp(-(t-p(3))/p(2))).^2));
        model = fmodel(t, p);
        model_amp = p(1);
        model_pvel = abs(1000*p(1)/4/p(2));
        model_vel = fvel(t, p);
    case 'gumbel'
        % a, b, alpha, beta
        % b > 0
        % p = [1; 4; 0; 0];
        fmodel = @(t, p) p(1)*exp(-exp((-t+p(3))/p(2))) + p(4);
        fvel = @(t, p) 1000*((p(1)/p(2))*exp((-t+p(3))/p(2)).*exp(-exp((-t+p(3))/p(2))));
        model = fmodel(t, p);
        model_amp = p(1);
        model_pvel = 1000*p(1)*exp(-1)/p(2);
        model_vel = fvel(t, p);
    case 'gaussian'
        % a, b, alpha, beta
        fmodel = @(t, p) p(1)/2*(1+erf((t-p(3))/(p(2)*sqrt(2)))) + p(4);
        fvel = @(t, p) 1000*p(1)/(p(2)*sqrt(2*pi))*exp(-(t-p(3)).^2/(2*p(2)^2));
        model = fmodel(t, p);
        model_amp = p(1);
        model_pvel = 1000*p(1)/(p(2)*sqrt(2*pi));
        model_vel = fvel(t, p);
    case 'fisk'
        % a, b, c, alpha, beta
        % t > alpha, b > 0 and c > 0
        % p = [1; 5; 6; 0; 0];
        fmodel = @(t, p) p(1)*(1 + ((t-p(4))/p(2)).^(-p(3))).^(-1) + p(5);
        fvel = @(t, p) 1000*p(1)*((p(3)/p(2))*(((t-p(4))/p(2)).^(p(3)-1)).*(1+((t-p(4))/p(2)).^p(3)).^(-2));
        model = fmodel(t, p);
        
        tp = p(4)+p(2)*((p(3)-1)/(p(3)+1))^(1/p(3));
        model_amp = p(1);
        model_pvel = fvel(tp, p);
        model_vel = fvel(t, p);
    case 'burr'
        error('to be changed!')
        % a, b, c, alpha, beta
        % p = [10; 6; 1; 0; 0];
        fmodel = @(t, p) 1 - (1+((t-p(4))/p(1)).^p(2)).^(-p(3)) + p(5);
        fvel = @(t, p) 1000*((p(2)*p(3)/p(1))*(((t-p(4))/p(1)).^(p(2)-1))*((1+((t-p(4))/p(1)).^p(2)).^(-p(3)-1)));
        model = fmodel(t, p);
        
        tp = p(1)*((p(2)-1)/(p(2)*p(3)+1))^(1/p(2));
        model_amp = 1;
        model_pvel = 1000*fvel(tp, p);
        model_vel = fvel(t, p);
end


end

