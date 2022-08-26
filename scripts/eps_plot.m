close all; clear all;

Eps = [ 0.0003 0.0005 0.00075 0.001 0.005 0.01 ]; % Desired accuracies

N_test = 100000; % Number of paths for the level check
N0 = 3000; % Number of initial samples on each level for MLMC

Lmin=5;
Lmax = 15;

alpha0=2;
beta0=0;
gamma0=0;

n=5; % Number of samples to use for the average of error

% Compute the theoretical value of the call option
y = call_heston_cf(20, 2, 0.9,0.5, 0.05, 0.05, 0, 1, 20);

real_err(1:length(Eps)) = 0;
Nl_(1:length(Eps),1:Lmax+1) = 0;

for i=1:length(Eps)
    eps = Eps(i);
    real_av(1:n) = 0; % array to store result for each simulation
    L=0;
    
    for j=1:n
        
        if L==0
            L = level_check(@new_all_a,eps,N_test,2,Lmax,'Call',2);
        end
        
        % MLMC estimation
        [P, Nl, C] = mlmc_a(N0,eps,L,@new_l,@new_La,2);
        
        % Compute emprirical error
        real_av(j) = abs(P-y);
    end    
    Nl_(i,:)=[Nl zeros(1, Lmax+1 - length(Nl))]';
    real_err(i) = sum(real_av)/n;    % take the average
end


set(0,'DefaultAxesColorOrder',[0 0 0]);
loglog([0.0001 Eps 0.02], [0.0001 Eps 0.02], ':',  Eps, real_err, 'o-')
xlabel('accuracy $\epsilon$','Interpreter','latex'); 
ylabel('$\bar{E}_{new}$','Interpreter','latex');
current = axis; 

savefig('error.fig')
  
  
% Below functions are from the paper "An Analysis of the Heston Stochastic Volatility Model: Implementation and Calibration using Matlab" by Ricardo Crisostomo

function y = call_heston_cf(s0, v0, vbar, a, vvol, r, rho, t, k)
% Heston call value using characteristic functions.
% y = call_heston_cf(s0, v0, vbar, a, vvol, r, rho, t, k)
% Inputs:
% s0: stock price
% v0: initial volatility (v0^2 initial variance)
% vbar: long-term variance mean
% a: variance mean-reversion speed
% vvol: volatility of the variance process
% r: risk-free rate
% rho: correlation between the Weiner processes of the stock price and its variance
% t: time to maturity
% k: option strike
% chfun_heston: Heston characteristic function

% 1st step: calculate pi1 and pi2
 % Inner integral 1
int1 = @(w, s0, v0, vbar, a, vvol, r, rho, t, k) real(exp(-1i.*w*log(k)).*chfun_heston(s0,...
v0, vbar, a, vvol, r, rho, t, w-1i)./(1i*w.*chfun_heston(s0, v0, vbar, a, vvol, r, rho, t,...
-1i))); % inner integral1
int1 = integral(@(w)int1(w,s0, v0, vbar, a, vvol, r, rho, t, k),0,100); % numerical integration
pi1 = int1/pi+0.5; % final pi1
 % Inner integral 2:
int2 = @(w, s0, v0, vbar, a, vvol, r, rho, t, k) real(exp(-1i.*w*log(k)).*chfun_heston(s0, ...
v0, vbar, a, vvol, r, rho, t, w)./(1i*w));
int2 = integral(@(w)int2(w,s0, v0, vbar, a, vvol, r, rho, t, k),0,100);int2 = real(int2);
pi2 = int2/pi+0.5; % final pi2
% 2nd step: calculate call value
y = s0*pi1-exp(-r*t)*k*pi2;

function y = chfun_heston(s0, v0, vbar, a, vvol, r, rho, t, w)
% Heston characteristic function.
% Inputs:
% s0: stock price
% v0: initial volatility (v0^2 initial variance)
% vbar: long-term variance mean
% a: variance mean-reversion speed
% vvol: volatility of the variance process
% r : risk-free rate
% rho: correlation between the Weiner processes for the stock price and its variance
% w: points at which to evaluate the function
% Output:
% Characteristic function of log(St) in the Heston model

% Interim calculations
alpha = -w.*w/2 - 1i*w/2;
beta = a - rho*vvol*1i*w;
gamma = vvol*vvol/2;
h = sqrt(beta.*beta - 4*alpha*gamma);
rplus = (beta + h)/vvol/vvol;
rminus = (beta - h)/vvol/vvol;
g=rminus./rplus;
% Required inputs for the characteristic function
C = a * (rminus * t - (2 / vvol^2) .* log((1 - g .* exp(-h*t))./(1-g)));
D = rminus .* (1 - exp(-h * t))./(1 - g .* exp(-h*t));
% Characteristic function evaluated at points w
y = exp(C*vbar + D*v0 + 1i*w*log(s0*exp(r*t))); 
end
end
