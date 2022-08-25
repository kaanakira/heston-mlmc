% This function is a modified version of mlmc.m from https://people.maths.ox.ac.uk/~gilesm/mlmc/
%
% function [P, Nl, Cl] = mlmc_a(N0,eps,L,mlmc_l,mlmc_L, varargin)
%
% accelerated multi-level Monte Carlo estimation 
%
% P     = value
% Nl    = number of samples at each level
% Cl    = complexity
%
% N0    = initial number of samples         > 0
% eps   = desired accuracy (rms error)      > 0 
% L     = the optimal refinement level L*   >= 2 
%
% varargin = optional additional user variables to be passed to mlmc_l
%
% mlmc_l = function for level l estimator 
% mlmc_L = function for level L estimator (acceleration)
%
% [sums, cost] = mlmc_fn(l,N, varargin)     low-level routine
%
% inputs:  l = level
%          N = number of samples
%          varargin = optional additional user variables
%
% output: sums(1) = sum(Y)
%         sums(2) = sum(Y.^2)
%         where Y are iid samples with expected value:
%         E[P_0]           on level 0
%         E[P_l - P_{l-1}] on level l>0
%         cost = cost of N samples

function [P, Nl, Cl] = mlmc_a(N0,eps,L,mlmc_l,mlmc_L, varargin)

%
% check input parameters
%
  if (L<2)
    error('error: needs L >= 2');
  end

  if (N0<=0 || eps<=0)
    error('error: needs N0>0, eps>0 \n');
  end
  
%
% initialisation
%
Nl(1:L+1)       = 0;
suml(1:2,1:L+1) = 0;
dNl(1:L+1)      = N0;
costl(1:L+1)    = 0;

while sum(dNl)>0 
%
% update sample sums
%
    for l=0:L
      if dNl(l+1) > 0
        if l<L
            [sums, cost] = mlmc_l(l,dNl(l+1), varargin{:});
        else
            [sums, cost] = mlmc_L(l,dNl(l+1), varargin{:});
        end
        Nl(l+1)     = Nl(l+1)     + dNl(l+1);
        suml(1,l+1) = suml(1,l+1) + sums(1);
        suml(2,l+1) = suml(2,l+1) + sums(2);
        costl(l+1)  = costl(l+1)  + cost;
      end
    end
    
%
% compute absolute average, variance and cost
%
    ml = abs(   suml(1,:)./Nl);
    Vl = max(0, suml(2,:)./Nl - ml.^2);
    Cl = costl./Nl;
    
    tht = 0.25;
    
%     
%   set optimal number of additional samples
%     
    Ns  = ceil( sqrt(Vl./Cl)*sum(sqrt(Vl.*Cl)) / ((1-tht)*eps^2) );
    dNl = max(0, Ns-Nl);
                
end

%
% finally, evaluate multilevel estimator
%
P = sum(suml(1,:)./Nl);

end
