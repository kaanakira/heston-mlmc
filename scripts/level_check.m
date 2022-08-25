%function last_level = level_check(mlmcl,eps,N,alpha,Lmax,varargin)
%
% Function to estimate the final level given RMSE bound epsilon.
% We perform this test until the condition is satisfied twice in a row for robustness
%
% inputs:  mlmc_l = function for level l estimator 
%          eps   = desired accuracy (rms error)      
%          N = number of samples to be used for the test
%          alpha = weak convergence rate of the numerical solver
%          Lmax  = maximum level of refinement  >= 2
%          varargin = optional additional user variables to be passed to mlmc_l
%
% output:  last_level = optimal refinement level L^*


function last_level = level_check(mlmcl,eps,N,alpha,Lmax,varargin)

checks = 1; % return if the result holds twice
l = 1;
theta = 0.25;

% L = 2 case is different
for i=1:2
    [sums, ~] = mlmcl(1, N, varargin{:});
    p1 = abs(sums(1)/N);
    [sums, ~] = mlmcl(2, N, varargin{:});
    p2 = abs(sums(1)/N);
    
    ml = [p1 p2];
    
    range = 0:1;
    rem = max(fliplr(ml) ./ 2.^(range*alpha)) / (2^alpha - 1);
    
    if rem<sqrt(theta)*eps
        if checks == 2
            last_level = 2;
            return;
        else
            checks = checks+1;
        end
    else
        break;
    end
end

% L > 2 case

while 1
    checks = 1;
    for k=1:2
        [sums, ~] = mlmcl(l, N, varargin{:});
        p1 = abs(sums(1)/N);
        [sums, ~] = mlmcl(l+1, N, varargin{:});
        p2 = abs(sums(1)/N);
        [sums, ~] = mlmcl(l+2, N, varargin{:});
        p3 = abs(sums(1)/N);
        
        ml = [p1 p2 p3];
        
        range = 0:2;
        rem = max(fliplr(ml) ./ 2.^(range*alpha)) / (2^alpha - 1);
        
        if rem<sqrt(theta)*eps
            if checks == 2
                last_level = l+2;
                return;
            else
                checks = checks+1;
            end
        else
            l = l+1;
            break;
        end
    end
    if Lmax == l+1 % break if l+2 exceeded Lmax (Lmax+1==l+2)
        last_level = Lmax;
        fprintf(1,'*** failed to achieve weak convergence *** \n');
        break;
    end
end
end

