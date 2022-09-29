% function [P, Cl] = nv_mlmc(Nl)
%
% multi-level Monte Carlo estimation (no-area Milstein + Ninomiya-Victoir)
% using L = length(Nl)-1 and corresponding number of samples Nl(l+1) for each level l = 0,1,...,L.
%
% P     = value
% Cl    = complexity
%
% Nl    = mumber of samples to be used at each level
% 

function [P, Cl] = nv_mlmc(Nl)

%
% Compute the finest level
%
L = length(Nl)-1;

%
% initialisation
%
suml(1:L+1) = 0;
costl(1:L+1)    = 0;

for l=0:L
    if Nl(l+1) > 0
        if l<L
            [sums, cost] = na_milstein_l(l,Nl(l+1), 2);
        else
            [sums, cost] = nv_La(l,Nl(l+1));
        end
        suml(l+1) = suml(l+1) + sums(1);
        costl(l+1)  = costl(l+1)  + cost;
    end
end
Cl = costl./Nl;

%
% evaluate multilevel estimator
%
P = sum(suml(1,:)./Nl);

end


