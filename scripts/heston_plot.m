% This code is adapted from https://people.maths.ox.ac.uk/~gilesm/mlmc/matlab/mlmc_test.m

close all; clear all;

Lmin  = 3;      % minimum refinement level
Lmax  = 15;     % maximum refinement level

L = 7; % number of levels to observe for the convergences of mean and variance

Eps   = [ 0.0001 0.0003 0.0005 0.00075 0.001 0.005]; % desired accuracies 

nvert = 2;

alpha0=2;
beta0=0;
gamma0=0;

% For Variance and mean

del1 = [];
del2 = [];
var1 = [];
var2 = [];
cost = [];

var1_m = [];
del1_m = [];
del2_m = [];

N =2000000; % Number of paths for each level for the convergence test

for l = 0:L
  sums = 0;
  cst  = 0;
  sums_m = 0;
  
  parfor j=1:100
    RandStream.setGlobalStream( ...
    RandStream.create('mrg32k3a','NumStreams',100,'StreamIndices',j));

    if l<L
        [sums_j_m, cst_j_m] = na_milstein_l(l, N/100, 2);
        [sums_j, cst_j] = new_l(l, N/100, 2);
    else
        [sums_j_m, cst_j_m] = na_milstein_l(l, N/100, 2);
        [sums_j, cst_j] = new_La(l, N/100);
    end
    
    sums = sums + sums_j/N;
    cst  = cst  + cst_j/N;
    
    sums_m = sums_m + sums_j_m/N;
    
  end  


  cost = [cost cst(1)];
  del1 = [del1 sums(1)];
  del2 = [del2 sums(5)];
  var1 = [var1 sums(2)-sums(1)^2 ];
  var2 = [var2 sums(6)-sums(5)^2 ];
  var2 = max(var2, 1e-10);
  
  del1_m = [del1_m sums_m(1)];
  del2_m = [del2_m sums_m(5)];
  var1_m = [var1_m sums_m(2)-sums_m(1)^2 ];
  
end

% regression for alpha, beta
L1 = 2;
L2 = L+1;
pa    = polyfit(L1:L2-1,log2(abs(del1(L1:L2-1))),1);  alpha = -pa(1);
pa_m    = polyfit(L1:L2,log2(abs(del1_m(L1:L2))),1);  alpha_m = -pa_m(1);
pb    = polyfit(L1:L2-1,log2(abs(var1(L1:L2-1))),1);  beta  = -pb(1);
pb_m  = polyfit(L1:L2,log2(abs(var1_m(L1:L2))),1);  beta_m  = -pb_m(1);

fprintf('\n alpha = %f  (exponent for MLMC weak convergence)\n',alpha);
fprintf('alpha_m = %f  (exponent for MLMC weak convergence)\n',alpha_m);
fprintf('beta  = %f  (exponent for MLMC variance) \n',beta);
fprintf('beta_m  = %f  (exponent for MLMC variance) \n',beta_m);

% reset random number generators
reset(RandStream.getGlobalStream);
spmd
  RandStream.setGlobalStream( ...
  RandStream.create('mrg32k3a','NumStreams',numlabs,'StreamIndices',labindex));
end


figs(1) = figure; 
pos=get(gcf,'pos'); pos(3:4)=pos(3:4).*[1.0 0.75*nvert]; set(gcf,'pos',pos);

% Figure 1 Variance

% reference line for slope -2
refline = log2(var1(2))+5;
for i=1:L
    refline = [refline -2*i+log2(var1(2))+5];
end

% reference line for slope -4
refline_ = log2(var1(2))+5;
for i=1:L
    refline_ = [refline_ -4*i+log2(var1(2))+5];
end

set(0,'DefaultAxesColorOrder',[0 0 0]);

set(0,'DefaultAxesLineStyleOrder',':o|-^|-*');

subplot(nvert,2,1)
plot(0:L,log2(var2),1:L,log2(var1_m(2:end)),1:L,log2(var1(2:end)),0:L,refline, ':', 0:L,refline_,':')
xlabel('level $\ell$','Interpreter','latex'); 
ylabel('$\log_2$ variance','Interpreter','latex');
current = axis; axis([ 0 L current(3:4) ]);
legend({'$P^{\mbox{new}}_\ell$', '$P_\ell\!-\!P_{\ell-1}$',...
    '$P^{\mbox{new}}_\ell\!-\! P_{\ell-1}^{\mbox{new}}$'}, ...
        'Interpreter','latex','Location','SouthWest')

% Figure 2 Mean

set(0,'DefaultAxesLineStyleOrder',':o|-^|-*');

% reference line for slope -2
refline = log2(abs(del1(2)));
for i=1:L
    refline = [refline -2*i+log2(abs(del1(2)))];
end

% reference line for slope -1
refline2 = log2(abs(del1(2)));
for i=1:L
    refline2 = [refline2 -i+log2(abs(del1(2)))];
end

subplot(nvert,2,2)
plot(0:L,log2(abs(del2)),1:L,log2(abs(del1_m(2:end))), 1:L,log2(abs ...
    (del1(2:end))), 0:L, refline, ':', 0:L, refline2, ':')
xlabel('level $\ell$','Interpreter','latex'); 
ylabel('$\log_2 |\mbox{mean}|$','Interpreter','latex');
current = axis; axis([ 0 L current(3:4) ]);
legend({'$P_\ell^{\mbox{new}}$','$P_\ell\!-\! P_{\ell-1}$', ...
    '$P_\ell^{\mbox{new}}\!-\! P_{\ell-1}^{\mbox{new}}$'}, ...
        'Interpreter','latex','Location','SouthWest')


% For cost

std_cost = zeros(1,length(Eps));
new_cost = zeros(1,length(Eps));
mc_cost = zeros(1,length(Eps));
nv_cost = zeros(1,length(Eps));

Nl_nv = zeros(Lmax+1,length(Eps));
Nl_new = zeros(Lmax+1,length(Eps));
Nl_a = zeros(Lmax+1,length(Eps));

theta = 0.25;

N_test = 500000; % number of samples to determine L
N0 = 3000;    % initial number of samples on coarsest levels

for i=1:length(Eps)
    eps = Eps(i);
    
    % Determine the finest levels (L's)
    L_ = [level_check(@na_milstein_l,eps,N_test,1,Lmax,2)  ...
                       level_check(@nv_all,eps,N_test,2,Lmax) ...
                       level_check(@new_all_a,eps,N_test,2,Lmax,'Call',2)];
    
    % MLMC simulations
    [P , Nl , Cl ] = mlmc_a(N0,eps,L_(1),@na_milstein_l,@na_milstein_l, 2);
    [P_nv , Nl_nv_ , Cl_nv ] = mlmc_a(N0,eps,L_(2),@na_milstein_l,@nv_La, 2);
    [P_new , Nl_new_ , Cl_new ] = mlmc_a(N0,eps,L_(3),@new_l,@new_La, 2);

    std_cost(i) = sum(Nl.*Cl);
    new_cost(i) = sum(Nl_new_.*Cl_new);
    nv_cost(i) = sum(Nl_nv_.*Cl_nv);

    Nl_a(:,i) = [Nl zeros(1, Lmax+1 - length(Nl))]';
    Nl_new(:,i) = [Nl_new_  zeros(1, Lmax+1 - length(Nl_new_))]';
    Nl_nv(:,i) = [Nl_nv_  zeros(1, Lmax+1 - length(Nl_nv_))]';
end

l = repmat(0:Lmax+1,length(Eps),1);

% Figure 3

set(0,'DefaultAxesLineStyleOrder',':o|-x|:d|-*|:s|-^');
[l_new, Nl_new] = trimmer(l, Nl_new);

subplot(nvert,2,3)
semilogy(l_new, Nl_new)
xlabel('level $\ell$','Interpreter','latex'); 
ylabel('$N_\ell$ (New)','Interpreter','latex'); 
current = axis; axis([ 0 size(Nl_new,1)-1 current(3:4) ]);
for i=1:length(Eps)
  labels{i} = num2str(Eps(i));
end
legend(labels,'Location','NorthEast')

% Figure 4

set(0, 'DefaultAxesLinestyleOrder', '-o|--^|-*');

subplot(nvert,2,4)
loglog(Eps,Eps.^2.*std_cost(:)', Eps,Eps.^2.*nv_cost(:)', Eps,Eps.^2.*new_cost(:)')
xlabel('accuracy $\varepsilon$','Interpreter','latex'); 
ylabel('$\varepsilon^2$ Cost','Interpreter','latex');
current = axis; axis([ Eps(1) Eps(end) current(3:4) ]);
legend('Antithetic', 'NV', 'New')

savefig('heston.fig')

function [l_new, N_new] = trimmer(l, N)
idx = 0;
for i=0:length(N)-1
    b = N(end-i, :);
    if any(b(:) ~=	0)
        break
    end
    idx = idx+1;
end
l_ = l';
N_new = N(1:end-idx,:);
l_new = l_(1:end-idx-1,:);
end
