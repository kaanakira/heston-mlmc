%function [sums, cost] = nv_all(l,N)
%
% function for level l estimator using Ninomiya-Victoir method
%
% inputs:  l = level
%          N = number of paths to be generated
%
% output: sums(1) = sum(Pf-Pc)
%         sums(2) = sum((Pf-Pc).^2)
%         sums(3) = sum((Pf-Pc).^3)
%         sums(4) = sum((Pf-Pc).^4)
%         sums(5) = sum(Pf)
%         sums(6) = sum(Pf.^2)
%         cost = cost of N samples

function [sums, cost] = nv_all(l,N)

u0  = log(20);   % initial asset value 
v0 = 2;
K   = 20;   % strike
T   = 1;     % maturity
r   = 0.05;  % risk-free interest rate
sig = 0.05;   % volatility

kappa =  0.5;
theta =  0.9;

xi = theta - (sig^2)/(4*kappa);

M  = 2;

nf = M^l;
nc = nf/M;

hf = T/nf;
hc = T/nc;

sums(1:6) = 0;

for N1 = 1:10000:N
  N2 = min(10000,N-N1+1);
  
  % initialize u and v
  uf = u0*ones(1,N2);
  uc = uf;

  vf = v0*ones(1,N2);
  vc = vf;

  if l==0
    dWf1 = sqrt(hf)*randn(1,N2);
    dWf2 = sqrt(hf)*randn(1,N2);
    
    eta = (rand(1,N2)<.5)*2 - 1;
    
    vf_ = xi + (max(0,vf)-xi)*exp(-0.5*kappa*hf);
    uf = uf + 0.5*hf*(r-xi/2) + 0.5*(1/kappa)*(max(0,vf)-xi) ...
        * (exp(-0.5*kappa*hf)-1);
    
    vf = (sqrt(vf_)+0.5*sig*dWf2).^2;
    uf = uf + sqrt(vf_.*(eta==1) + vf.*(eta==-1)).*dWf1;
    
    uf = uf + 0.5*hf*(r-xi/2) + 0.5*(1/kappa)*(vf-xi) ...
        * (exp(-0.5*kappa*hf)-1);
    vf = xi + (vf-xi)*exp(-0.5*kappa*hf);
    
  else
      for n = 1:nc
          dWc1 = zeros(1,N2);
          dWc2 = zeros(1,N2);
          
          dWf1 = sqrt(hf)*randn(M,N2);
          dWf2 = sqrt(hf)*randn(M,N2);

          eta = (rand(M,N2)<.5)*2 - 1;
          
          for m = 1:M
              dWc1 = dWc1 + dWf1(m,:);
              dWc2 = dWc2 + dWf2(m,:);
              
              % Ninomiya-Victoir
              
              vf_ = xi + (max(0,vf)-xi)*exp(-0.5*kappa*hf);
              uf = uf + 0.5*hf*(r-xi/2) + 0.5*(1/kappa)*(max(0,vf)-xi) ...
                  * (exp(-0.5*kappa*hf)-1);

              vf = (sqrt(vf_)+0.5*sig*dWf2(m,:)).^2;
              uf = uf + sqrt(vf_.*(eta(m,:)==1) + vf.*(eta(m,:)==-1)).*dWf1(m,:);
              
              uf = uf + 0.5*hf*(r-xi/2) + 0.5*(1/kappa)*(vf-xi) ...
                  * (exp(-0.5*kappa*hf)-1);
              vf = xi + (vf-xi)*exp(-0.5*kappa*hf);
          end

          eta = eta(1,:); % use the first eta for coarse level
         
          vc_ = xi + (max(0,vc)-xi)*exp(-0.5*kappa*hc);
          uc = uc + 0.5*hc*(r-xi/2) + 0.5*(1/kappa)*(max(0,vc)-xi) ...
              * (exp(-0.5*kappa*hc)-1);
          
          vc = (sqrt(vc_)+0.5*sig*dWc2).^2;
          uc = uc + sqrt(vc_.*(eta==1) + vc.*(eta==-1)).*dWc1;
          
          uc = uc + 0.5*hc*(r-xi/2) + 0.5*(1/kappa)*(vc-xi) ...
              * (exp(-0.5*kappa*hc)-1);
          vc = xi + (vc-xi)*exp(-0.5*kappa*hc);

      end
  end

  Pf  = exp(-r*T)*max(0,exp(uf)-K);
  if l>0
    Pc  = exp(-r*T)*max(0,exp(uc)-K);
  else
    Pc = zeros(1,N2); 
  end
  
  sums(1) = sums(1) + sum(Pf-Pc);
  sums(2) = sums(2) + sum((Pf-Pc).^2);
  sums(3) = sums(3) + sum((Pf-Pc).^3);
  sums(4) = sums(4) + sum((Pf-Pc).^4);
  sums(5) = sums(5) + sum(Pf);
  sums(6) = sums(6) + sum(Pf.^2);
end
cost = 3*N*nf;
end
