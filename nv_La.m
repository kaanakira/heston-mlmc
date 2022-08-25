function [sums, cost] = nv_La(l,N,~)
%-------------------------------------------------------
%
% level L, acceleration
%
% Ninomiya-Victoir on the fine grid and no-area Milstein on coarse grid
%
u0  = log(20);   % initial asset value log1 = 0
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

if l==0
   disp('error: needs l > 0');
end

for N1 = 1:10000:N
  N2 = min(10000,N-N1+1);
  
  % initialize u and v
  uf = u0*ones(1,N2);
  uc = uf;
  uf_a = uf;

  vf = v0*ones(1,N2);
  vc = vf;
  
  vf = [vf vf vf vf];
  uf = [uf uf uf uf];
  
  for n = 1:nc
      dWc1 = zeros(1,N2);
      dWc2 = zeros(1,N2);
      
      dWf1 = sqrt(hf)*randn(M,N2);
      dWf2 = sqrt(hf)*randn(M,N2);
      
      nu = (rand(M,N2)<.5)*2 - 1;
      
      % antithetic treatment
      dWf1 = [dWf1 dWf1 dWf1(M:-1:1,:) dWf1(M:-1:1,:)];
      dWf2 = [dWf2 dWf2 dWf2(M:-1:1,:) dWf2(M:-1:1,:)];
      nu = [nu -nu nu -nu];
      
      for m = 1:M
          dWc1 = dWc1 + dWf1(m,1:N2);
          dWc2 = dWc2 + dWf2(m,1:N2);
          
          vf_ = xi + (max(0,vf)-xi)*exp(-0.5*kappa*hf);
          uf = uf + 0.5*hf*(r-xi/2) + 0.5*(1/kappa)*(max(0,vf)-xi) ...
              * (exp(-0.5*kappa*hf)-1);
          
          vf = (sqrt(vf_)+0.5*sig*dWf2(m,:)).^2;
          uf = uf + sqrt(vf_.*(nu(m,:)==1) + vf.*(nu(m,:)==-1)).*dWf1(m,:);
          
          uf = uf + 0.5*hf*(r-xi/2) + 0.5*(1/kappa)*(vf-xi) ...
              * (exp(-0.5*kappa*hf)-1);
          vf = xi + (vf-xi)*exp(-0.5*kappa*hf);
      end
      
      uc = uc + (r-0.5*vc)*hc + sqrt(max(0,vc)).*dWc1  ...
          + 0.25*sig*(dWc1.*dWc2);
      
      vc = vc + kappa*(theta-vc)*hc + sig*sqrt(max(0,vc)).*dWc2 ...
          + 0.25*(sig^2)*(dWc2.^2-hc);
      
  end
  
  Pf  = exp(-r*T)*max(0,exp(uf)-K);
  Pf = 0.25*(Pf(1:N2)+Pf(1+N2:2*N2)+Pf(1+2*N2:3*N2)+Pf(1+3*N2:end));
  Pc  = exp(-r*T)*max(0,exp(uc)-K);
  
  sums(1) = sums(1) + sum(Pf-Pc);
  sums(2) = sums(2) + sum((Pf-Pc).^2);
  sums(3) = sums(3) + sum((Pf-Pc).^3);
  sums(4) = sums(4) + sum((Pf-Pc).^4);
  sums(5) = sums(5) + sum(Pf);
  sums(6) = sums(6) + sum(Pf.^2);
  
end

cost = 3*N*nf;
end