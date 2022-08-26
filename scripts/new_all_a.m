%function [sums, cost] = new_all_a(l,N,typ,option)
%
% l-level estimator using the new method
% This function is used to estimate the final level given RMSE bound epsilon.
%
% inputs:  l = level
%          N = number of paths to be generated
%          typ = 'Call' to compute Euroopean call option and 'S' to compute
%          the price (e^u) itself
%          option = 1 or 2
%                   option 2 uses the antithetic treatment
%                   option 1 is without the antithetic treatment
%
% output: sums(1) = sum(Pf-Pc)
%         sums(2) = sum((Pf-Pc).^2)
%         sums(3) = sum((Pf-Pc).^3)
%         sums(4) = sum((Pf-Pc).^4)
%         sums(5) = sum(Pf)
%         sums(6) = sum(Pf.^2)
%         cost = cost of N samples

function [sums, cost] = new_all_a(l,N,typ,option)
%-------------------------------------------------------
%
% level l > 0
%

typ = lower(typ);

u0  = log(20);   % initial log asset value
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
  
  uf = [uf uf];
  vf = [vf vf];
  
  uc = [uc uc];
  vc = [vc vc];
  
  if l==0
    dWf1 = sqrt(hf)*randn(1,N2);
    dWf2 = sqrt(hf)*randn(1,N2);
    
    A = randn(1,N2).*sqrt(1/12 * hf^2 + 1/12*hf*((dWf1).^2+(dWf2).^2));
    
    A = [A -A];
    
    dWf1 = [dWf1 dWf1];
    dWf2 = [dWf2 dWf2];
    
    vf_ = xi + (max(0,vf)-xi)*exp(-0.5*kappa*hf);
    
    uf = uf + 0.5*hf*(r-xi/2) + 0.5*(1/kappa)*(max(0,vf)-xi) ...
        * (exp(-0.5*kappa*hf)-1);
    
    vf = (sqrt(vf_)+0.5*sig*dWf2).^2;
    
    uf = uf + sqrt(vf_).*dWf1+(sig/4)*(dWf1.*dWf2)+0.5*sig*A;
    
    uf = uf + 0.5*hf*(r-xi/2) + 0.5*(1/kappa)*(vf-xi) ...
        * (exp(-0.5*kappa*hf)-1);
    
    vf = xi + (vf-xi)*exp(-0.5*kappa*hf);
           
  else

      for n = 1:nc
          dWc1 = zeros(1,N2);
          dWc2 = zeros(1,N2);
          
          Ac = zeros(1,N2);
          
          dWf1 = sqrt(hf)*randn(M,N2);
          dWf2 = sqrt(hf)*randn(M,N2);
          
          A = randn(M,N2).*sqrt(1/12 * hf^2 + 1/12*hf*...
                                ((dWf1).^2+(dWf2).^2)); %For Lévy area;
          
          dWf1 = [dWf1 dWf1(M:-1:1,:)];
          dWf2 = [dWf2 dWf2(M:-1:1,:)];
                                   
          A = [A -A(M:-1:1,:)];

          for m = 1:M
              dWc1 = dWc1 + dWf1(m,1:N2);
              dWc2 = dWc2 + dWf2(m,1:N2);
              
              Ac = Ac + A(m,1:N2);  % Computation of coarse Levy Area
              
              % Ninomiya-Victoir + log-ODE
              
              vf_ = xi + (max(0,vf)-xi)*exp(-0.5*kappa*hf);
              
              uf = uf + 0.5*hf*(r-xi/2) + 0.5*(1/kappa)*(max(0,vf)-xi) ...
                  * (exp(-0.5*kappa*hf)-1);
              
              vf = (sqrt(vf_)+0.5*sig*dWf2(m,:)).^2;

              uf = uf + sqrt(vf_).*dWf1(m,:)+(sig/4)*(dWf1(m,:).*dWf2(m,:))+0.5*sig*A(m,:);

              uf = uf + 0.5*hf*(r-xi/2) + 0.5*(1/kappa)*(vf-xi) ...
                  * (exp(-0.5*kappa*hf)-1);

              vf = xi + (vf-xi)*exp(-0.5*kappa*hf);
          end
          
          Ac = Ac + 0.5*(dWf1(1,1:N2).*dWf2(2,1:N2) - dWf2(1,1:N2).*dWf1(2,1:N2)); % Computation of coarse Levy Area
          Ac = [Ac -Ac];
          
          dWc1 = [dWc1 dWc1];
          dWc2 = [dWc2 dWc2];
                    
          vc_ = xi + (max(0,vc)-xi)*exp(-0.5*kappa*hc);
          uc = uc + 0.5*hc*(r-xi/2) + 0.5*(1/kappa)*(max(0,vc)-xi) ...
              * (exp(-0.5*kappa*hc)-1);
          
          vc = (sqrt(vc_)+0.5*sig*dWc2).^2;
          uc = uc + sqrt(vc_).*dWc1+(sig/4)*(dWc1.*dWc2)+0.5*sig*Ac;
          
          uc = uc + 0.5*hc*(r-xi/2) + 0.5*(1/kappa)*(vc-xi) ...
              * (exp(-0.5*kappa*hc)-1);
          vc = xi + (vc-xi)*exp(-0.5*kappa*hc);

      end
  end
  
  
  if  strcmp(typ,'call')
      Pf  = exp(-r*T)*max(0,exp(uf)-K);
      Pc  = exp(-r*T)*max(0,exp(uc)-K);
  elseif strcmp(typ,'s')
      Pf  = exp(uf);
      Pc  = exp(uc);
  else
      disp('error: needs typ = Call or S');
  end
  
  if option == 2
      Pf = 0.5*(Pf(1:N2)+Pf(N2+1:end));
      if l>0
        Pc = 0.5*(Pc(1:N2)+Pc(N2+1:end));
      else
        Pc = 0;
      end
  elseif option == 1
      Pf = Pf(1:N2);
      if l>0
        Pc = Pc(1:N2);
      else
        Pc = 0;
      end
  else
      disp('error: needs option = 1 or 2');
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
