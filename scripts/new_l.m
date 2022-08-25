function [sums, cost] = new_l(l,N,option)
%-------------------------------------------------------
%
% level l, truncated version of the new scheme
%
% option: 1 (without) or 2 (with antithetic treatment)
% 
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


for N1 = 1:10000:N
  N2 = min(10000,N-N1+1);
  
  % initialize u and v
  uf = u0*ones(1,N2);
  uc = uf;

  vf = v0*ones(1,N2);
  vc = vf;
  
  Pc = zeros(1,N2);

  if l==0
    dWf1 = sqrt(hf)*randn(1,N2);
    dWf2 = sqrt(hf)*randn(1,N2);
    
    vf_ = xi + (max(0,vf)-xi)*exp(-0.5*kappa*hf);
    uf = uf + 0.5*hf*(r-xi/2) + 0.5*(1/kappa)*(max(0,vf)-xi) ...
        * (exp(-0.5*kappa*hf)-1);
    vf = (sqrt(max(vf_,0))+0.5*sig*dWf2).^2;
    uf = uf + sqrt(max(vf_,0)).*dWf1+(sig/4)*(dWf1.*dWf2);
    uf = uf + 0.5*hf*(r-xi/2) + 0.5*(1/kappa)*(max(0,vf)-xi) ...
        * (exp(-0.5*kappa*hf)-1);
    vf = xi + (max(0,vf)-xi)*exp(-0.5*kappa*hf);
           
  else
      uf = [uf uf];
      vf = [vf vf];

      for n = 1:nc
          dWc1 = zeros(1,N2);
          dWc2 = zeros(1,N2);
          
          dWf1 = sqrt(hf)*randn(M,N2);
          dWf2 = sqrt(hf)*randn(M,N2);
          
          dWf1 = [dWf1 dWf1(M:-1:1,:)];
          dWf2 = [dWf2 dWf2(M:-1:1,:)];

          for m = 1:M
              dWc1 = dWc1 + dWf1(m,1:N2);
              dWc2 = dWc2 + dWf2(m,1:N2);
                            
              vf_ = xi + (max(0,vf)-xi)*exp(-0.5*kappa*hf);
              uf = uf + 0.5*hf*(r-xi/2) + 0.5*(1/kappa)*(max(0,vf)-xi) ...
                  * (exp(-0.5*kappa*hf)-1);
              
              vf = (sqrt(max(vf_,0))+0.5*sig*dWf2(m,:)).^2;
              uf = uf + sqrt(max(vf_,0)).*dWf1(m,:)+(sig/4)*(dWf1(m,:) ... 
                  .*dWf2(m,:));

              uf = uf + 0.5*hf*(r-xi/2) + 0.5*(1/kappa)*(max(0,vf)-xi) ...
                  * (exp(-0.5*kappa*hf)-1);
              vf = xi + (max(0,vf)-xi)*exp(-0.5*kappa*hf);
          end
                    
          vc_ = xi + (max(0,vc)-xi)*exp(-0.5*kappa*hc);
          uc = uc + 0.5*hc*(r-xi/2) + 0.5*(1/kappa)*(max(0,vc)-xi) ...
              * (exp(-0.5*kappa*hc)-1);
          
          vc = (sqrt(max(vc_,0))+0.5*sig*dWc2).^2;
          uc = uc + sqrt(max(vc_,0)).*dWc1+(sig/4)*(dWc1.*dWc2);
          
          uc = uc + 0.5*hc*(r-xi/2) + 0.5*(1/kappa)*(max(0,vc)-xi) ...
              * (exp(-0.5*kappa*hc)-1);
          vc = xi + (max(0,vc)-xi)*exp(-0.5*kappa*hc);

      end
  end

  Pf  = exp(-r*T)*max(0,exp(uf)-K);
  
  if l>0
    if option == 2
        Pf = 0.5*(Pf(1:N2)+Pf(N2+1:end));
    elseif option == 1
        Pf = Pf(1:N2);
    else
        disp('error: needs option = 1 or 2')
    end
    Pc  = exp(-r*T)*max(0,exp(uc)-K);
  else
    Pc=0;
  end
 
  sums(1) = sums(1) + sum(Pf-Pc);
  sums(2) = sums(2) + sum((Pf-Pc).^2);
  sums(3) = sums(3) + sum((Pf-Pc).^3);
  sums(4) = sums(4) + sum((Pf-Pc).^4);
  sums(5) = sums(5) + sum(Pf);
  sums(6) = sums(6) + sum(Pf.^2);
end
cost = 2*N*nf;
end