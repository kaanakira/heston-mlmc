function [sums, cost] = na_milstein_l(l,N,option)
%-------------------------------------------------------
%
% level l, no-area Milstein 
%
% option: 1 (without) or 2 (with antithetic treatment)
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
      uf = uf+(r-0.5*vf)*hf+sqrt(max(0,vf)).*dWf1+0.25*sig*(dWf1.*dWf2);
      
      vf = vf + kappa*(theta-vf)*hf + sig*sqrt(max(0,vf)).*dWf2 ...
          +0.25*(sig^2)*(dWf2.^2-hf);
           
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
                            
              uf = uf + (r-0.5*vf)*hf + sqrt(max(0,vf)).*dWf1(m,:)  ...
                  + 0.25*sig*(dWf1(m,:).*dWf2(m,:));

              vf = vf + kappa*(theta-vf)*hf + sig*sqrt(max(0,vf)) ...
                  .*dWf2(m,:)+ 0.25*(sig^2)*(dWf2(m,:).^2-hf);
          end
                    
          uc = uc + (r-0.5*vc)*hc + sqrt(max(0,vc)).*dWc1  ...
              + 0.25*sig*(dWc1.*dWc2);
          
          vc = vc + kappa*(theta-vc)*hc + sig*sqrt(max(0,vc)).*dWc2 ...
              + 0.25*(sig^2)*(dWc2.^2-hc);

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