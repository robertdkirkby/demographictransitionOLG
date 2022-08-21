function F=DemogTransOLG_ReturnFn(h,aprime,a,z,agej,r,A,delta,alpha,sigma,psi,eta,Jr,pension,tau,kappa_j,J,warmglowparam1,warmglowparam2,AccidentBeq, eta1,eta2)

KdivL=((r+delta)/(alpha*A))^(1/(alpha-1));
w=A*(1-alpha)*(KdivL^alpha); % wage rate (per effective labour unit)

% Progressive income tax
if agej<Jr
    Income=w*kappa_j*h+r*a; % Income is labor income and capital income
else
    Income=r*a;
end
IncomeTax=0;
if Income>0
    IncomeTax=eta1+eta2*log(Income)*Income;
end
% Note: Have made pensions exempt from income tax.

F=-Inf;

if agej<Jr
    c=(1+r)*a+(1-tau)*w*kappa_j*h-IncomeTax+(1+r)*AccidentBeq-aprime; % Working age
else
    c=(1+r)*a+pension+(1+r)*AccidentBeq-aprime; % Retired
end

if c>0
    F=(c^(1-sigma))/(1-sigma) -psi*(h^(1+eta))/(1+eta);
end

% Warm-glow bequest
if agej==J % Final period
    warmglow=warmglowparam1*(aprime^(1-warmglowparam2))/(1-warmglowparam2);
    F=F+warmglow;
end
% Notice that we have modelled the warm-glow in such a way that you only
% get it if you die in the last period (j=J). But we know there is a 1-sj
% risk of dying every period. So we might prefer to model that we get the
% warm glow bequest if we die at any age. The following commented out two lines
% implement this alternative. [note: need to add sj to inputs of ReturnFn to use it]
% warmglow=warmglowparam1*(aprime^(1-warmglowparam2))/(1-warmglowparam2); % Note: same formula as above
% F=F+(1-sj)*warmglow
% Note: if using this, have to make sure sj=0 for j=J.
% Comment: I am not aware of any study saying which of these two
% alternatives works better empirically.

end
