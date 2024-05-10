function y = BSCH(S,T,K,r,sigma)
normal = @(x) (1+erf(x/sqrt(2)))/2;
d1 = (log(S/K)+(r+.5*sigma^2)*T)/sigma/sqrt(T);
d2 = (log(S/K)+(r-.5*sigma^2)*T)/sigma/sqrt(T);
y = S*normal(d1)-K*exp(-r*T)*normal(d2);
