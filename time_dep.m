clear all, clc, close all

S0 = 30; % initial stock price
K = 4; % Strike price
T = 4; % time to expiration

%time dependent r and sigma:
r = @(t) 0.2*t + 2; % risk-free rate (linjär)
sigma = @(t) (t-(T/2)).^2; % true volatility (parabolisk)
sigma_mis = @(t) 0.2*t + 1; % mis-specified volatility (antar att den är linjär)

M = 1e6; % number Monte Carlo sims
N = 1e2; % number of timesteps
dt = T/N;
t_span = 0:dt:T;
randn("state",0);

S = zeros(M,1);
mis_S = zeros(M,1);

for i = 1:M
    %generate the Brownian motion
    dW = sqrt(dt) * randn(N,1);

    %calculate integrals for stock price
    f = @(t) r(t) - 1/2*sigma(t);
    I = T/6 * ( f(0) + 4*f(T/2) + f(T) ); %simpsons w/ n = 2

    f_mis = @(t) r(t) - 1/2*sigma_mis(t);
    I_mis = T/6 * ( f_mis(0) + 4*f_mis(T/2) + f_mis(T) ); %simpsons w/ n = 2

    II = sigma(t_span(1:end-1))*dW; %Euler-Maruyama scheme
    II_mis = sigma_mis(t_span(1:end-1))*dW; %Euler-Maruyama scheme
    
    %calculates the stock price at t = T
    S(i) = S0*exp(I+II); 
    mis_S(i) = S0*exp(I_mis + II_mis);S = S0*ones(M,1); % S(0) for all realizations
    W = zeros(M,1); % W(0) for all realizations
    for j=1:N
      dW = sqrt(dt)*randn(M,1);
      S = S + S.*(r*dt+sigma*dW); % processes at next time step
      W = W + dW; % Brownian paths at next step
    end
end

III = T/6 * ( r(0) + 4*r(T/2) + r(T) ); %simpsons w/ n = 2

v = exp(-III) * mean(max(S-K,0));
mis_v = exp(-III) * mean(max(mis_S-K,0));

disp(v)
disp(mis_v)

