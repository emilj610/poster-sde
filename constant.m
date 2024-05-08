clear all, clc, close all

S0 = 10; % initial stock price
K = 4; % Strike price
r = 1; % risk-free rate
T = 1; % time to expiration
sigma = 0.3; % true volatility
mis_sigma = 1.2; % mis-specified volatility

M = 1e6; % number Monte Carlo sims
N = 1e2; % number of timesteps
dt = T/N;
randn("state",0);

S = zeros(M,1);
mis_S = zeros(M,1);

for i = 1:M
    %generate the Brownian motion
    dW = sqrt(dt) * randn(N,1);
    Wt = sum(dW);
    
    S(i) = S0*exp((r-sigma^2/2)*T+sigma*Wt);
    mis_S(i) = S0*exp((r-mis_sigma^2/2)*T+mis_sigma*Wt);
end

v = exp(-r*T) * mean(max(S-K,0));
mis_v = exp(-r*T) * mean(max(mis_S-K,0));

disp(v)
disp(mis_v)

