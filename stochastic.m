clear all, clc, close all

S0 = 3; % initial stock price
K = 4; % Strike price
T = 1; % time to expiration

r = 1; % risk-free rate
sigma0 = 0.3; % true volatility
mis_sigma0 = 3; % mis-specified volatility

M = 1e6; % number Monte Carlo sims
N = 1e2; % number of timesteps
dt = T/N;
randn("state",0);

S = zeros(M,1);
mis_S = zeros(M,1);

%% Det är inte klart!

sigma = sigma0*ones(M,1); % S(0) for all realizations
W = zeros(M,1); % W(0) for all realizations
for j=1:N
    dW = sqrt(dt)*randn(M,1);
    S = S + S.*(r*dt+sigma*dW); % processes at next time step
    W = W + dW; % Brownian paths at next step
end

for i = 1:M
    %generate the Brownian motion
    dW = sqrt(dt) * randn(N,1);
    Wt = sum(dW);
    
    S(i) = S0*exp((r-sigma^2/2)*T+sigma*Wt);
    mis_S(i) = S0*exp((r-mis_sigma^2/2)*T+mis_sigma*Wt);
end

v = exp(-r*T) * mean(max(S-K,0));
mis_v = exp(-r*T) * mean(max(mis_S-K,0));

