clear all, clc, close all

S0 = 100; % initial stock price
K = 110; % Strike price
r = 0.05; % risk-free rate
T = 1; % time to expiration
sigma = 0.5; % true volatility
mis_sigma = 1; % mis-specified volatility

M = 1e6; % number Monte Carlo sims
N = 1e2; % number of timesteps
dt = T/N;
randn("state",0);

its = 10;
prices = zeros(its,1);
values = zeros(its,1);
values_var = zeros(its,1);

price = BSCH(S0,T,K,r,sigma);

mis_sigma = 0.45;

%% mis calculated value
for k = 1:its
mis_S = zeros(M,1);

mis_price = BSCH(S0,T,K,r,mis_sigma);
prices(k) = mis_price;

for i = 1:M
    %generate the Brownian motion
    dW = sqrt(dt) * randn(N,1);
    Wt = sum(dW);

    mis_S(i) = S0*exp((r-mis_sigma^2/2)*T+mis_sigma*Wt);
end

mis_v = exp(-r*T) * max(mis_S-K,0);

values(k) = mean(mis_v);
values_var(k) = var(mis_v);

mis_sigma = mis_sigma + 0.01;
end

%% true value
S = zeros(M,1);
for i = 1:M
    %generate the Brownian motion
    dW = sqrt(dt) * randn(N,1);
    Wt = sum(dW);
    
    S(i) = S0*exp((r-sigma^2/2)*T+sigma*Wt);
end

v = exp(-r*T) * max(S-K,0);


true_var = var(v);
true_value = mean(v);
