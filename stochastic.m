clear all, clc, close all

S0 = 50; % initial stock price
K = 4; % Strike price
T = 4; % time to expiration
r = 0.05; % risk-free rate

%paramters för SDE:n
%for exakt
th1 = 0.15; %rate
lt1 = 0.4; %long term volatility
xhi1 = 0.3; %volatility of volatility
sig0 = 0.5;

%for mis-specified
th2 = 0.4; %rate
lt2 = 0.25; %long term volatility
xhi2 = 0.25; %volatility of volatility
mis_sig0 = 1;

M = 1e6; % number Monte Carlo sims
N = 1e2; % number of timesteps

dt = T/N;
t_span = 0:dt:T;
randn("state",0);

t = 0;

S = zeros(M,N+1);
mis_S = zeros(M,N+1);
sig = zeros(M,N+1);
mis_sig = zeros(M,N+1);

S(:,1) = S0;
mis_S(:,1) = S0;
sig(:,1) = sig0;
mis_sig(:,1) = mis_sig0;

for i = 1:N

    t = t + dt;
    
    dW = sqrt(dt) * randn(M,1);
    dW2 = sqrt(dt) * randn(M,1);

    %taking the sigma steps with Euler Maruyama scheme
    sig(:,i+1) = sig(:,i) + th1*(lt1-sig(:,i))*dt + xhi1*sqrt(sig(:,i)).*dW2;
    mis_sig(:,i+1) = mis_sig(:,i) + th2*(lt2-mis_sig(:,i))*dt + xhi2*sqrt(mis_sig(:,1)).*dW2;
    
    %Euler Maruyama scheme for stock price S
    S(:,i+1) = S(:,i) + r*S(:,i)*dt + sig(:,i).*S(:,i).*dW;
    mis_S(:,i+1) = mis_S(:,i) + r*mis_S(:,i)*dt + mis_sig(:,i).*mis_S(:,i).*dW;
end

%%

v = exp(-r*T) * mean(max(S(:,end)-K,0));
mis_v = exp(-r*T) * mean(max(mis_S(:,end)-K,0));

disp(v)
disp(mis_v)

%% plots

hold on
for i = 1:M/1e5-1
    plot(t_span,S(i*1e4,:))
end
plot2 = plot(t_span, S(end,:), 'r', 'LineWidth', 3); % Highlighted plot

for i = 1:M/1e5-1
    plot(t_span,mis_S(i*1e4,:))
end
plot3 = plot(t_span, mis_S(end,:), 'b', 'LineWidth', 3); % Highlighted plot
legend([plot2, plot3], {'Exact Stock Price', 'Mispecified Stock Price'}); % Correct usage

title('Simulated stock prices')
hold off




