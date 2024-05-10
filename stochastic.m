clear all, clc, close all

S0 = 50; % initial stock price
K = 40; % Strike price
T = 4; % time to expiration
r = 0.1; % risk-free rate

%paramters för SDE:n
%for exakt
th1 = 0.25; %rate
lt1 = 0.4; %long term volatility
xhi1 = 0.3; %volatility of volatility

sig0 = 0.5;
mis_sig0 = 0.5;

params = [0.25, 0.2, 0.7; 0.25, 0.2, 0.1; 0.4, 0.4, 0.7; 0.4, 0.4, 0.1; 0.4, 0.7, 0.7; 0.4, 0.7, 0.1];

all_v = zeros(1,length(params));
all_mis_v = zeros(1,length(params));
tracking_errs = zeros(1,length(params));
std = zeros(1,length(params));

for j = 1:length(params)

th2 = params(j,1); %rate
lt2 = params(j,2); %long term volatility
xhi2 = params(j,3); %volatility of volatility


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

v = exp(-r*T) * max(S(:,end)-K,0);
mis_v = exp(-r*T) * max(mis_S(:,end)-K,0);

all_v(j) = mean(v);
all_mis_v(j) = mean(mis_v);

std(j) = sqrt(var(mis_v));

tracking_errs(j) = all_mis_v(j) - all_v(j);

end

%%

disp(all_v)
disp(all_mis_v)

%% table

T = table(params(:,1), params(:,2), params(:,3), all_mis_v', real(tracking_errs)', std', ...
    'VariableNames', {'Rate', 'Long Term Vol', 'Vol of Vol', 'Mispecified Value', 'Tracking Error', 'Standard Deviation'});
disp(T);

%% plots

hold on
S_bar = zeros(1,N+1);
for i = 1:N+1
    S_bar(i) = mean(S(:,i));
end
plot2 = plot(t_span, S_bar, 'r', 'LineWidth', 2); % Highlighted plot

mis_S_bar = zeros(1,N+1);
for i = 1:N+1
    mis_S_bar(i) = mean(mis_S(:,i));
end
plot3 = plot(t_span, mis_S_bar, 'b', 'LineWidth', 2); % Highlighted plot
%legend('Mean of Exact Stock Price', 'Mean of Mispecified Stock Price'); % Correct usage

title('Mean of simulated stock prices')
hold off




