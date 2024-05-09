clear all, clc, close all

S0 = 50; % initial stock price
K = 4; % Strike price
T = 4; % time to expiration

%time dependent r and sigma:
r = @(t) 0.2*t-0.1; % risk-free rate (linjär)
sigma = @(t) 0.3*(t-(T/2)).^2; % true volatility (parabolisk)
sigma_mis = @(t) 0.2*t + 0.5; % mis-specified volatility (antar att den är linjär)

M = 1e6; % number Monte Carlo sims
N = 1e2; % number of timesteps

dt = T/N;
t_span = 0:dt:T;
randn("state",0);

t = 0;
II = zeros(M,1);
II_mis = zeros(M,1);

S = zeros(M,N+1);
S_mis = zeros(M,N+1);

S(:,1) = S0;
S_mis(:,1) = S0;

for i = 1:N

    t = t + dt;
    
    dW = sqrt(dt) * randn(M,1);
    
    %calculate integral I 
    f = @(t) r(t) - 1/2*sigma(t);
    f_mis = @(t) r(t) - 1/2*sigma_mis(t);
    I = integral(f,0,t);
    I_mis = integral(f_mis,0,t);

    %computes integral II with Euler-Maruyama scheme
    II = sigma(t)*dW;
    II_mis = sigma_mis(t)*dW;

    %computes stock price
    S(:,i+1) = S0*exp(I+II); 
    S_mis(:,i+1) = S0*exp(I_mis + II_mis);

end

III = T/6 * ( r(0) + 4*r(T/2) + r(T) ); %simpsons w/ n = 2

v = exp(-III) * mean(max(S(:,end)-K,0));
mis_v = exp(-III) * mean(max(S_mis(:,end)-K,0));

disp(v)
disp(mis_v)

%% plots

hold on
for i = 1:M/1e5-1
    plot(t_span,S(i*1e4,:))
end
plot2 = plot(t_span, S(end,:), 'r', 'LineWidth', 3); % Highlighted plot

for i = 1:M/1e5-1
    plot(t_span,S_mis(i*1e4,:))
end
plot3 = plot(t_span, S_mis(end,:), 'b', 'LineWidth', 3); % Highlighted plot
legend([plot2, plot3], {'Exact Stock Price', 'Mispecified Stock Price'}); % Correct usage

title('Simulated stock prices')
hold off




