clear all, clc, close all

S0 = 60; % initial stock price
K = 20; % Strike price
T = 4; % time to expiration

%time dependent r and sigma:
r = @(t) 0.1; % risk-free rate (constant)

%true volatility
a1 = 0.3;
sigma = @(t) a1*(t-(T/2)).^2; 

a2s = [0.1, 0.2, 0.4, 0.5, 0.6];

all_v = zeros(1,length(a2s));
all_mis_v = zeros(1,length(a2s));
tracking_errs = zeros(1,length(a2s));


for j = 1:length(a2s)

%mis-specified volatility
a2 = a2s(j);
sigma_mis = @(t) a2*(t-(T/2.2)).^2;

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
    f = @(t) r(t) - 1/2*sigma(t).^2;
    f_mis = @(t) r(t) - 1/2*sigma_mis(t).^2;
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

disp(mean(S(:,end)));
disp(mean(S_mis(:,end)))

v =  exp(-III) * max(S(:,end)-K,0);
mis_v = exp(-III) * max(S_mis(:,end)-K,0);

all_v(j) = mean(v);
all_mis_v(j) = mean(mis_v);

std(j) = sqrt(var(mis_v));

tracking_errs(j) = all_mis_v(j) - all_v(j);

end

%%

disp(all_v)
disp(all_mis_v)

%% table

T = table(a2s', all_mis_v', tracking_errs', std', 'VariableNames', {'a', 'Mispecified Value', 'Tracking Error', 'Standard Deviation'})

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

figure()
hold on
plot(t_span,sigma(t_span),'r')
plot(t_span,sigma_mis(t_span),'b')


