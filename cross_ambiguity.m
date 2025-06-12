% In this MATLAB code we try to plot the cross-ambiguity function.

clc; clear; close all;


% Parameters
T = 1;
f = linspace(-4, 4, 200);     % Frequency range
t = linspace(0, 2*T, 200);    % Time range

% Meshgrid for evaluating function over f and t
[F, Tt] = meshgrid(f, t);

% Initialize cross-ambiguity function
A = zeros(size(F));

% Compute A(f, t) based on piecewise definition
j = 1i;
idx1 = Tt < T;
idx2 = ~idx1;

% Equation (3.6)
A(idx1) = (j ./ (2 * pi .* F(idx1))) .* (1 - exp(-j * 2 * pi .* F(idx1) .* Tt(idx1)));
A(idx2) = (1 ./ (2 * pi .* F(idx2))) .* (1 - exp(-j * 2 * pi .* F(idx2) .* (2*T - Tt(idx2))));

% Avoid division by zero (handle f=0 separately using L'Hôpital's rule)
A(F == 0) = Tt(F == 0) / T;

% Plot magnitude
figure;
mesh(F, Tt, abs(A));
xlabel('f');
ylabel('t');
zlabel('|A_{g_{tx}, g_{tx}}(f, t)|');
title('Cross-Ambiguity Function Magnitude of Square Pulse');
grid on;


%% Plotting |A(f0, t)| vs. t for fixed f = f0
% Parameters
T = 1;
t = linspace(0, 2*T, 500);  % Time vector
j = 1i;

% Frequencies to evaluate
f_values = [0.5, 1, 2];

% Plot
figure;
hold on;
for f0 = f_values
    A = zeros(size(t));
    
    % Compute A(f0, t) piecewise
    idx1 = t < T;
    idx2 = ~idx1;
    
    A(idx1) = (j ./ (2 * pi * f0)) .* (1 - exp(-j * 2 * pi * f0 .* t(idx1)));
    A(idx2) = (1 ./ (2 * pi * f0)) .* (1 - exp(-j * 2 * pi * f0 .* (2*T - t(idx2))));
    
    plot(t, abs(A), 'DisplayName', ['f = ' num2str(f0)]);
end

xlabel('Time t');
ylabel('|A(f_0, t)|');
title('Cross-Ambiguity vs. Time for fixed f_0');
legend show;
grid on;


%% Plotting |A(f,t0)| vs. f for fixed t = t0
% Parameters
T = 1;
f = linspace(-4, 4, 500);   % Frequency vector
j = 1i;

% Times to evaluate
t_values = [0.25, 0.75, 1.5];

% Plot
figure;
hold on;
for t0 = t_values
    A = zeros(size(f));
    
    idx1 = t0 < T;
    idx2 = ~idx1;
    
    if idx1
        A = (j ./ (2 * pi .* f)) .* (1 - exp(-j * 2 * pi .* f * t0));
    else
        A = (1 ./ (2 * pi .* f)) .* (1 - exp(-j * 2 * pi .* f * (2*T - t0)));
    end
    
    % Handle f = 0 case
    A(f == 0) = t0 / T;
    
    plot(f, abs(A), 'DisplayName', ['t = ' num2str(t0)]);
end

xlabel('Frequency f');
ylabel('|A(f, t_0)|');
title('Cross-Ambiguity vs. Frequency for fixed t_0');
legend show;
grid on;
