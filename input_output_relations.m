% In this simulation we try to observe the simulate the input-output
% relations in different domains, as in: the time domain, the
% time-frequency domain, the time-delay domain and the delay-Doppler
% domian. Across the code, we consider the transmitter and the receiver
% pulse shaping matrix to be identity matrix, meaning that we are using the
% rectangular pulse for pulse-shaping.

clc; clear; close all;

%% OTFS FRAME PARAMETERS
N = 8;                                                  % number of Doppler bins (time-slots)
M = 8;                                                  % number of delay bins (sub-carriers)
delta_f = 15e3;                                         % sub-carrier spacing
T = 1/delta_f;                                          % block duration
fc = 4e9;                                               % carrier frequency, 4GHz
c = 3e8;                                                % speed of electromagnetic waves

% defining the OTFS grid delay and Doppler resolutions
delay_resolution = 1/(M*delta_f);                                      
Doppler_resolution = 1/(N*T);


%% MOBILE USER PARAMETERES AND CHANNEL PARAMETERS
% In this part of the code, we aim to generate the three parameters we need
% for our analysis: complex channel gain (g_i), (integer) delays (l_i) and
% the Doppler taps (k_i).

max_user_speed = 100;                                   % maximum speed of the user in kmph
speed_in_ms = max_user_speed * 1000 / 3600;             % same speed in m/s
max_Doppler = speed_in_ms*fc / c;                       % maximum possible Doppler shift

% maximum normalized Doppler spread
k_max = max_Doppler / Doppler_resolution;

% defining the delays. Note that 'pdp' below refers to the Power Delay Profile, which
% tells us how the power of a signal arriving at the receiver is
% distribuited onver various delays.
delays = [0, 30, 70, 90, 110, 190, 410]*1e-9;
pdp = [0.0, -1.0, -2.0, -3.0, -8.0, -17.2, -20.8];      % pdp in dB
pdp_linear = 10.^(pdp/10);                              % dB to linear scale
pdp_linear = pdp_linear / sum(pdp_linear);              % normalization
taps = length(pdp);

% generating the channel coefficients (assuming Rayleigh fading). Note that
% g_i(i) represents the complex fading coefficient for the i-th path.
g_i = sqrt(pdp_linear).*(sqrt(1/2) * (randn(1, taps) + 1i*randn(1, taps)));

% generating delay taps (assuming integer delay taps)
l_i = round(delays./delay_resolution);

% generating Doppler taps (assuming Jakes spectrum). Jakes spectrum is a
% statistical model used to simulate the Doppler effects based on the
% following assumptions:
% 1. The mobile receiver is surrounded by an infinite number of obstacles
% 2. The incoming wave components are uniformly distribuited in angle.
% 3. The resulting Doppler shifts follow a specific cosine-shifted
% probability distribution.
k_i = (k_max * cos(2*pi*rand(1, taps)));


%% OTFS MODULATION
% as usual, we first generate the symbols according to some modulation
% scheme and place them in the OTFS / delay-Doppler grid
mod_size = 4;                                           
symbols_per_frame = N*M;                                % number of symbols per frame
bits_per_frame = symbols_per_frame * log2(mod_size);    % number of bits per frame
tx_bits = randi([0,1], bits_per_frame, 1);              % generating random bits

% considering QAM modulation
tx_symbols = qammod(tx_bits, mod_size, 'gray', 'InputType', 'bit');
X = reshape(tx_symbols, M, N);                          % generating the MxN OTFS delay-Doppler frame
x = reshape(X.', N*M, 1);                               % vectorized version of the frame information symbols

% let us now perform the OTFS modulation to get the transmitted signal
% vector s. But first, let us get the normalized DFT matrices
Fn = dftmtx(N);
Fn = Fn/norm(Fn);                                       
Fm = dftmtx(M);
Fm = Fm/norm(Fm);

X_tilda = X * Fn';

% let us now define the pulse-shaping matrix G_tx
t_samples = linspace(0, 1, M);

% considering rectangular pulse
g_tx_rect = ones(1, M);
G_tx_rect = diag(g_tx_rect);

% obtaining the symbol vector
s = reshape((G_tx_rect * X_tilda), N*M, 1);


%% DISCRETE_TIME CHANNEL COEFFICIENTS AND MATRIX
z = exp(1i*2*pi/(N*M));
delay_spread = max(l_i) - min(l_i);

% generating discrete-time baseband channel in TDL form
gs = zeros(delay_spread+1, N*M);
for q = 0 : N*M-1
    for i = 1 : taps
        gs(l_i(i)+1, q+1) = gs(l_i(i)+1, q+1) + g_i(i)*z^(k_i(i)*(q-l_i(i)));
    end
end

% generating the discrete-time baseband channel matrix
G = zeros(N*M, N*M);
for q = 0 : N*M-1
    for l = 0:delay_spread
        if q>=l
            G(q+1, q-l+1) = gs(l+1, q+1);
        end
    end
end


%% TIME DOMAIN ANALYSIS
% In this domain, the relation between the received symbol and the
% transmitted symbol vector is given simply as r = Gs + w, where w is an
% AWGN noise vector. Thus, let us for create an AWGN noise vector
w = randn(N*M, 1);

r = G*s + w;                                            % ....eq(1)


%% TIME-FREQUENCY DOMAIN ANALYSIS
% First, let us define the time-frequency domain symbols vectors as y_hat
% and x_hat respectively from r and s as follows:
x_hat = kron(eye(N), Fm) * s;
y_hat = kron(eye(N), Fm) * r;

% substituting the above two in eq(1), we get the time-frequency channel
% matrix H_hat and the transformed noise vector as
H_hat = kron(eye(N), Fm) * G * kron(eye(N), Fm');
w_hat = kron(eye(N), Fm) * w;


%% TIME-DELAY DOMAIN ANALYSIS
% Again, let us define the relationship between the time-delay
% symbols y_tilda and x_tilda respectively from the vectors r and s. But
% first, let define the row-column interleaver matrix / the permuation
% matrix P as follows:
P = zeros(N*M, M*N);
for j = 1:M
    for i = 1:N
        E = zeros(M,N);
        E(i,j) = 1;
        P((j-1)*M+1:j*M, (i-1)*N+1:i*N) = E;
    end
end

% now, we have
x_tilda = P' * s;
y_tilda = P' * r;

% substituing the above two results in eq(1), we get the delay-time channel
% matrix and the transformed noise vector as
H_tilda = P' * G * P;
w_tilda = P' * w;


%% DELAY-DOPPLER DOMAIN ANALYSIS
% Let us define a new matrix A as:
A = kron(eye(M), Fn) * P';

% the relationship between the delay-Doppler symbol vectors y and x is then
% given as follows:
y = A * r;
x = A * s;

% substituing the above results in eq(1) we get the delay-Doppler matrix H
% and the transformed delay-Doppler noise vector as follows:
H = A * G * A';
z = A * w;