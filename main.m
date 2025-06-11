% In this we code we are trying to implement the OTFS transmitter-receiver
% architecture as per the book 'Delay-Doppler Communications' by Emanuele
% Viterbo, Yi Hong and Tharaj Thaj.

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


%% TRANSMITTER SIDE OPERATIONS (OTFS MODULATION)
% At the transmitter side, given the modulation schem, we first generate
% the OTFS frame symbols. Given these symbols, we have to first convert
% them into the time domain and this can be done in two ways: first obtain
% the time-frequncy symbols and then get the delay-time symbols. This can
% be done by perfoming ISFFT (inverse symplectic fast Fourier transform)
% and then applying the inverse Fourier transform w.r.t the delay domain.
% The other way is to apply the IDZT (inverse discrete Zak transform) on
% the delay-Doppler symbol grid to get the delay-time symbols.

% generating the OTFS symbols
mod_size = 4;                                           
symbols_per_frame = N*M;                                % number of symbols per frame
bits_per_frame = symbols_per_frame * log2(mod_size);    % number of bits per frame
tx_bits = randi([0,1], bits_per_frame, 1);              % generating random bits

% considering QAM modulation
tx_symbols = qammod(tx_bits, mod_size, 'gray', 'InputType', 'bit');
X = reshape(tx_symbols, M, N);                          % generating the MxN OTFS delay-Doppler frame
x = reshape(X.', N*M, 1);                               % vectorized version of the frame information symbols

% OTFS Modulation
Fn = dftmtx(N);
Fn = Fn/norm(Fn);                                       % normalized DFT matrices
Fm = dftmtx(M);
Fm = Fm/norm(Fm);

% % Method 1: ISFFT to get the samples in the time-frequency domain, and then
% % IFFT along the delay domain to get the delay-time domain samples
% X_tf = Fm * X * Fn';                                    % Fourier transform w.r.t delay and inverse Fourier transform w.r.t time domain
% X_t = Fm' * X_tf;                                       % Inverse Fourier transform w.r.t delay domain

% Method 2: IDZT
X_t = X * Fn';                                          % IDZT (invesse Fourier transform w.r.t time domain)

% let us now define the pulse-shaping matrix G_tx
t_samples = linspace(0, 1, M);

% Pulse 1: rectangular pulse
g_tx_rect = ones(1, M);
G_tx_rect = diag(g_tx_rect);

% Pulse 2: raised cosine
rolloff = 0.25;
g_tx_rcos = rcosdesign(rolloff, 6, M, 'normal');
G_tx_rcos = diag(g_tx_rcos(21:21+M-1)); % change what samples we are choosing for pulse shaping.. note that 1 h

% choosing the pulse-shaping matrix
G_tx = G_tx_rcos;

% obtaining the symbol vector
s = reshape((G_tx * X_t), N*M, 1);


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


%% GENERATING r BY PASSING THE Tx SIGNAL THROUGH THE CHANNEL
% % Method 1: using the TDL model
% r = zeros(N*M, 1);
% for q = 0 : N*M-1
%     for l = 0 : delay_spread-1
%         if q >= l
%             r(q+1) = r(q+1) + gs(l+1, q+1)*s(q-l+1);
%         end
%     end
% end

% Method 2: using the time-domain channel matrix (G)
r = G * s;

% adding AWGN
Es = mean(abs(qammod(0:mod_size-1, mod_size).^2));      % average QAM symbol energy (i.e. energy per symbol)
SNR_dB = 0:5:30;
SNR = 10.^(SNR_dB/10);
num_SNR = length(SNR_dB);

for idx = 1:num_SNR
    % computing noise power for current SNR
    sigma_w2 = Es / SNR(idx);

    % generating complex Gaussian noise
    noise = sqrt(sigma_w2/2) * (randn(N*M,1) + 1i*randn(N*M,1));

    % adding noise to received signal
    r_noisy = r + noise;
end


%% RECEIVER SIDE OPERATIONS (OTFS DEMODULATION)
% At the receiver end, we receive the continuous signal r(t), which we
% sample at the time instances Ts = T/M. The received vector signal vector
% r is then used to get the delay-time matrix Y_t,  

G_rx_rect = G_tx_rect;
G_rx_rcos = G_tx_rcos;

Y_t = reshape(r, M, N);
Y_t_demod = G_rx_rcos * Y_t;
Y = Y_t_demod * Fn;

y = reshape(Y.', N*M, 1);                               % vectorizing Y


%% OTFS delay-Doppler LMMSE detection
% x_hat = (H' * H+sigma_w_2 * eye(M*N))^(-1)*(H' * y);
% x_hat = qammod(x_hat, mod_size, 'gray');                % QAM demodulation
% 
% 
% s_hat = (G' * G+sigma_w_2 * eye(M*N))^(-1)*(G' * r);    % estimated time domain samples
% X_hat_td = reshape(s_hat, M, N);                
% X_hat = X_hat_td * Fn;
% x_hat = reshape(X_hat.', N*M, 1);
% x_hat = qammod(x_hat, mod_size, 'gray');               % QAM demodulation
% 
% 
% Vectorize received r_t
% r_vec = r_t(1:M*N);  % Make sure we use only M*N samples
% 
% Construct S_matrix again for LMMSE (Heisenberg operator)
% G = S_matrix;  % Already defined earlier
% 
% Estimate transmitted signal in time domain using LMMSE
% s_hat = (G' * G + (10^(-SNR_dB/10)) * eye(M*N)) \ (G' * r_vec);
% 
% Reshape to M x N (time-frequency domain)
% X_hat_tf = reshape(s_hat, M, N);
% 
% Inverse SFFT: Recover estimated Delay-Doppler domain signal
% X_hat_DD = sqrt(M*N) * ifft(fft(X_hat_tf, [], 1), [], 2);
% 
% Demodulate QAM symbols
% x_hat_indices = qamdemod(X_hat_DD(:), M_order, 'UnitAveragePower', true);
% received_bits = de2bi(x_hat_indices, log2(M_order), 'left-msb');
% received_bits = received_bits(:);  % Make into a bit vector



