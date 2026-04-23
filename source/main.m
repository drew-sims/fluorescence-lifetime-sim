%% Parameters
tau = 1.0;      % ns
TL = 12.5;     % ns
Ts = 0.2;      % ns
t_end = 100;   % ns
sample_period = Ts/1000;

% TODO: Weird behaviour that I cannot explain, had to subtract one sample
% period otherwise it will double count, which makes sense in terms of the
% periodicity of the FFT, though get_excitement are not acting as I expect.
t = 0 : sample_period : (t_end - sample_period);

%% Part A: Time-Domain Signals
x = get_excitation(t, TL);
hf = get_fluorescence_irf(t, tau);

t0 = 6;
sigma = 1;
hd = get_detector_irf(t, t0, sigma); % t0 = 2ns, sigma = 0.5ns

f = ifft(fft(x) .* fft(hf));
d = ifft(fft(f) .* fft(hd));
[y_n, t_n] = get_sampling(t,d,Ts);

figure(1);
subplot(4,1,1); plot(t, x); title('Excitation x(t)'); ylabel('Amp'); grid on;
subplot(4,1,2); plot(t, f); title('Emitted Fluorescence f(t)'); ylabel('Amp'); grid on;
subplot(4,1,3); plot(t, d); title('Detector Output d(t) and y(t)'); ylabel('Amp'); grid on;
subplot(4,1,4); stem(t_n, y_n, 'MarkerSize', 2); title('Digitized y[n]'); xlabel('Time (ns)'); grid on;

%% Part B: Frequency Response
fs = 1/sample_period;
N = length(hf);
freqs = (0:N-1)*(fs/N);

H_theory = tau ./ (1 + 1j * 2*pi*freqs * tau);

H_num = fft(hf) * sample_period;

figure(2);
subplot(2,1,1);
loglog(freqs, abs(H_theory), 'k', freqs, abs(H_num), 'r--');
title('Fluorescence Magnitude Response'); ylabel('|H(f)|'); legend('Theory','FFT'); grid on;
subplot(2,1,2);
semilogx(freqs, angle(H_theory), 'k', freqs, angle(H_num), 'r--');
title('Fluorescence Phase Response'); ylabel('Phase (rad)'); xlabel('Frequency (GHz)'); grid on;

%% Problem 3: Phasor Analysis

i_info = load("FLIMhistogram.mat");

t_limit = (0:97) * (TL/98);   % 1x98, matches 3rd dim of H
H = i_info.FLIMhistogram;
max(H)
plot_phasor(TL, t_limit, H);