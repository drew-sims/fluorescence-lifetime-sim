%% Problem 1

tau = 1.0;      % ns
TL = 12.5;     % ns
Ts = 0.2;      % ns
t_end = 50;   % ns
sample_period = Ts/1000;

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

%% Part C: Parameter Sweep

% Sweep over sigma and t0 values
sigma_values = [0.1, 1.0, 10]; % Example values for sigma
t0_values = [2, 4, 6]; % Example values for t0
results = zeros(length(sigma_values), length(t0_values));

figure(3);

for i = 1:length(sigma_values)
    for j = 1:length(t0_values)
        hd = get_detector_irf(t, t0_values(j), sigma_values(i));
        f = ifft(fft(x) .* fft(hf));
        d = ifft(fft(f) .* fft(hd));

        results(i, j) = max(abs(d)); % Store the maximum output for each parameter combination
        subplot(length(sigma_values), length(t0_values), (i-1)*length(t0_values) + j);
        plot(t, d); % Plot the detector output for current parameters
        title(['t0 = ', num2str(t0_values(j)), ', \sigma = ', num2str(sigma_values(i))]);
        ylabel('Amp'); grid on;
    end
end



%% Problem 3: Phasor Analysis

i_info = load("../input/FLIMhistogram.mat");

t_limit = (0:97) * (TL/98);   % 1x98, matches 3rd dim of H
H = i_info.FLIMhistogram;
max(H)
plot_phasor(TL, t_limit, H);