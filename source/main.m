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

%% Problem 1.A: Time-Domain Signals
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

%% Problem 1.B: Frequency Response
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

%% Problem 2.B

tau2 = tau;              

f_L = 0.1;                
w_L = 2*pi*f_L;          

fs_ADC = 1e-4;            
f_Nyq = fs_ADC/2;        

% Choose g(t) to shift 100 MHz down to 10 kHz
f_LF = 1e-5;               
f_g = f_L - f_LF;         
w_g = 2*pi*f_g;            
w_LF = 2*pi*f_LF;          

M = 1 / sqrt(1 + (w_L*tau2)^2);
phi = atan(w_L*tau2);

t_rf = 0:0.02:200;        

x_t = 1 + cos(w_L*t_rf);
f_t = 1 + M*cos(w_L*t_rf - phi);
g_t = cos(w_g*t_rf);
y_t = f_t .* g_t;

figure(3);
subplot(2,1,1);
plot(t_rf, x_t);
title('x(t) = 1 + cos(\omega_L t)');
xlabel('Time (ns)');
ylabel('x(t)');
grid on;

subplot(2,1,2);
plot(t_rf, f_t);
title('f(t)');
xlabel('Time (ns)');
ylabel('f(t)');
grid on;

figure(4);
subplot(2,1,1);
plot(t_rf, g_t);
title('g(t) = cos(\omega_g t)');
xlabel('Time (ns)');
ylabel('g(t)');
grid on;

subplot(2,1,2);
plot(t_rf, y_t);
title('y(t) = f(t)g(t)');
xlabel('Time (ns)');
ylabel('y(t)');
grid on;

dt_fft = 1;                  
t_fft = 0:dt_fft:(1e7 - dt_fft); 

f_fft = 1 + M*cos(w_L*t_fft - phi);
g_fft = cos(w_g*t_fft);
y_fft = f_fft .* g_fft;

N = length(t_fft);
fs_fft = 1/dt_fft;         

Y = fft(y_fft)/N;

Nhalf = floor(N/2);
freq = (0:Nhalf)*(fs_fft/N);    
Y_adjusted = Y(1:Nhalf+1);

Ymag = 2*abs(Y_adjusted);
Ymag(1) = abs(Y_adjusted(1));         

freq_kHz = freq * 1e6;

figure(5);
plot(freq_kHz, Ymag);
xlim([0 100]);
title('Magnitude of Y(j\omega)');
xlabel('Frequency (kHz)');
ylabel('|Y(j\omega)|');
grid on;

%% Problem 2.C

x_fft = 1 + cos(w_L*t_fft);          
f_fft = 1 + M*cos(w_L*t_fft - phi);  
g_fft = cos(w_g*t_fft);              

c_fft = x_fft .* g_fft;              
y_fft = f_fft .* g_fft;              

C = fft(c_fft)/N;
Y = fft(y_fft)/N;

C_adjusted = C(1:Nhalf+1);
Y_adjusted = Y(1:Nhalf+1);

Cmag = 2*abs(C_adjusted);
Ymag = 2*abs(Y_adjusted);

Cmag(1) = abs(C_adjusted(1));
Ymag(1) = abs(Y_adjusted(1));

Cphase = angle(C_adjusted);
Yphase = angle(Y_adjusted);

[~, idx_LF] = min(abs(freq - f_LF));

A_copy = Cmag(idx_LF);
phase_copy = Cphase(idx_LF);

A_f = Ymag(idx_LF);
phase_f = Yphase(idx_LF);

M = A_f / A_copy;

delta_phase = angle(exp(1j*(phase_f - phase_copy)));

phi = -delta_phase;

tau_from_phase = tan(phi) / w_L;

fprintf('M = %.4f\n', M);

fprintf('Phase difference = %.4f rad\n', delta_phase);

fprintf('Lifetime from phase = %.4f ns\n', tau_from_phase);
fprintf('True lifetime = %.4f ns\n', tau2);

%% Problem 3: Phasor Analysis

i_info = load("FLIMhistogram.mat");

t_limit = (0:97) * (TL/98);   % 1x98, matches 3rd dim of H
H = i_info.FLIMhistogram;
max(H)
plot_phasor(TL, t_limit, H);