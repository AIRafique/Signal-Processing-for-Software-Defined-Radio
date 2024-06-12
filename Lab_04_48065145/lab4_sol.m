clear
clc

%% Task 1
% a)
n = 0:127;
x = cos(pi/36*n) + cos(1.5*pi/36*n);
X = fft(x, 4096); % apply zero padding
figure
subplot(211)
stem(n, x)
subplot(212)
plot(0:4095, abs(X))

% b)
x1 = zeros(1, 128);
x1(1:4:end) = x(1:4:end);
X1 = fft(x1, 4096);
figure
subplot(211)
stem(n, x1, 'r')
subplot(212)
plot(0:4095, abs(X1), 'r')

% c)
xd = x(1:4:end);
Xd = fft(xd, 4096);
figure
subplot(211)
stem(0:length(xd)-1, xd, 'g')
subplot(212)
plot(0:4095, abs(Xd), 'g')

xu = zeros(1, 512);
xu(1:4:end) = x;
Xu = fft(xu, 4096);
figure
subplot(211)
stem(0:length(xu)-1, xu, 'm')
subplot(212)
plot(0:4095, abs(Xu), 'm')


%% Task 2
% a)
M = 128; % order of the FIR filter
n = -M/2:M/2;
omega_c = 300e3/6e6*2*pi; % convert 300 kHz to normalised angular frequency
beta = 6;
gain = 5;
h_id = sin(omega_c*n)./n/pi; % h_id is the ideal filter
h_id(M/2+1) = omega_c/pi; % replace sin(0)/0
h = h_id.*kaiser(M+1, beta)'*gain;
figure
freqz(h)

% b)
A = 1;
theta = 0;
fc = 200e3;
fs = 1.2e6;
t = 0:1/fs:0.1; % the signal s(t) was sampled for a duration of 0.1 s
s = A*exp(j*2*pi*fc*t + theta);
x = downsample(s, 2);
fs_x = fs/2;
t_x = [0:length(x)-1]/fs_x;
y = upsample(x, 10);
fs_y = fs_x*10;
t_y = [0:length(y)-1]/fs_y;
z = filter(h, 1, y);
fs_z = fs_y;
t_z = t_y;

figure % plot the time-domain sampled signals
subplot(421)
stem(t, real(s))
xlim([0, 1e-4])
subplot(422)
stem(t, imag(s))
xlim([0, 1e-4])
subplot(423)
stem(t_x, real(x), 'r')
xlim([0, 1e-4])
subplot(424)
stem(t_x, imag(x), 'r')
xlim([0, 1e-4])
subplot(425)
stem(t_y, real(y), 'g')
xlim([0, 1e-4])
subplot(426)
stem(t_y, imag(y), 'g')
xlim([0, 1e-4])
subplot(427)
stem(t_z, real(z), 'm')
xlim([0, 1e-4])
subplot(428)
stem(t_z, imag(z), 'm')
xlim([0, 1e-4])

S = fft(s);
f = linspace(-fs/2, fs/2, length(s)+1);
X = fft(x);
f_x = linspace(-fs_x/2, fs_x/2, length(x)+1);
Y = fft(y);
f_y = linspace(-fs_y/2, fs_y/2, length(y)+1);
Z = fft(z);
f_z = f_y;

figure % plot the spectra
subplot(421)
stem(f(2:end), real(fftshift(S)))
subplot(422)
stem(f(2:end), imag(fftshift(S)))
subplot(423)
stem(f_x(2:end), real(fftshift(X)), 'r')
subplot(424)
stem(f_x(2:end), real(fftshift(X)), 'r')
subplot(425)
stem(f_y(2:end), real(fftshift(Y)), 'g')
subplot(426)
stem(f_y(2:end), real(fftshift(Y)), 'g')
subplot(427)
stem(f_z(2:end), real(fftshift(Z)), 'm')
subplot(428)
stem(f_z(2:end), real(fftshift(Z)), 'm')

% c)
fc = 500e3;
s = A*exp(j*2*pi*fc*t + theta);
x = downsample(s, 2);
y = upsample(x, 10);
z = filter(h, 1, y);
S = fft(s);
X = fft(x);
Y = fft(y);
Z = fft(z);
figure % plot the sampled signals
subplot(421)
stem(t, real(s))
xlim([0, 1e-4])
subplot(422)
stem(t, imag(s))
xlim([0, 1e-4])
subplot(423)
stem(t_x, real(x), 'r')
xlim([0, 1e-4])
subplot(424)
stem(t_x, imag(x), 'r')
xlim([0, 1e-4])
subplot(425)
stem(t_y, real(y), 'g')
xlim([0, 1e-4])
subplot(426)
stem(t_y, imag(y), 'g')
xlim([0, 1e-4])
subplot(427)
stem(t_z, real(z), 'm')
xlim([0, 1e-4])
subplot(428)
stem(t_z, imag(z), 'm')
xlim([0, 1e-4])
figure % plot the spectra
subplot(421)
stem(f(2:end), real(fftshift(S)))
subplot(422)
stem(f(2:end), imag(fftshift(S)))
subplot(423)
stem(f_x(2:end), real(fftshift(X)), 'r')
subplot(424)
stem(f_x(2:end), real(fftshift(X)), 'r')
subplot(425)
stem(f_y(2:end), real(fftshift(Y)), 'g')
subplot(426)
stem(f_y(2:end), real(fftshift(Y)), 'g')
subplot(427)
stem(f_z(2:end), real(fftshift(Z)), 'm')
subplot(428)
stem(f_z(2:end), real(fftshift(Z)), 'm')


%% Task 3
% a)
f = 44.1e3;
L = 160;
M = 147;
% single-stage interpolation
f2 = f*L;
fc = f2/L/2;
fp = fc-1e3;
fs = fc+1e3;
Ap = 0.05; % passband ripple
As = 40; % stopband attenuation in dB
% As an example, use Kaiser window to design the FIR filter via filterDesigner
% It is equivalent to the codes below
[N, Wn, beta, Ftype] = kaiserord([fp fs], [1 0], [0.05, db2mag(-As)], f2);
B = fir1(N, Wn, Ftype, kaiser(N+1, beta));
length(B) % the number of coefficients in the FIR filter

% b)
H = fft(B, 2^14);
figure
subplot(211)
plot(0:2^14-1, mag2db(abs(H)))
xlim([0 2^13])
set(gca, 'XTick', [0 2^13], 'XTickLabel', {'0', '\pi'})
subplot(212)
plot(0:2^14-1, mag2db(abs(H)))
xlim([0 2^10])
set(gca, 'XTick', [0, 2^13/L, 2^10], 'XTickLabel', {'0', '\pi/L', '\pi/8'})

% c)
% two-stage interpolation
L1 = 10;
L2 = 16;
f1 = f*L1;
fc1 = f1/L1/2;
fp1 = fc1-1e3;
fs1 = fc1+1e3;
% Use Kaiser window to design the FIR filter via filterDesigner
% It is equivalent to the codes below
[N1, Wn1, beta1, Ftype1] = kaiserord([fp1 fs1], [1 0], [0.025, db2mag(-As)], f1);
B1 = fir1(N1, Wn1, Ftype1, kaiser(N1 + 1, beta1));

fc2 = f2/L2/2;
fp2 = fc2 - 1e3*L1;
fs2 = fc2 + 1e3*L1;
% Use Kaiser window to design the FIR filter via filterDesigner
% It is equivalent to the codes below
[N2, Wn2, beta2, Ftype2] = kaiserord([fp2 fs2], [1 0], [0.025, db2mag(-As)], f2);
B2 = fir1(N2, Wn2, Ftype2, kaiser(N2 + 1, beta2));

length(B1)+length(B2) % total number of coefficients in the two FIR filters

H1 = fft(B1, 2^14);
figure
subplot(211)
plot(0:2^14-1, mag2db(abs(H1)), 'r')
xlim([0 2^13])
set(gca, 'XTick', [0 2^13], 'XTickLabel', {'0', '\pi'})
subplot(212)
plot(0:2^14-1, mag2db(abs(H1)), 'r')
xlim([0 2^10])
set(gca, 'XTick', [0, 2^13/L1, 2^10], 'XTickLabel', {'0', '\pi/L_{1}', '\pi/8'})

H2 = fft(B2, 2^14);
figure
subplot(211)
plot(0:2^14-1, mag2db(abs(H2)), 'm')
xlim([0 2^13])
set(gca, 'XTick', [0 2^13], 'XTickLabel', {'0', '\pi'})
subplot(212)
plot(0:2^14-1, mag2db(abs(H2)), 'm')
xlim([0 2^10])
set(gca, 'XTick', [0, 2^13/L2, 2^10], 'XTickLabel', {'0', '\pi/L_{2}', '\pi/8'})
