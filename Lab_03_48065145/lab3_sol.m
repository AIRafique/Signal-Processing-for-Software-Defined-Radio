clear
clc

%% Task 1
% a)
N = 256;
x = zeros(N, 1);
x(1:8) = 1;

n = 0:N-1;
X1 = (exp(-1i*2*pi/length(x)).^(n'*n))*x; % via matrix computation

X2 = fft(x);

figure
subplot(221)
stem(n, real(X1))
subplot(222)
stem(n, real(X2), 'r')
subplot(223)
stem(n, imag(X1))
subplot(224)
stem(n, imag(X2), 'r')

% b)
f1 = 6.5e3;
f2 = 7e3;
f3 = 9e3;
fs = 32e3;

n1 = 0:15;
x1 = cos(2*pi*f1*n1/fs) + cos(2*pi*f2*n1/fs) + cos(2*pi*f3*n1/fs);
X1 = fft(x1);

x2 = x1;
X2 = fft(x1, 256);

n3 = 0:255;
x3 = cos(2*pi*f1*n3/fs) + cos(2*pi*f2*n3/fs) + cos(2*pi*f3*n3/fs);
X3 = fft(x3);

figure
subplot(311)
stem(n1, abs(X1))
xlim([0 16])
subplot(312)
stem(0:255, abs(X2), 'r')
xlim([0 256])
subplot(313)
stem(n3, abs(X3), 'm')
xlim([0 256])

% c)
% In order to clearly distinguish two adjacent frequency components,
% we can let the main lobe of the rectangular window equal the frequency gap.
% The former is 2*pi/N, or fs/N, and the latter is 500 Hz
N = round(fs/(f2-f1));
n = 0:N-1;
x = cos(2*pi*f1*n/fs) + cos(2*pi*f2*n/fs) + cos(2*pi*f3*n/fs);
X = fft(x, 256); % use zero padding to increase the sampling density of the spectrum
figure
plot(abs(X))
% N can be further reduced to determine the "minimum" length of 
% the sampled data when adjacent frequency components can still be
% distinguished


%% Task 2
% a)
N = 16;
wvtool(rectwin(N), bartlett(N), hann(N), hamming(N), blackman(N))
N = 32;
wvtool(rectwin(N), bartlett(N), hann(N), hamming(N), blackman(N))
N = 64;
wvtool(rectwin(N), bartlett(N), hann(N), hamming(N), blackman(N))

% b)
f1 = 6e3;
f2 = 6.5e3;
f3 = 8e3;
fs = 32e3;
t = 0:1/fs:0.1;
x = cos(2*pi*f1*t) + cos(2*pi*f2*t) + cos(2*pi*f3*t);
% filterDesigner
% You can generate MATLAB codes from the filterDesigner, through
% 'File'-->'Generate MATLAB Code'-->'Filter Design Function'
% Alternatively,
fc1 = (f1+f2)/2;
fc2 = (f2+f3)/2;
M = 100;
win = hamming(M+1);
B = fir1(M, [fc1 fc2]/(fs/2), 'stop', win);

% c)
figure
stem(0:M, B)
figure
freqz(B)

% d)
y = filter(B, 1, x);
X = fft(x, 4096);
Y = fft(y, 4096);
figure
subplot(221)
plot(t, x)
xlim([0 0.01])
subplot(222)
plot(t, y, 'r')
xlim([0 0.01])
subplot(223)
plot(0:4095, abs(X))
xlim([0 2048])
set(gca, 'XTick', 0:512:2048, 'XTickLabel', num2cell(0:4e3:16e3))
subplot(224)
plot(0:4095, abs(Y), 'r')
xlim([0 2048])
set(gca, 'XTick', 0:512:2048, 'XTickLabel', num2cell(0:4e3:16e3))

% e)
% e.g. via cascading of the filter to increase the stopband attenuation
h2 = conv(B, B);
figure
freqz(h2)
y2 = filter(h2, 1, x);
Y2 = fft(y2, 4096);
figure
plot(0:4095, abs(Y2), 'm')
xlim([0 2048])
set(gca, 'XTick', 0:512:2048, 'XTickLabel', num2cell(0:4e3:16e3))


%% Task 3
% a)
L = 512;
M = length(B);
N = L+M-1;
P = ceil(length(x)/L);
yy = zeros(1, P*L+M-1);
xx = zeros(1, P*L);
xx(1:length(x)) = x;
H = fft(B, N);
figure
hold
for k = 1 : P
    xk = xx((k-1)*L+1: k*L);
    XK = fft(xk, N);
    YK = XK.*H;
    yk = ifft(YK, N);
    plot((k-1)*L+1:(k-1)*L+N, yk)
    yy((k-1)*L+1:(k-1)*L+N) = yy((k-1)*L+1:(k-1)*L+N) + yk;
end
figure
subplot(211)
plot(yy)
y2 = conv(B, x);
subplot(212)
plot(y2, 'r')
xlim([0 4000])

% b)
% overlap-save
yy = zeros(1, P*L);
zi = zeros(1, length(B)-1);
for k = 1 : P
    xk = xx((k-1)*L+1: k*L);
    [yk, zf] = filter(B, 1, xk, zi);
    yy((k-1)*L+1:k*L) = yk;
    zi = zf;
end
figure
plot(yy, 'g')
