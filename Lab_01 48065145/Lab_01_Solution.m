clear
clc

%% Task 1
% the self-defined functions are at the end of the file
N = 6;
A_a = mymat_a(N)
A_b = mymat_b(N)
A_c = mymat_c(N)


%% Task 2
% a)
e = randn(1000, 1)*2;
mean(e)
std(e)

% b)
x = rand(1000, 1)*10;
y = 2.4*x + e;
% note x and y are both column vectors

% c)
figure
scatter(x, y)
[~, idx] = sort(y);
idx(end-9)
y(idx(end-9))
hold
plot(x(idx(end-9)), y(idx(end-9)), 'm+')

% d)
beta = inv(x'*x)*x'*y 
xx = 0:0.1:10;
plot(xx, beta*xx, 'r', 'LineWidth', 2)

% e)
x\y

% f)
y = 1.3 + 1.6*x - 0.5*x.^2 + e;
figure
scatter(x, y)
hold
beta = (x.^(0:2))\y
yy = beta(1) + beta(2)*xx + beta(3)*xx.^2;
plot(xx, yy, 'r', 'LineWidth', 2)


%% Task 3
% a)
r1 = roots([1 0 0 0 -1])
poly(r1(2:3)) % pick up proper indices for complex-conjugated roots

% b)
r2 = roots([1 -6 11 -6])

% c)
r3 = roots([1 2.8 1.7 0.9 2.3 2.3 1])
poly(r3(2:3)) % pick up proper indices for complex-conjugated roots
poly(r3(5:6)) % pick up proper indices for complex-conjugated roots


%% Task 4
% a)
g = 10;
v = 5;
theta = pi/4;
t = 0:1e-2:1;
x = v*cos(theta)*t;
y = 2 + v*sin(theta)*t - 0.5*g*t.^2;
figure
plot(x, y)

% b), using a graphic approach
xt = 10;
yt = 3.5;
theta = pi/4;
v = 5:0.1:15;
y_v = zeros(size(v)); % an empty vector to store the y coordinates at xt for varied v
figure
plot(xt, yt, 'ro')
hold
for k = 1:length(v)
    t = 0:1e-2:xt/(v(k)*cos(theta));
    x = v(k)*cos(theta)*t;
    y = 2 + v(k)*sin(theta)*t - 0.5*g*t.^2;
    plot(x, y)
    y_v(k) = y(end);
end
[~, idx_v] = min(abs(y_v - yt));
v(idx_v)

% c)
% solve the math equations first, e.g. to eliminate theta
% then apply the graphic approach
t = 0:1e-2:10;
v = ((xt^2 + (yt - 2 + 0.5*g*t.^2).^2)./t.^2).^0.5;
figure
plot(t, v)
[v_min, idx] = min(v);
v_min
t_min = t(idx);
theta_min = acos(xt/v_min/t_min);
theta_min/pi*180 % in degree

% d)
tt = 0:1e-2:t_min;
x = v_min*cos(theta_min)*tt;
y = 2 + v_min*sin(theta_min)*tt - 0.5*g*tt.^2;
figure
plot(x, y, 'color', [1 1 1]) % for fixing the plot area
hold
for k = 1:length(tt)
    plot(x(k), y(k), 'b.')
    drawnow
    M(k) = getframe;
end
movie(M, 1, 100)


%% Task 5
% a)
r = [0.3, 1.8, 2.2, 2.5, 2.7];
N = 50;
figure
for p = 1:5
    x = zeros(N+1, 1);
    x(1) = 0.1;
    for k = 1:N
        x(k+1) = x(k) + r(p)*x(k) - r(p)*x(k)^2;
    end
    subplot(2,3,p)
    plot(0:N, x)
    title(['r = ' num2str(r(p))])
end

% b)
figure
for p = 1:3
    x = zeros(N+1, 2);
    x(1,:) = [0 1] + 1e-3;
    for k = 1:N
        x(k+1,:) = x(k,:) + r(p)*x(k,:) - r(p)*x(k,:).^2;
    end
    subplot(2,3,p)
    plot(0:N, x(:,1))
    title(['r = ' num2str(r(p))])
    subplot(2,3,p+3)
    plot(0:N, x(:,2))
end

% c)
figure
hold
for r = 1.5:0.01:3
    N = 200;
    x = zeros(N+1, 1);
    x(1) = 0.1;
    for k = 1:N
        x(k+1) = x(k) + r*x(k) - r*x(k)^2;
    end
    plot(r*ones(100,1), x(102:end), '.')
end


%% self-defined functions for Task 1
% method a using loop
function A = mymat_a(N)
A = zeros(N);
for k = 1:N
    for m = 1:N
        A(k, m) = N*(m-1) + k - 1;
    end
end
end

% method b using reshape
function A = mymat_b(N)
A = reshape(0:N^2-1, N, []);
end

% method c using matrix computation
function A = mymat_c(N)
A = (0:N-1)'*ones(1,N) + ones(N,1)*N*(0:N-1);
end
