%% TASK 1 A
f_s1 = 1200; f_s2 = 6000; f = 440; duration = 0.02; duration_sound = 2;
t = 0:1/(100*f):duration_sound; x = sin(2*pi*f*t);
t_s1 = 0:1/f_s1:duration; t_s2 = 0:1/f_s2:duration;
x_s1 = sin(2*pi*f*t_s1); x_s2 = sin(2*pi*f*t_s2);
figure (1);
subplot(2,1,1); plot(t_s1, x_s1, 'r.-', 'LineWidth', 1); ylim([-1.2 1.2]);
subplot(2,1,2); plot(t_s2, x_s2, 'b.-', 'LineWidth', 1); ylim([-1.2 1.2]);
sound(x_s1, f_s1); pause(duration_sound); sound(x_s2, f_s2); 
%% TASK 1 B
f_s = 1500; duration = 0.02; duration_sound = 2; t = 0:1/(100*f_s):duration_sound;
x1 = cos(2*pi*700*t); x2 = cos(2*pi*800*t);
figure (2)
subplot(2,1,1); plot(t, x1, 'r.-', 'LineWidth', 1); ylim([-1.2 1.2]); xlim([0 0.02]);
subplot(2,1,2); plot(t, x2, 'b.-', 'LineWidth', 1); ylim([-1.2 1.2]); xlim([0 0.02]);
% sound(x1, f_s); pause(duration_sound); sound(x2, f_s);
x1_new = cos(2*pi*700*t + pi/2); x2_new = cos(2*pi*800*t + pi/2);
figure (3)
subplot(2,1,1); plot(t, x1_new, 'r.-', 'LineWidth', 1); ylim([-1.2 1.2]); xlim([0 0.02]);
subplot(2,1,2); plot(t, x2_new, 'b.-', 'LineWidth', 1); ylim([-1.2 1.2]); xlim([0 0.02]);
% sound(x1_new, f_s); pause(duration_sound); sound(x2_new, f_s);
%% TASK 1 C
A = 128; f0 = 1200; theta = pi/4; t_range = [0 0.001]; fs1 = 12000; fs2 = 42000;
t_continuous = linspace(t_range(1), t_range(2), 10000); x_continuous = A * sin(2*pi*f0*t_continuous + theta);
t1 = linspace(t_range(1), t_range(2), fs1*(t_range(2)-t_range(1))); x_sampled1 = A * sin(2*pi*f0*t1 + theta);
t2 = linspace(t_range(1), t_range(2), fs2*(t_range(2)-t_range(1))); x_sampled2 = A * sin(2*pi*f0*t2 + theta);
x_sampled2_int = int8((x_sampled2/A)*127);
figure; plot(t_continuous, x_continuous, 'k-', 'LineWidth', 1); hold on;
stem(t1, x_sampled1, 'r', 'LineWidth', 1, 'Marker', 'o', 'MarkerFaceColor', 'r'); 
stem(t2, x_sampled2, 'b', 'LineWidth', 1, 'Marker', 's', 'MarkerFaceColor', 'b');
legend('Continuous Signal', 'Sampled at 12 kHz', 'Sampled at 42 kHz', 'Location', 'best');
xlabel('Time (s)'); ylabel('Amplitude'); title('Comparison of Sampled Signals'); grid on; hold off;
sound(x_sampled1, fs1); pause(0.5); sound(x_sampled2, fs2); pause(0.5);


