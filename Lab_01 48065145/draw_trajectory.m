function draw_trajectory(v, theta, x0, y0)
    g = 10; % Acceleration due to gravity (m/s^2)

    % Compute time of flight (T) using the vertical motion equation
    T = max(roots([-0.5*g, v*sind(theta), y0])); % maximum positive root

    % Time intervals for plotting
    t_intervals = linspace(0, T, 1000);

    % Compute x and y positions at each time interval
    x_positions = x0 + v * cosd(theta) * t_intervals;
    y_positions = y0 + v * sind(theta) * t_intervals - 0.5 * g * t_intervals.^2;

    % Plot the trajectory
    plot(x_positions, y_positions);
    xlabel('Horizontal Distance (m)');
    ylabel('Vertical Distance (m)');
    title('Projectile Trajectory');
    grid on;
end
