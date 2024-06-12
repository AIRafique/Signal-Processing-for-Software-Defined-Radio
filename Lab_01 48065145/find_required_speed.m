function v_required = find_required_speed(theta, xt, yt, y0)
    g = 10; % Acceleration due to gravity (m/s^2)

    % Convert launch angle from degrees to radians
    theta_rad = deg2rad(theta);

    % Compute the horizontal distance to the target (range)
    R = xt;

    % Compute the required initial velocity using the modified range formula
    v_required = sqrt(((R * g - (yt - y0)^2) / sin(2 * theta_rad)));
end