function x = logistic_map(r, x0, num_iterations)
    % Initialize array to store iterates
    x = zeros(1, num_iterations + 1);
    
    % Set initial condition
    x(1) = x0;
    
    % Iterate to compute x[k]
    for k = 1:num_iterations
        x(k + 1) = r * x(k) * (1 - x(k));
    end
end