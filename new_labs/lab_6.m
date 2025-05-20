a = 0;
b = 1.5;
epsilon = 1e-5;

f = @(x) (2*x).^3 .* cos(x);
f_prime = @(x) 6*(2*x).^2 .* cos(x) - (2*x).^3 .* sin(x); % First derivative
f_double_prime = @(x) 12*(2*x) .* cos(x) - 12*(2*x).^2 .* sin(x) - (2*x).^3 .* cos(x); % Second derivative

fprintf("[    a   ,          b   ] = значение интеграла\n");

function total_integral = SimpsonRuleAdaptiveIterative(f, f_double_prime, a, b, global_epsilon)

    interval_list = [a, b]; % Initial interval
    total_integral = 0;

    while ~isempty(interval_list)
        a_current = interval_list(1);
        b_current = interval_list(2);
        interval_list(1:2) = []; % Remove current interval from the list

        h = b_current - a_current;

        % Exit if the interval is too small
        if h <= 0
            continue;
        end

        c = (a_current + b_current) / 2;
        fa = f(a_current);
        fb = f(b_current);
        fc = f(c);

        % Calculate the integral estimates
        I1 = (h / 6) * (fa + 4 * fc + fb); % Single estimate

        % Intermediate points
        d = (a_current + c) / 2;
        e = (c + b_current) / 2;
        fd = f(d);
        fe = f(e);

        I2 = (h / 12) * (fa + 4 * fd + 2 * fc + 4 * fe + fb); % Double estimate

        R = abs(I2 - I1) / 15;

        % Calculate the curvature at the midpoint
        curvature = abs(f_double_prime(c));

        % Adjust tolerance based on curvature
        adjusted_epsilon = global_epsilon / (1 + curvature);
        % adjusted_epsilon = (global_epsilon * h) / 1.5;

        if R < adjusted_epsilon
            total_integral = total_integral + I2;
            fprintf("[%.6f,\t%.6f] = %.5f; Точность: %.10f\n", a_current, b_current, I2, adjusted_epsilon);
        else
            % Inject new subintervals based on curvature
            interval_list = [interval_list; a_current; c; c; b_current]; % Add new intervals
        end
    end

end

total_integral = SimpsonRuleAdaptiveIterative(f, f_double_prime, a, b, epsilon);
fprintf("\nИтоговое значение интеграла: %.5f\n", total_integral);