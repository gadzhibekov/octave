a = 0;
b = 1.5;
epsilon = 1e-5;
f = @(x) (2*x).^3 .* cos(x);

% Вторая производная функции
f2_deriv = @(x) -8*x.^3 .* cos(x) - 48*x.^2 .* sin(x) + 48*x .* cos(x);

fprintf("[    a   ,          b   ] = значение интеграла\n");

function total_integral = OptimizedSimpsonRule(f, f2_deriv, a, b, global_epsilon)
    % Инициализация массива для отрезков
    interval_list = [a, b];
    total_integral = 0;
    
    while ~isempty(interval_list)
        a_current = interval_list(1);
        b_current = interval_list(2);
        interval_list(1:2) = [];
        
        h = b_current - a_current;
        current_epsilon = global_epsilon * h / (b - a);
        
        if h <= 0
            continue;
        end
        
        c = (a_current + b_current) / 2;
        fa = f(a_current);
        fb = f(b_current);
        fc = f(c);
        I1 = (h / 6) * (fa + 4 * fc + fb);
        
        d = (a_current + c) / 2;
        e = (c + b_current) / 2;
        fd = f(d);
        fe = f(e);
        I2 = (h / 12) * (fa + 4 * fd + 2 * fc + 4 * fe + fb);
        
        R = abs(I2 - I1) / 15;
        
        if R < current_epsilon
            total_integral = total_integral + I2;
            fprintf("[%.6f,\t%.6f] = %.5f; Точность: %.10f\n", a_current, b_current, I2, current_epsilon);
        else
            % Оцениваем кривизну на текущем отрезке
            curvature = max(abs(f2_deriv(linspace(a_current, b_current, 5))));
            
            % Если кривизна большая, разбиваем на 2 подотрезка
            % Если кривизна маленькая, разбиваем на 1 подотрезок (т.е. не разбиваем)
            if curvature > 10 % Порог кривизны можно подобрать экспериментально
                interval_list = [interval_list; a_current; c; c; b_current];
            else
                interval_list = [interval_list; a_current; b_current];
            end
        end
    end
end

total_integral = OptimizedSimpsonRule(f, f2_deriv, a, b, epsilon);
fprintf("\nИтоговое значение интеграла: %.5f\n", total_integral);