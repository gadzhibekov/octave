fprintf('Задание 1\n');
% 1.1
f = @(x) tan(x);
a = 0.5;
b = 1.5;
n = 2;
m = 2;
epsilon = 1e-5;

function y = lagrange_poly(x, nodes, values)
    n = length(nodes);
    y = 0;
    for i = 1:n
        L = 1;
        for j = 1:n
            if j != i
                L = L * (x - nodes(j)) / (nodes(i) - nodes(j));
            end
        end
        y = y + values(i) * L;
    end
end

while true
    % 1.2
    nodes = linspace(a, b, n + 1);
    values = f(nodes);

    x_plot = linspace(a, b, 1000);
    f_plot = f(x_plot);
    L_plot = arrayfun(@(x) lagrange_poly(x, nodes, values), x_plot);

    % figure;
    % plot(x_plot, f_plot, 'b-', 'LineWidth', 2, 'DisplayName', 'tan(x)');
    % hold on;
    % plot(x_plot, L_plot, 'r--', 'LineWidth', 2, 'DisplayName', sprintf('L(x) при n=%d', n));
    % scatter(nodes, values, 100, 'k', 'filled', 'DisplayName', 'Узлы интерполяции');
    % xlabel('x'); ylabel('y');
    % title(sprintf('Интерполяция tan(x) (n=%d)', n));
    % legend('show'); grid on; hold off;

    % 1.3
    df = @(x) 1 ./ (cos(x).^2);
    h = 1e-5;
    dL = @(x) (lagrange_poly(x + h, nodes, values) - lagrange_poly(x - h, nodes, values)) / (2 * h);

    % 1.4
    x_mid = (a + b) / 2;
    R = abs(df(x_mid) - dL(x_mid));
    fprintf('n = %d,\tR = %.8f\n', n, R);

    % 1.5
    if R > epsilon
        n = n + 1;
    else
        % 1.6
        R_plot = arrayfun(@(x) abs(df(x) - dL(x)), x_plot);
        figure;
        plot(x_plot, R_plot, 'm-', 'LineWidth', 2);
        xlabel('x'); ylabel('R(x)');
        title(sprintf('Погрешность первых производных (n=%d)', n));
        grid on;
        fprintf('Точность первых производных достигнута при n = %d\n\n', n);
        break;
    end
end

while true
    nodes = linspace(a, b, m + 1);
    values = f(nodes);

    % 1.7
    d2f = @(x) 2 * sin(x) ./ (cos(x).^3);
    
    h = 1e-5;
    d2L = @(x) (lagrange_poly(x + h, nodes, values) - 2*lagrange_poly(x, nodes, values) + lagrange_poly(x - h, nodes, values)) / (h^2);

    % 1.8
    R2 = abs(d2f(x_mid) - d2L(x_mid));
    fprintf('m = %d,\tR2 = %.8f\n', m, R2);

    % 1.9
    if R2 > epsilon
        m = m + 1;
    else
        nodes = linspace(a, b, m);
        values = f(nodes);
        
        % 1.10
        R2_plot = arrayfun(@(x) abs(d2f(x) - d2L(x)), x_plot);
        figure;
        plot(x_plot, R2_plot, 'g-', 'LineWidth', 2);
        xlabel('x'); ylabel('R2(x)');
        title(sprintf('Погрешность вторых производных (m=%d)', m));
        grid on;
        fprintf('Точность вторых производных достигнута при m = %d\n\n', m);
        break;
    end
end

fprintf('Задание 2\n');
% 2.1
a = 0;
b = 1.5;
epsilon = 1e-5;

f = @(x) (2*x).^3 .* cos(x);

i = 0;
a_i = a;
b_i = b;
h_i = b_i - a_i;

function integral = SimpsonRule(f, a, b, h)
    integral = (h/6) * (f(a) + 4*f((a + b) / 2) + f(b));
end

total_integral = 0;

while true
    % 2.2
    I_h_i = SimpsonRule(f, a_i, b_i, h_i);
    
    % 2.3
    h_i_plus_1 = h_i / 2;

    % 2.4
    b_i_plus_1 = (a_i + b_i) / 2;

    I_h_i_plus_1 = (h_i_plus_1 / 6) * (f(a_i) + ...
        4*f(a_i + 0.5*h_i_plus_1) + ...
        2*f(a_i + h_i_plus_1) + ...
        4*f(b_i_plus_1) + ...
        f(b_i));

    % 2.5
    R_i = abs(I_h_i - I_h_i_plus_1) / 15;

    % 2.6
    if R_i > (epsilon * h_i) / (b - a)
        b_i = b_i_plus_1;
    else
        total_integral += I_h_i_plus_1;
        a_i = b_i_plus_1;
        b_i = b;
    end

    if a_i >= b_i
        break;
    end

    h_i = b_i - a_i;
end

fprintf("Значение интеграла: %d\n", total_integral);