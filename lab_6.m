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

fprintf("\nЗадание 2\n");
a = 0;
b = 1.5;
epsilon = 1e-5;

f = @(x) (2*x).^3 .* cos(x);

fprintf("[    a   ,          b   ] = значение интеграла\n");
function integral = SimpsonRuleAdaptivePrint(f, a, b, epsilon)
  c = (a + b) / 2;
  h = b - a;

  fa = f(a);
  fb = f(b);
  fc = f(c);

  % Однократное вычисление интеграла с шагом h
  I1 = (h / 6) * (fa + 4*fc + fb);

  d = (a + c) / 2;
  e = (c + b) / 2;

  fd = f(d);
  fe = f(e);

  % Разбиение на два подинтервала
  I2 = (h / 12) * (fa + 4*fd + 2*fc + 4*fe + fb);

  R = abs(I2 - I1) / 15;

  if R < epsilon
    integral = I2;
    fprintf("[%.6f,\t%.6f] = %.5f\n", a, b, integral);
  else
    integral_left = SimpsonRuleAdaptivePrint(f, a, c, epsilon/2);
    integral_right = SimpsonRuleAdaptivePrint(f, c, b, epsilon/2);
    integral = integral_left + integral_right;
  end
end

total_integral = SimpsonRuleAdaptivePrint(f, a, b, epsilon);
fprintf("\nИтоговое значение интеграла: %.5f\n", total_integral);
