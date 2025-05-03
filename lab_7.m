% 1
f = @(t, y) -y./(2*t) + t.^2;
y_exact = @(t) (2/7)*t.^3 + (5/7)*t.^(-1/2);

% 2
function [t, y] = euler_explicit(f, t0, y0, tf, h)
    t = t0:h:tf;
    N = length(t);
    y = zeros(1, N);
    y(1) = y0;
    for i = 1:N-1
        y(i+1) = y(i) + h * f(t(i), y(i));
    end
end

% 3
function [t, y] = euler_cauchy(f, t0, y0, tf, h)
    t = t0:h:tf;
    N = length(t);
    y = zeros(1, N);
    y(1) = y0;
    for i = 1:N-1
        y_predict = y(i) + h * f(t(i), y(i));
        y(i+1) = y(i) + h * 0.5 * (f(t(i), y(i)) + f(t(i+1), y_predict));
    end
end

% 4
function [t, y] = use_ode45(f, t0, y0, tf, h)
    tspan = t0:h:tf;
    [t, y] = ode45(f, tspan, y0);
    y = y(:)';
    t = t(:)';
end

% 5
t0 = 1; tf = 2;
h1 = 0.1; h2 = 0.2;
y_0 = 1;

[t_euler1, y_euler1] = euler_explicit(f, t0, y_0, tf, h1);
[t_euler2, y_euler2] = euler_explicit(f, t0, y_0, tf, h2);

[t_eulerc1, y_eulerc1] = euler_cauchy(f, t0, y_0, tf, h1);
[t_eulerc2, y_eulerc2] = euler_cauchy(f, t0, y_0, tf, h2);

[t_rk1, y_rk1] = use_ode45(f, t0, y_0, tf, h1);
[t_rk2, y_rk2] = use_ode45(f, t0, y_0, tf, h2);

y_ex_1 = y_exact(t_euler1);
y_ex_2 = y_exact(t_euler2);
y_ex_ode1 = y_exact(t_rk1);

figure
plot(t_euler1, y_ex_1, 'k-', 'LineWidth', 2), hold on
plot(t_euler1, y_euler1, 'ro-', 'DisplayName', 'Эйлер явный, h=0.1')
plot(t_eulerc1, y_eulerc1, 'gs-', 'DisplayName', 'Эйлер-Коши, h=0.1')
plot(t_rk1, y_rk1, 'b^-', 'DisplayName', 'ode45, h=0.1')
legend('Точное решение', 'Эйлер', 'Эйлер-Коши', 'ode45')
xlabel('t'), ylabel('y(t)'), grid on
title('Графики решения задачи Коши')

% 6
eps_euler  = max(abs(y_ex_1 - y_euler1));
eps_eulerc = max(abs(y_ex_1 - y_eulerc1));
eps_rk     = max(abs(y_ex_ode1 - y_rk1));

p_euler  = 1;
p_eulerc = 2;
p_rk     = 4;

eps_runge_euler  = max(abs(y_euler1(1:2:end)  - y_euler2)  / (2^p_euler  - 1));
eps_runge_eulerc = max(abs(y_eulerc1(1:2:end) - y_eulerc2) / (2^p_eulerc - 1));
eps_runge_rk     = max(abs(y_rk1(1:2:end)     - y_rk2)     / (2^p_rk     - 1));

fprintf('Максимальная погрешность (Эйлер): %e\n', eps_euler);
fprintf('Погрешность по правилу Рунге (Эйлер): %e\n\n', eps_runge_euler);

fprintf('Максимальная погрешность (Эйлер-Коши): %e\n', eps_eulerc);
fprintf('Погрешность по правилу Рунге (Эйлер-Коши): %e\n\n', eps_runge_eulerc);

fprintf('Максимальная погрешность (ode45): %e\n', eps_rk);
fprintf('Погрешность по правилу Рунге (ode45): %e\n\n', eps_runge_rk);
