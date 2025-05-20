% Параметры задачи
t0 = 1; tf = 2; y0 = 1;
h1 = 0.1; h2 = 0.2;

% Функция правой части
f = @(t,y) -y./(2*t) + t.^2;

% Аналитическое решение
y_exact_fun = @(t) t.^3/3.5 + (5/7)*t.^(-0.5);

% Сетка для h1
t1 = t0:h1:tf; N1 = length(t1);
y_euler1 = zeros(1,N1);
y_cauchy1 = zeros(1,N1);
y_euler1(1) = y0; y_cauchy1(1) = y0;

% Сетка для h2
t2 = t0:h2:tf; N2 = length(t2);
y_euler2 = zeros(1,N2);
y_cauchy2 = zeros(1,N2);
y_euler2(1) = y0; y_cauchy2(1) = y0;

% Явный метод Эйлера и метод Эйлера–Коши
for i = 1:N1-1
    % Эйлер
    y_euler1(i+1) = y_euler1(i) + h1*f(t1(i), y_euler1(i));
    % Эйлер–Коши
    y_temp = y_cauchy1(i) + h1*f(t1(i), y_cauchy1(i));
    y_cauchy1(i+1) = y_cauchy1(i) + h1 * (f(t1(i), y_cauchy1(i)) + f(t1(i+1), y_temp)) / 2;
end

for i = 1:N2-1
    % Эйлер
    y_euler2(i+1) = y_euler2(i) + h2*f(t2(i), y_euler2(i));
    % Эйлер–Коши
    y_temp = y_cauchy2(i) + h2*f(t2(i), y_cauchy2(i));
    y_cauchy2(i+1) = y_cauchy2(i) + h2 * (f(t2(i), y_cauchy2(i)) + f(t2(i+1), y_temp)) / 2;
end

% Решение встроенной функцией ode45
[t_ode, y_ode] = ode45(f, [t0 tf], y0);

% Оценка погрешности по правилу Рунге
p_euler  = 1; p_cauchy = 2;
err_runge_euler  = max(abs(y_euler1(1:2:end)  - y_euler2)/(2^p_euler  - 1));
err_runge_cauchy = max(abs(y_cauchy1(1:2:end) - y_cauchy2)/(2^p_cauchy - 1));

% Точное решение в узлах t1
y_exact1 = y_exact_fun(t1);

% Оценка максимальной погрешности
err_euler1  = max(abs(y_exact1 - y_euler1));
err_cauchy1 = max(abs(y_exact1 - y_cauchy1));
err_ode     = max(abs(y_exact1 - interp1(t_ode, y_ode, t1)));

% Вывод погрешностей
fprintf('Макс. ошибки:\n Эйлер: %g\n Эйлер–Коши: %g\n ode45: %g\n', ...
        err_euler1, err_cauchy1, err_ode);
fprintf('Ошибка Рунге:\n Эйлер: %g\n Эйлер–Коши: %g\n', ...
        err_runge_euler, err_runge_cauchy);

% Построение графиков
figure; hold on; grid on;
plot(t1,    y_exact1,   'k-',  'LineWidth',1.5);
plot(t1,    y_euler1,   'ro--');
plot(t1,    y_cauchy1,  'bs--');
plot(t_ode, y_ode,      'm-.', 'LineWidth',1);
legend('Аналитич.', 'Эйлер', 'Эйлер–Коши', 'ode45');
xlabel('t'); ylabel('y'); title('Сравнение решений');