x = [2.5, 2.75, 3, 3.25, 3.5, 3.75, 4, 4.25, 4.5, 4.75, 5];
y = [6.109, 2.615, -0.157, -2.010, -2.697, -3.615, -3.478, -2.250, 0.193, 2.086, 5.882];


% 1.1 - 1.3
n = length(x);
m_values = 0:n-1;
sigma = zeros(size(m_values));

figure;
hold on;
plot(x, y, 'ro', 'MarkerFaceColor', 'r');

for m = m_values
    coeffs = polyfit(x, y, m);
    
    P_m = polyval(coeffs, x);
    sigma(m + 1) = sqrt(sum((P_m - y).^2) / (n - m));
    
    plot(x, P_m, 'DisplayName', sprintf('P_%d(x)', m));
end

hold off;
xlabel('x');
ylabel('y');
title('Многочлены наилучшего среднеквадратичного приближения');
legend show;


% 1.4
[~, idx] = min(sigma);
r = idx - 1;
coeffs_r = polyfit(x, y, r);

fprintf('Найденный многочлен степени %d:\n', r);
fprintf('P_%d(x) = ', r);
fprintf('%g', coeffs_r(1));
for i = 2:length(coeffs_r)
    fprintf(' + %gx^%d', coeffs_r(i), i-1);
end
fprintf('\n');


% 2.1 - 2.2
c = abs(min(y)) + 1;
Y = log(y + c); 
X = log(x);


% 2.3
coefficients = polyfit(X, Y, 1);
a = coefficients(1);
b = coefficients(2);


% 2.4
figure;
plot(X, Y, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'r', 'DisplayName', 'Преобразованные данные');
hold on;
Y_fitted = a * X + b;
plot(X, Y_fitted, 'g-', 'LineWidth', 2, 'DisplayName', 'Линейная зависимость');
xlabel('ln(X)');
ylabel('ln(Y)');
title('Линейная зависимость после преобразования');
legend show;
grid on;


% 2.5
fprintf('Линейная зависимость имеет вид: Y = %.4f * X + %.4f\n', a, b);