eps = 10^(-4);

% 1.1
x = 0.8:0.01:1;
y = x.^4 .* 3.^x - 2;

figure();
plot(x, y);
xlabel('x');
ylabel('y');
title('График функции y = x^4 * 3^x - 2');
grid on;
legend("f(x)");

a = 0.8;
b = 1;


% 1.2
phi1 = (2 ./ 3.^x).^0.25;
phi2 = log(2 ./ x.^4) ./ log(3);


% 1.3
z = x;
figure();
plot(x, y, x, z, x, phi1, x, phi2);
grid on;
legend("y = f(x)", "y = x", "phi1", "phi2");


% 1.4
mdphi1 = abs(-log(3) ./ (2 .* (8 .^ 0.25) .* (3 .^ (0.25 .* x))));
mdphi2 = abs(-4 ./ (x .* log(3)));
y_1_4 = ones(size(x));                                                      % //  y_1_4 = 1

figure();
plot(x, y_1_4, x, mdphi1, x, mdphi2);
grid on;
legend("y = 1", "mdphi1", "mdphi2");


% 1.5
q = 0.248177; % // q = |-log(3) / (2 * (8 ^ 0.25) * (3 ^ (0.25 * x)))| при x = 1; (mdphi1 возрастает на отрезке [a, b])


% 1.6 - 1.7
x0 = a;
k1 = 0;
check = true;
eps_0 = ((1 - q) / q) * eps;

fprintf("\n\n\nЗадание 1:\n");
fprintf("Номер итерации | Корень  | Разность | Точность \n");
while check
    k1 = k1 + 1;
    x0last = x0;
    x0 = (2 / 3 ^ x0) ^ 0.25;

    if abs(x0 - x0last) <= eps_0
        check = false;
    end

    fprintf("___\n");
    fprintf(" %-14d | %7f | %9f | %9f \n", k1, x0, abs(x0 - x0last), eps_0);
end

% 1.8
fprintf("\nКорень, найденный с точностью не меньше чем 𝜀() \n");
x0


% 2.1
df = ((4 .* x .^ 3) .* (3 .^ x)) + ((log(3) .* x .^ 4) .* (3 .^ x));
figure();
plot(x, df);
grid on;
legend("df");


% 2.2
m = ((4 .* a .^ 3) .* (3 .^ a)) + ((log(3) .* a .^ 4) .* (3 .^ a));
M = ((4 .* b .^ 3) .* (3 .^ b)) + ((log(3) .* b .^ 4) .* (3 .^ b));


% 2.3
alpha = 2 / (M + m);
new_q = (M - m) / (M + m);

% 2.4 - 2.5
new_x0 = a;
check = true;
k2 = 0;
eps_0_new = (1 - new_q) / new_q * eps;

fprintf("\n\n\nЗадание 2:\n");
fprintf("Номер итерации | Корень  | Разность | Точность \n");
while check
    k2 = k2 + 1;
    newx0last = new_x0;
    new_x0 = newx0last - alpha * (newx0last^4 * 3^newx0last - 2);
    if abs(new_x0 - newx0last) <= eps_0_new
        check = false;
    end
    fprintf("___\n");
    fprintf(" %-14d | %7f | %9f | %9f \n", k2, new_x0, abs(new_x0 - newx0last), eps_0_new);
end


% 2.6
fprintf("\nКорень, найденный с точностью не меньше чем 𝜀() \n");
new_x0


% 3.1 - 3.2
check = true;
k3 = 0;
nnew_x0 = a;

fprintf("\n\n\nЗадание 3:\n");
fprintf("Номер итерации | Корень  | Разность | Точность \n");
while check
    nnew_x0_last = nnew_x0; 
    k3 = k3 + 1;
    f_val = nnew_x0_last^4 * 3^nnew_x0_last - 2;
    df_val = ((4 * nnew_x0_last^3) * (3^nnew_x0_last)) + ((log(3) * nnew_x0_last^4) * (3^nnew_x0_last));
    nnew_x0 = nnew_x0_last - f_val / df_val;
    if abs(nnew_x0 - nnew_x0_last) <= eps
        check = false;
    end
    fprintf("___\n");
    fprintf(" %-14d | %7f | %9f | %9f \n", k3, nnew_x0, abs(nnew_x0 - nnew_x0_last), eps);
end


% 3.3
fprintf("\nКорень, найденный с точностью не меньше чем 𝜀() \n");
nnew_x0


% 3.4
fprintf("\nКоличество итераций в трех таблицах (k1, k2, k3) сответственно\n");
k1, k2, k3