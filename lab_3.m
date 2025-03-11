eps = 10^(-4);

x2 =-0.2:0.01:0;

% f1 = cos(x1 + 0.5) + x2 - 0.8;
% f2 = -2 .* x1 + sin(x2) - 1.6;

% 1.1
g1 = -acos(0.8 - x2) - 0.5;
g2 = (1.6 - sin(x2)) ./ (-2);

figure();
plot(x2, g1, x2, g2);
xlabel('x2');
ylabel('x1');
title('Графики функций: g1 и g2');
grid on;
legend("g1", "g2");


% 1.2
phi1 = -cos(x2 + 0.5) + 0.8;
phi2 = (1.6 - sin(x2)) ./ (-2);

% 1.3
x1_0 = -0.2;
x2_0 = 1.9;

J = [
        -sin(x1_0 + 0.5), 0
        0                 , cos(x2_0) / 2
    ];

q = norm(J, 1);

% 1.4
eps_0 = ((1 - q) / q) * eps;


[x, fval, info] = fsolve (@f, [x1_0; x2_0]);
fprintf('Сошлось к %f\n',x);

x_0 = [x1_0, x2_0];



printf("\n\n___________________________________________________________\n");
printf("| n |    x1     |     x2    |   difference     | epsilon_0 |\n");

k1 = 0;
x_prev = x_0;

while true

    k1 += 1;
    x_next = [(1.6 - sin(x_prev(2))) / (-2), -cos(x_prev(1) + 0.5) + 0.8];
    difference = norm(x_next - x_prev);
    x_prev = x_next;

    if (difference <= eps_0)
        printf("| %d | %f | %f |   %f       | %f  |\n", k1, x_next(1), x_next(2), difference, eps_0);
        break;
    end

    printf("| %d | %f | %f |   %f       | %f  |\n", k1, x_next(1), x_next(2), difference, eps_0);
endwhile
printf("___________________________________________________________\n");

% 1.5
printf("Номер итерации на которой была достигнута заданная точность: %d\n", k1);


% 2.1
printf("\n\n___________________________________________________________\n");
printf("| n |    x1     |     x2    |   difference     | epsilon_0 |\n");

x_prev = x_0;
k2 = 0;

while true
    k2 += 1;

    eq_1 = x_prev;
    eq_2 = [-sin(x_prev(1) + 0.5), 1; -2, cos(x_prev(2))];
    eq_3 = [cos(x_prev(1) + 0.5) + x_prev(2) - 0.8; -2 * x_prev(1) + sin(x_prev(2)) - 1.6];

    x_next = eq_1 - inv(eq_2) * eq_3;
    difference = norm(x_next - x_prev);
    x_prev = x_next;

    if (difference <=  eps)
        printf("| %d | %f | %f |   %f       | %f  |\n", k2, x_next(1), x_next(2), difference, eps);
        break
    end

    printf("| %d | %f | %f |   %f       | %f  |\n", k2, x_next(1), x_next(2), difference, eps);

endwhile
printf("___________________________________________________________\n");


% 2.2
printf("Номер итерации на которой была достигнута заданная точность: %d\n", k2);