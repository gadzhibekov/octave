f = @(x) x .* log(sqrt(x - 2));

a = 2.5;
b = 3.1;

x = a:0.1:b;

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


% 2.1
function y = newton_poly(x, nodes, values)   
    n = length(nodes);
    
    div_diff = zeros(n, n);
    div_diff(:, 1) = values; % // div_diff(:, 1) = values'; ------ <'>
    
    for j = 2:n
        for i = 1:(n - j + 1)
            div_diff(i, j) = (div_diff(i + 1, j - 1) - div_diff(i, j - 1)) / (nodes(i + j - 1) - nodes(i));
        end
    end
    

    y = div_diff(1, 1);
    for i = 2:n
        term = div_diff(1, i);
        for j = 1:(i - 1)
            term = term * (x - nodes(j));
        end
        y = y + term;
    end
end



% 1.1
i = 0:2;
x_i = a + (i * (b - a)) / 3;
f_i_values = f(x_i);


% 1.2
lagrange = lagrange_poly(a, x_i, f_i_values);


% 1.3
x_j = [(5*a + b)/6, (a + b)/2, (a + 5*b)/6];
f_j_values = f(x_j);
lagrange_r = [0, 0, 0];

for iter = 0:2
    lagrange_r(iter + 1) = lagrange_poly(x_j(iter + 1), x_i, f_i_values);
end


% 1.4
figure();
plot(x, f(x), x, lagrange, x_i, f_i_values);
xlabel('x');
ylabel('y');
title('Задание 1.4');
legend('Функция f(x)', 'Многочлен Лагранжа', 'Узлы интерполяции');
grid on;
hold off;




% 1.5
x_i = ((a + b)/2) + ((b - a)/2) * cos(((2 * i + 1) * pi)/8);
f_i_values = f(x_i);
lagrange = lagrange_poly(a, x_i, f_i_values);


% 1.6 - 1.7
lagrange_c = [0, 0, 0];

for iter = 0:2
    lagrange_c(iter + 1) = lagrange_poly(x_j(iter + 1), x_i, f_i_values);
end


% 1.8
figure();
plot(x, f(x), x, lagrange, x_i, f_i_values);
xlabel('x');
ylabel('y');
title('Задание 1.8');
legend('Функция f(x)', 'Многочлен Лагранжа', 'Узлы интерполяции');
grid on;
hold off;



% 1.9
printf("\n\n_____________________________________________________________________________________________________________________\n");
printf("| j |    f(x_j)    |     L_r(x_j)    |   abs(f(x_j) - L_r(x_j))     | L_c(x_j) |      abs(f(x_j) - L_c(x_j))        |\n");

for iter = 0:2;
    
    printf("| %d |    %f |     %f   |   %f                   | %f |             %f               |\n", 
    iter, f_j_values(iter + 1), lagrange_r(iter + 1), abs(f_j_values(iter + 1) - lagrange_r(iter + 1)), 
                                                lagrange_c(iter + 1), abs(f_j_values(iter + 1) - lagrange_c(iter + 1)));

end
printf("_____________________________________________________________________________________________________________________\n");



% 2.2

i = 0:2;
x_i = a + (i * (b - a)) / 3;
f_i_values = f(x_i);
newton = newton_poly(a, x_i, f_i_values);

newton_poly_r = [0, 0, 0];

for iter = 0:2
    newton_poly_r(iter + 1) = newton_poly(x_j(iter + 1), x_i, f_i_values);
end


% 2.3
figure();
plot(x, f(x), x, newton, x_i, f_i_values);
xlabel('x');
ylabel('y');
title('Задание 2.3');
legend('Функция f(x)', 'Многочлен Лагранжа', 'Узлы интерполяции');
grid on;
hold off;



% 2.4 - 2.5
x_i = ((a + b)/2) + ((b - a)/2) * cos(((2 * i + 1) * pi)/8);
f_i_values = f(x_i);
newton = newton_poly(a, x_i, f_i_values);


newton_poly_c = [0, 0, 0];

for iter = 0:2
    newton_poly_c(iter + 1) = newton_poly(x_j(iter + 1), x_i, f_i_values);
end



% 2.6
figure();
plot(x, f(x), x, newton, x_i, f_i_values);
xlabel('x');
ylabel('y');
title('Задание 2.6');
legend('Функция f(x)', 'Многочлен Лагранжа', 'Узлы интерполяции');
grid on;
hold off;




% 2.7
printf("\n\n_____________________________________________________________________________________________________________________\n");
printf("| j |    f(x_j)    |     P_r(x_j)    |   abs(f(x_j) - P_r(x_j))     | P_c(x_j) |      abs(f(x_j) - P_c(x_j))        |\n");

for iter = 0:2;
    
    printf("| %d |    %f |     %f   |   %f                   | %f |             %f               |\n", 
    iter, f_j_values(iter + 1), newton_poly_r(iter + 1), abs(f_j_values(iter + 1) - newton_poly_r(iter + 1)), 
                                                newton_poly_c(iter + 1), abs(f_j_values(iter + 1) - newton_poly_c(iter + 1)));

end
printf("_____________________________________________________________________________________________________________________\n");