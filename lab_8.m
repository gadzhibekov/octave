a = 0.1;
b = 0.8;
l = b - a;
A = 3;
B = 1;
tau = 0.05; 
h = 0.1;
T = 1;

x = a:h:b;
Nx = length(x);
nt = round(T/tau);
tvec = 0:tau:T;

u = zeros(Nx, nt+1);
for i = 1:Nx
    u(i,1) = ((A - B) * (x(i) - a)) / l + A;
end

k = @(x) cos(x);
f = @(x) 10 * sin(x);

k_half = zeros(Nx-1, 1);
for i = 1:(Nx-1)
    k_half(i) = (k(x(i)) + k(x(i+1))) / 2;
end

for n = 1:nt
    u_new = u(:, n);
    current_t = tvec(n);
    
    for i = 2:(Nx-1)
        d2udx = ( k_half(i) * (u(i+1, n) - u(i, n)) - k_half(i-1) * (u(i, n) - u(i-1, n)) ) / h^2;
        u_new(i) = u(i, n) + tau * ( d2udx + f(x(i)) * (1 - exp(-current_t)));
    end
    
    u_new(1) = A;
    u_new(Nx) = B;
    
    u(:, n+1) = u_new;
end

u_halfTau = (u(:,1) + u(:,2)) / 2;

index_10tau = 10 + 1;
u_10tau = u(:, index_10tau);

index_20tau = 20 + 1;
u_20tau = u(:, index_20tau);

figure;
plot(x, u_halfTau, 'r-o', 'LineWidth', 1.5); hold on;
plot(x, u_10tau, 'b-s', 'LineWidth', 1.5);
plot(x, u_20tau, 'g-^', 'LineWidth', 1.5);
xlabel('x');
ylabel('u(x,t)');
legend('t = 0.5\tau', 't = 10\tau', 't = 20\tau');
title('Приближенное решение задачи методом явной разностной схемы');
grid on;
