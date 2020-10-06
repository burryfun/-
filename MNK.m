clear all;

DEGREE = 5; % polynom degree

% INPUT POINTS
% x_p = [0.5, 1, 1.5, 2, 2.5, 3];
% y_p = [0.35, 0.8, 1.7, 1.85, 3.51, 1.02];
x_p = -1:0.1:3;
y_p = 0.5*x_p.^5-0.9*x_p.^4+0.13*x_p.^3-6*x_p.^2+0.4*x_p+5;

% Add normal distribution for initial function y_p
y_p_random = zeros(1, length(y_p));
for i = 1:length(y_p)
    r1 = random('Normal',0,1);
    y_p_random(i) = y_p(i) + r1;
end

% get interpolation function
[x,f] = getInterpFunction(x_p, y_p_random, DEGREE);

% Display all functions
p = plot(x_p,y_p,'--r', x_p, y_p_random,'ok', x, f);
p(1).LineWidth = 2;
legend('f(x)', 'f(x)+random', 'interpolation');

% S coefficient
function S = get_S(vec_x, degree)
    n = degree + 1;
    m = length(vec_x);
    S = zeros(n, n);
    for i = 0:n-1
        for j = i:n-1+i
            s = 0;
            for k = 1:m
                s = s + vec_x(k)^j;
            end
            S(i+1, j-i+1) = s;
        end
    end
end

% T coefficient
function T = get_T(vec_x, vec_y, degree)
    n = degree + 1;
    m = length(vec_x);
    T = zeros(n, 1);
    for i = 0:n-1
        t = 0;
        for k = 1:m
            t = t + vec_x(k)^i * vec_y(k);
        end
        T(i+1,1) = t;
    end
end

function [x,f] = getInterpFunction(vec_x, vec_y, degree)
    s = get_S(vec_x, degree);
    t = get_T(vec_x, vec_y, degree);
    A = linsolve(s, t);
    x = vec_x(1) : (vec_x(length(vec_x))-vec_x(1))/100 : vec_x(length(vec_x));
    f = 0;
    for i = 0:length(A)-1
        f = f + A(i+1)*x.^(i);
    end
end