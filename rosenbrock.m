function f = rosenbrock(x,a)
f = a.a*(x(2) - x(1)^2)^2 + (a.b - x(1))^2;