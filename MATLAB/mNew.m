function x0 = mNew(x0, tol)
    f = inline('x + sin(x) - 1', 'x');
    x1 = (cos(1 / 2) * x0 - sin(x0) + 1) / (1 + cos(1 / 2));
    while abs(x1 - x0) >= tol || abs(f(x0)) >= tol
        x1 = (cos(1 / 2) * x0 - sin(x0) + 1) / (1 + cos(1 / 2));
        x0 = x1;
        x0;
    end
end

