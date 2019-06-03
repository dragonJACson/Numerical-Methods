function [k, xnew] = secant(x0, x1, f, tol)
    xnew = x1;
    k = 0;
    while(abs(x0 - xnew) >= tol || abs(f(xnew)) >= tol)
        xnew = x1 - ((x1 - x0) / (f(x1) - f(x0))) * f(x1);
        x0 = x1;
        x1 = xnew;
        k = k + 1;
    end
end
