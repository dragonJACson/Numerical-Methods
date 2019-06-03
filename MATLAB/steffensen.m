function [k, xnew] = steffensen(x0, f, tol)
    k = 1;
    xnew = x0 - f(x0) * f(x0) / (f(x0) - f(x0 - f(x0)));
    while(abs(xnew - x0) >= tol || abs(f(xnew)) >= tol)
        xnew = x0 - f(x0) * f(x0) / (f(x0) - f(x0 - f(x0)));
        x0 = xnew;
        k = k + 1;
    end
end
