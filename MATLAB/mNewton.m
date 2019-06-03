function x0 = mNewton(f, x0, tol)
    syms x;
    f1 = matlabFunction(diff(f, x));
    f0 = matlabFunction(f);
    x1 = x0 - f0(x0) / f1(x0);
    while abs(x1 - x0) >= tol || abs(f0(x0)) >= tol
        x1 = x0 - f0(x0) / f1(x0);
        x0 = x1;
        x0;
    end
    fprintf('The result is %f \n', x0);
end

