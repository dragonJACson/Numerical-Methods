function [x1, n] = picard(f1, f2, f3, x0, tol)
    n = 1;
    x1(1) = f1(x0(2), x0(3)); x1(2) = f2(x0(1), x0(3)); x1(3) = f3(x0(1), x0(2));
    while(norm(x1 - x0) >= tol)
        x0 = x1;
        x1(1) = f1(x0(2), x0(3)); x1(2) = f2(x0(1), x0(3)); x1(3) = f3(x0(1), x0(2));
        n = n + 1;
    end
end
