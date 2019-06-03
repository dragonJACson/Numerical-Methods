function [x1, k] = newton(f1, f2, f3, x0, tol)
    k = 1;
    syms l m n;
    r = f1(l, m, n);
    p = f2(l, m, n);
    q = f3(l, m, n);
    J = jacobian([r; p; q], [l m n]);
    P = double(subs(J, {l, m, n}, {x0(1), x0(2), x0(3)}));
    x1(1, 1) = f1(x0(1), x0(2), x0(3)); x1(2, 1) = f2(x0(1), x0(2), x0(3)); x1(3, 1) = f3(x0(1), x0(2), x0(3));
    x1 = x0 - P \ x1;
    while(norm(x1 - x0) >= tol)
        x0 = x1;
        P = double(subs(J, {l, m, n}, {x0(1), x0(2), x0(3)}));
        x1(1, 1) = f1(x0(1), x0(2), x0(3)); x1(2, 1) = f2(x0(1), x0(2), x0(3)); x1(3, 1) = f3(x0(1), x0(2), x0(3));
        x1 = x0 - P \ x1;
        k = k + 1;
    end
end
