function radauiia(a, b, y10, y20, d, h)
    syms X11 X12 X21 X22;
    f1 = @(t, y1, y2)2 * t * y1 * log(piecewise(y2 < 10^(-3), 10^(-3), y2 > 10^(-3), y2));
    f2 = @(t, y1, y2)(-2) * t * y2 * log(piecewise(y1 < 10^(-3), 10^(-3), y1 > 10^(-3), y1));
    realf1 = @(t, y1, y2)2 * t * y1 * log(max(y2, 10^(-3)));
    realf2 = @(t, y1, y2)(-2) * t * y2 * log(max(y1, 10^(-3)));
    df1 = @(t, y1, y2)2 * t * log(max(y2, 10^(-3)));
    df2 = @(t, y1, y2)(-2) * t * log(max(y1, 10^(-3)));
    % a 为 t 的初值，y0 为与 a 对应 y 的迭代初值，b 为区间长度，d 为阶数，h 为步长
    % a = 0, b = 10, y10 = 1, y20 = exp(1), d = 1, h = 0.01
    t0 = a; alpha = [0 -1 1]'; beta = [-1 8 5]' / 12;
    % A B C 为二级三阶 Radau IIA 方法的矩阵及向量，用于计算每次迭代的初值
    A=[5/12 -1/12; 3/4 1/4]; B = [3/4 1/4]; c = [1/3 1]';
    tc1 = t0 + c(1) * h; tc2 = t0 + c(2) * h;
    fun1 = [y10 + h * kron(A(1, 1:2), eye(d)) * [f1(tc1, X11, X21); f1(tc2, X12, X22)] - X11,...
        y20 + h * kron(A(1, 1:2), eye(d)) * [f2(tc1, X11, X21); f2(tc2, X12, X22)] - X21,...
        y10 + h * kron(A(2, 1:2), eye(d)) * [f1(tc1, X11, X21); f1(tc2, X12, X22)] - X12,...
        y20 + h * kron(A(2, 1:2), eye(d)) * [f2(tc1, X11, X21); f2(tc2, X12, X22)] - X22];
    func1 = matlabFunction(fun1, 'Vars', {[X11, X12, X21, X22]}, 'file', 'genFunc.m');
    solution = fsolve(func1, [0 0 0 0]);
    Y11 = solution(:,1); Y12 = solution(:,2); Y21 = solution(:,3); Y22 = solution(:,4);
    yy1 = y10 + h * kron(B(1:2), eye(d)) * [realf1(tc1, Y11, Y21); realf1(tc2, Y12, Y22)];
    yy2 = y20 + h * kron(B(1:2), eye(d)) * [realf2(tc1, Y11, Y21); realf2(tc2, Y12, Y22)];
    ans1 = [y10]; ans2 = [y20]; err = [0];
    for n = a + 2 * h:h:b
        t1 = t0 + h; t2 = t0 + 2 * h; tc1 = t0 + (c(1) + 1) * h;
        tc2 = t0 + (c(2) + 1) * h;
        Y11 = yy1; Y21 = yy2;
        fun2 = [yy1 + h * kron(A(1, 1:2), eye(d)) * [f1(tc1, X11, X21); f1(tc2, X12, X22)] - X11,...
            yy2 + h * kron(A(1, 1:2), eye(d)) * [f2(tc1, X11, X21); f2(tc2, X12, X22)] - X21,...
            yy1 + h * kron(A(2, 1:2), eye(d)) * [f1(tc1, X11, X21); f1(tc2, X12, X22)] - X12,...
            yy2 + h * kron(A(2, 1:2), eye(d)) * [f2(tc1, X11, X21); f2(tc2, X12, X22)] - X22];
        func2 = matlabFunction(fun2, 'Vars', {[X11, X12, X21, X22]}, 'file', 'genFunc2.m');
        solution = fsolve(func2, [Y11 Y12 Y21 Y22]);
        Y11 = solution(:,1); Y12 = solution(:,2); Y21 = solution(:,3); Y22 = solution(:,4);
        yy120 = yy1 + h * kron(B(1:2), eye(d)) * [realf1(tc1, Y11, Y21); realf1(tc2, Y12, Y22)];
        yy220 = yy2 + h * kron(B(1:2), eye(d)) * [realf2(tc1, Y11, Y21); realf2(tc2, Y12, Y22)];
        w1 = h * (beta(1) * realf1(t0, y10, y20) + beta(2) * realf1(t1, yy1, yy2)) - (alpha(1) * y10 + alpha(2) * yy1);
        w2 = h * (beta(1) * realf2(t0, y10, y20) + beta(2) * realf2(t1, yy1, yy2)) - (alpha(1) * y20 + alpha(2) * yy2);
        err11 = 1; err12 = 1; err21 = 1; err22 = 1;
        while err11 >= 10^(-12) && err12 >= 10^(-12) && err21 >= 10^(-12) && err22 >= 10^(-12)
            r1 = yy120 - h * beta(3) * realf1(t2, yy120, yy220) - w1;
            r2 = yy220 - h * beta(3) * realf2(t2, yy120, yy220) - w2;
            yy121 = yy120 - (eye(d) - h * beta(3) * df1(t2, yy120, yy220)) \ r1;
            yy221 = yy220 - (eye(d) - h * beta(3) * df2(t2, yy120, yy220)) \ r2;
            err11 = norm(yy121 - yy120); err12 = norm(r1);
            err21 = norm(yy221 - yy220); err22 = norm(r2);
            yy120 = yy121;
            yy220 = yy221;
        end
        t0 = t0 + h;
        yy12 = yy120; y10 = yy1; yy1 = yy12;
        yy22 = yy220; y20 = yy2; yy2 = yy22;
        err = [err; norm([yy12; yy22] - [exp(sin(n^2)); exp(cos(n^2))])];
        ans1 = [ans1; yy12];
        ans2 = [ans2; yy22];
    end
    plot(a+h:h:b, ans1);
    hold on;
    plot(a+h:h:b, ans2);
    figure;
    plot(a+h:h:b, err);
end
