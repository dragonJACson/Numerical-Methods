function useLagrange(fun, n, a, b)
    % Lagrange 插值调用函数
    % fun 为试验函数
    % n 为插值多项式的次数
    % [a, b] 为插值区间
    x0 = linspace(a, b, n + 1);y0 = feval(fun, x0); % 在插值区间上均匀地取 n + 1 个点，得到对应的节点和函数值
    x = a:0.1:b;y = Lagrange(x0, y0, x); % 以 0.1 为步长取 [a, b] 上的点，求出对应的插值多项式的值
    fplot(fun, [a b], 'r-'); % 画出试验函数在区间 [a, b] 上的图像，红色实线
    hold on;
    plot(x, y, 'b--'); % 画出插值函数图像，蓝色虚线
    xlabel('x');ylabel('y = f(x) o and y= Ln(x) --');
end
