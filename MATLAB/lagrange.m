function y = lagrange(x0, y0, x)
    %Lagrange 插值
    n = length(x0); % 取 n 为 x0 数组长，即插值点个数
    m = length(x); % 取 m 为 x 数组长，即 10 * (b - a) + 1
    for i = 1:m
        z = x(i);
        s = 0.0;
        for k = 1:n
            p = 1.0;
            for j = 1:n
                if(j ~= k)
                    p = p * (z - x0(j)) / (x0(k) - x0(j)); % 求出对应拉格朗日基函数的值
                end
            end
            s = s + p * y0(k); % 将插值多项式的每一项相加，循环结束后得到完整的插值多项式
        end
        y(i) = s; % 得到插值多项式在 x 处的值
    end
end
