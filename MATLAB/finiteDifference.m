f = @(x, y)(2 * pi^2 * exp(1)^(pi * (x + y)) * sin(pi * x + pi * y));
% U(i, j) 对应列数 64 * i + j - 64
% U(i, j) = xxx 对应列数 64 * i - 63
A = zeros(3969, 3969);
for i = 1:63
    for j = 1:63
        A(63 * i + j - 63, 63 * i + j - 63) = 4;
        if(i - 1 ~= 0)
            A(63 * i + j - 63, 63 * i + j - 126) = -1;
        end
        if(j - 1 ~= 0)
            A(63 * i + j - 63, 63 * i + j - 64) = -1;
        end
        if(i + 1 ~= 64)
            A(63 * i + j - 63, 63 * i + j) = -1;
        end
        if(j + 1 ~= 64)
            A(63 * i + j - 63, 63 * i + j - 62) = -1;
        end
        b(63 * i + j - 63) = (1 / 64)^2 / 4 * f(i / 64, j / 64);
        b = b';
    end
end



