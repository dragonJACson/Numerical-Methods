function o = shannonFano(p, o, s, e) % p 为概率，o 为部分编码，s 为起始元素位置，e 为终止元素位置
[nr, nc] = size(p(1:1, s:e));% 求当前总元素个数 nc
if(nc == 1)% 判断是否元素只剩下一个，如果是，直接返回
    return
end
he = ceil(nc/2); % 求当前总元素个数的一半（四舍五入）
m = s + he - 1; % m 定为当前起始位置 + 元素个数的一半 - 1
fsum = sum(p(1:1, s:e)); % 求当前所有元素的概率和
hsum = fsum/2 ; % 求概率和的一半
sumt = 0;
for n = s:m % 计算从 s 到 m 的概率和 sumt
    sumt = sumt + p(n);
    if(sumt >= hsum) % 如果在中间某元素时已经大于概率的一半了，令 m = n，退出循环
        m = n;
        break
    end
end

sume = fsum - sumt + p(m); % 概率总和减去 sumt，并加上第 m 项的概率为 sume
sume = round(sume, 2, 'significant'); % sume 四舍五入保留两位有效数字
sumt = round(sumt, 2, 'significant'); % sumt 四舍五入保留两位有效数字
if(sume < sumt) % 此时如果 sume 比 sumt 小，m 应该归于下方一组
    m = m - 1;
end

for n = s:m % 给处于上方的组每一元素加上 0
    o{n} = strcat(o{n}, '0');
end

am = m + 1;
for n  = am:e % 给处于下方的组每一元素加上 1
    o{n} = strcat(o{n}, '1');
end

o = ShannonFano(p, o, s, m); % 再次调用该函数，将下方的那组再进行分组
o = ShannonFano(p, o, am, e); % 再次调用该函数，将上方的那组再进行分组
return
end
