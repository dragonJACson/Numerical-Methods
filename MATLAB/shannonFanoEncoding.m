function shannonFanoEncoding(p) % 范诺编码主程序
[~, nc] = size(p); % nc 是概率的个数
x = 1:1:nc; % x 为 1 到 nc

o = {' ' , ' '}; % 创建一个空向量用于储存输出结果
for n = 1:nc
    o{n} = ' '; % 由概率个数决定向量维数
    p(n) = round(p(n), 4, 'significant'); % 将 p(n) 中的概率四舍五入到 4 位有效数字
end

[p, x] = bubbleSort(p, x); % 将 p 中的概率进行冒泡排序（从大到小）

o = ShannonFano(p, o, 1, nc); % 运行范诺编码算法程序得到结果

entropy = 0;
averagel = 0;
efficiency = 0;
for n = 1:nc
    l(n) = length(cell2mat(o(n))); % 得到每个消息符号对应的码长
    averagel = averagel + l(n) * p(n); % 用公式计算平均码长
    entropy = entropy - p(n) * log2(p(n)); % 用公式计算熵
end
efficiency = entropy / averagel;

for n = 1:nc
    output = [x(n), p(n), o(n), l(n)]; % 按概率从大到小的顺序，输出编号，概率，编码
    disp(output);
% disp(x);
% disp(o);
end

fprintf('Entropy: %.4f (bit) \n', entropy);
fprintf('Average length: %.4f \n', averagel);
fprintf('Efficiency: %2.2f%% \n', efficiency * 100)
