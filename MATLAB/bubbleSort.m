function [p, x] = bubbleSort(p, x) % 冒泡算法程序，p 为所有概率，x 为按概率输入顺序对概率编上的号码
    [~, nc] = size(p); % nc 为概率总数
    for i = 1:nc-1 % i 从 1 到 nc - 1
        for j = i+1:nc % j 从 i + 1 到 nc
            if(p(i) < p(j)) % 如果第 i 项元素比第 j 项元素小，将两项交换位置
                temp = p(i);
                p(i) = p(j);
                p(j) = temp;
                tempx = x(i); % 对应编号也要交换
                x(i) = x(j);
                x(j) = tempx;
            end
        end
    end
    return
end
