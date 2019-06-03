function [h, e] = huffmanCode(p)
    % 此处显示详细说明
    % p 为概率分布，此函数功能是进行哈夫曼编码
    % h 为各个元素的麻子
    % e 为输出的平均码长

    if length(find(p<0))~=0
        error('概率不应该小于 0！')
    end

    if abs(sum(p)-1)>10e-10
        error('概率之和大于 1，请检查输入！')
    end

    n=length(p);

    p=sort(p)
    q=p;
    m=zeros(n-1,n);
    for i=1:n-1
        [q,e]=sort(q);
        m(i,:)=[e(1:n-i+1),zeros(1,i-1)]; % 由数组 l 构建一个矩阵，该矩阵表明概率合并时的顺序，用于后面的编码
        q=[q(1)+q(2),q(3:n),1];
    end

    for i=1:n-1
        c(i,1:n*n)=blanks(n*n); % c 矩阵用于进行 huffman 编码
    end
        c(n-1,n)='1'; % 由于 a 矩阵的第 n-1 行的前两个元素为进行 huffman 编码加和运算时所得的最后两个概率（在本例中为 0.02、0.08），因此其值为 0 或 1
        c(n-1,2*n)='0';
    for i=2:n-1
        c(n-i,1:n-1)=c(n-i+1,n*(find(m(n-i+1,:)==1))-(n-2):n*(find(m(n-i+1,:)==1))); % 矩阵 c 的第 n-i 的第一个元素的 n-1 的字符赋值为对应于 a 矩阵中第 n-i+1 行中值为 1 的位置在 c 矩阵中的编码值
        c(n-i,n)='0';
        c(n-i,n+1:2*n-1)=c(n-i,1:n-1); % 矩阵 c 的第 n-i 的第二个元素的 n-1 的字符与第 n-i 行的第一个元素的前 n-1 个符号相同，因为其根节点相同
        c(n-i,2*n)='1';
        for j=1:i-1
            c(n-i,(j+1)*n+1:(j+2)*n)=c(n-i+1,n*(find(m(n-i+1,:)==j+1)-1)+1:n*find(m(n-i+1,:)==j+1));
                % 矩阵 c 中第 n-i 行第 j+1 列的值等于对应于 a 矩阵中第 n-i+1 行中值为 j+1 的前面一个元素的位置在 c 矩阵中的编码值
        end
    end
    for i=1:n
        h(i,1:n)=c(1,n*(find(m(1,:)==i)-1)+1:find(m(1,:)==i)*n); % 用 h 表示最后的 huffman 编码
        len(i)=length(find(abs(h(i,:))~=32)); % 计算每一个编码的长度
    end
    e=sum(p.*len); % 计算平均码长
end
