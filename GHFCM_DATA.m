function [H, center,  iter] = GHFCM_DATA(data, c, m, maxIter, epsilon)

    % 数据的样本数 n 和特征维度 d
    [n, d] = size(data);
    % % 初始化聚类中心(方法一)
    % center = zeros(c, d);
    % for j = 1:c
    %     l1 = round(size(data, 1) / c * (j - 1) + size(data, 1) / (2 * c));
    %     center(j, :) = data(l1, :);  % 四舍五入索引并赋值
    % end
    % 初始化聚类中心（方法二）
    center = zeros(c, d); % 初始化聚类中心
    index = randperm(n);   % 对样本序数随机排列
    center = data(index(1:c),:);  % 选取随机排列的序数的前cluster_n个
    %初始化隶属度矩阵
    % lambda = 1;  % 指数分布的参数（控制分布的形状）
    % H = 1 + exprnd(lambda, c, n);  % 生成c x n的矩阵，范围为[1, +∞)
    % H = H ./ sum(1./H); % 确保满足隶属度约束条件

    
%%
    %开始迭代
     for iter = 1:maxIter

        %计算data_i和center的算欧式距离
        distance = zeros(size(center,1),size(data,1));
        for k=1:size(center,1)
            for i=1:size(data,1)
                 distance(k,i) = (norm(data(i, :) - center(k, :)))^2 + 0.0001;
            end
        end

        %更新隶属度矩阵
        H = zeros(size(center,1),size(data,1));
        for k=1:size(center,1)
            for i=1:size(data,1)
                %隶属度_分母
                denominator_H = distance(k,i)^(1/(m+1));
                %隶属度_分子
                numerator_H = 0;
                for r=1:1:size(center,1)
                  numerator_H = numerator_H + distance(r,i)^(1/(m+1));
                end
                  H(k,i) = numerator_H/denominator_H;
            end
        end

        %更新聚类中心
        center_new = zeros(c, d);
        for k=1:size(center,1)
            num = 0;
            den = 0;
            for i =1:n
            num = num +( H(k,i)^m) .* data(i,:);
            den = den + H(k,i)^m;
            end
            center_new(k,:) = num/den;
        end

        %收敛性判断
        temp = 0.0;
        for k =1:c
            temp =temp + (center_new(k)-center(k))^2;
        end
        if temp < epsilon
            break;
        end
         center = center_new;
     end
end