function [H, center,  iter] = GHFKHM_DATA(data, c, m, maxIter, epsilon)

    % 数据的样本数 n 和特征维度 d
    [n, d] = size(data);
    % 初始化聚类中心(方法一)
    center = zeros(c, d);
    for j = 1:c
        l1 = round(size(data, 1) / c * (j - 1) + size(data, 1) / (2 * c));
        center(j, :) = data(l1, :);  % 四舍五入索引并赋值
    end
    % % 初始化聚类中心（方法二）
    % center = zeros(c, d); % 初始化聚类中心
    % index = randperm(n);   % 对样本序数随机排列
    % center = data(index(1:c),:);  % 选取随机排列的序数的前cluster_n个
    %初始化隶属度矩阵
    % lambda = 1;  % 指数分布的参数（控制分布的形状）
    % H = 1 + exprnd(lambda, c, n);  % 生成c x n的矩阵，范围为[1, +∞)
    % H = H ./ sum(1./H); % 确保满足隶属度约束条件

    %%
    %开始迭代
    for iter = 1:maxIter

        % 计算data_i和center的算欧式距离
        distance = zeros(size(center,1),size(data,1));
        for k=1:size(center,1)
            for i=1:size(data,1)
                distance(k,i) = (norm(data(i, :) - center(k, :)))^2 + 0.0001;
            end
        end

        % 更新隶属度矩阵
        H = zeros(size(center,1),size(data,1));
        for k=1:size(center,1)
            for i=1:size(data,1)
                % 分母部分
                denominator_H = (distance(k,i))^(1/(m-1));
                numerator_H = 0;
                % 分子部分
                for r=1:1:size(center,1)
                  numerator_H = numerator_H + (distance(r,i))^(1/(m-1));
                end
                % 计算H
                H(k,i) = numerator_H/denominator_H;
            end
        end
        % % 更新聚类中心
        % center_new = zeros(c, d);
        % for k=1:size(center,1)
        %     % 分母
        %     denominator = 0;
        %     % 分子
        %     numerator = 0;
        %     % 遍历每个数据点 i
        %     for i = 1:n
        %         % 计算分母和分子所需的临时变量 tep1 和 tep2
        %         tep1 = 0; % 分母的加和项
        %         tep2 = 0; % 分子的加和项
        % 
        %         % 遍历每个聚类中心 r
        %         for r = 1:size(center, 1)
        %             tep1 = tep1 + (H(k, i)^m * distance(k, i) / (H(r, i)^m * distance(r, i)));
        %         end
        % 
        %         % 累加到分母
        %         denominator = denominator + H(k, i)^m / tep1^2;
        % 
        %         % 计算分子部分（加权数据点）
        %         for r = 1:size(center, 1)
        %             tep2 = tep2 + (H(k, i)^m * distance(k, i) / (H(r, i)^m * distance(r, i)));
        %         end
        %         numerator = numerator + (H(k, i)^m / tep2^2) .* data(i, :);
        %     end
        % 
        %     center_new(k,:) = numerator/denominator;
        % end



       % 更新聚类中心
        center_new = zeros(c, d);
        for k=1:size(center,1)

            % 分母
            denominator = 0;
            for i =1:n
                tep1 = 0.0;
                for  r=1:size(center,1)
                    tep1 = tep1 +( H(k,i)^m * distance(k,i) /( H(r,i)^m * distance(r,i) ) );
                end
                denominator = denominator + H(k,i)^m/tep1^2;
            end

            % 分子
            numerator = 0;
            for i =1:n
                tep2 = 0.0;
                for  r=1:size(center,1)
                    tep2 = tep2 +( H(k,i)^m * distance(k,i)/( H(r,i)^m * distance(r,i)) );
                end
                numerator = numerator + H(k,i)^m/tep2^2.* data(i,:);
            end
            center_new(k,:) = numerator/denominator;
        end

        % 收敛性判断
        temp = 0.0;
        for k =1:c
            temp =temp + norm(center_new(k,:)-center(k,:))^2;
        end
        % disp(temp)
        if temp < epsilon
            break;
        end
        center = center_new;
    end
end