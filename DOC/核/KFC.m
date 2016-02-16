clear all
clc;

load glass.txt;
cluster_n=6;

data = glass;
data(:,1) = [];
data(:,size(data,2)) = [];
data_n = size(data, 1); %数据的个数    size(A,1)返回A的行数，size(A,2)返回A的列数
in_n = size(data, 2);% 数据维数

default_options = [10;	% 隶属度函数的幂次方
		300;	% 最大迭代次数
		1e-5;	%步长
		1];	% 判定条件

options = default_options;


expo = options(1);		% Exponent for U 隶属度函数的幂次方
max_iter = options(2);		% Max. iteration 最大迭代次数
min_impro = options(3);		% Min. improvement  最小进化步长
display = options(4);		% Display info or not 显示信息与否


obj_fcn = zeros(max_iter, 1);	% Array for objective function
%初始化隶属度矩阵
U = rand(cluster_n, data_n); %rand()产生随机矩阵 
col_sum = sum(U);
U = U./col_sum(ones(cluster_n, 1), :);%归一化

for i = 1:max_iter,%迭代次数控制
	
    tic
    mf = U.^expo;      % MF matrix after exponential modification
    center = mf*data./((ones(size(data, 2), 1)*sum((mf')))'); % 建立新的聚类中心
   
       out = zeros(size(center, 1), size(data, 1));  %每个点到每个中心的距离，行数为中心数
     
        if size(center, 2) > 1,%样本的维数大于一执行以下程序
            for k = 1:size(center, 1),%给K赋值
 	             out(k, :) = sqrt(sum((((data-ones(size(data, 1), 1)*center(k, :)).^2)')));%得到欧氏距离
            end
       else	% 1-D data
            for k = 1:size(center, 1),
 	             out(k, :) = abs(center(k)-data)';
            end
        end
    % dist =out;
   %%%%%%%%%%%%%%%%%%%%
    obj_fcn(i) = sum(sum((out.^2).*mf));  % 目标函数
    tmp = out.^(-2/(expo-1));      % 根据新的隶属度矩阵求解公式求出
    U= tmp./(ones(cluster_n, 1)*sum(tmp));  % 新的隶属度矩阵

	if display, 
		fprintf('Iteration count = %d, obj. fcn = %f\n', i, obj_fcn(i));
        %输出迭代次数和函数的结果
	end
	% check termination condition
	if i > 1,  %进化步长控制
		if abs(obj_fcn(i) - obj_fcn(i-1)) < min_impro, break; end,
    end
    toc
end

iter_n = i;	% 实际的迭代次数 
obj_fcn(iter_n+1:max_iter) = [];

B = ones(data_n,1);
for i = 1:size(U,1)
    for j = 1:size(U,2)
        if U(B(j,1),j)< U(i,j)
            B(j,1) = i;
        end
    end
end

t = size(glass,2);
A = glass(:,t);

B = reshape(B,1,data_n);
A = reshape(A,1,data_n);
t = nmi(A,B);

plot(U(1,:),'-ro');
grid on
hold on
plot(U(2,:),'-g*');
plot(U(3,:),'-b+');
plot(U(4,:),'-m*');
plot(U(5,:),'-y+');
plot(U(6,:),'-k*');


ylabel('Membership degrees')
legend('FCM1','FCM2','FCM3','FCM4','FCM5','FCM6','location','northeast');

toc