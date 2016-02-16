clear all
clc;
%% 
load wine.txt
data = wine;
data(:,1) = [];

data_n = size(data,1);
M = size(data,2);

cluster_n = 3;
minr = 0.005;%在计算高斯核时的参数

default_options = [1.08;	% 隶属度函数的幂次方
		300;	% 最大迭代次数
		0.0001;	%步长
		1];	% 判定条件
    
options = default_options;

expo = options(1);		% Exponent for U 隶属度函数的幂次方
max_iter = options(2);		% Max. iteration 最大迭代次数
min_impro = options(3);		% Min. improvement  最小进化步长
display = options(4);

Uv = zeros(data_n,cluster_n);
membership = zeros(data_n,cluster_n);
obj_fcn  = zeros(max_iter,1);
K = zeros(M,data_n,data_n);
A = zeros(M,data_n,cluster_n);
B = zeros(M,1);
W = zeros(M,1);
D = zeros(data_n,cluster_n);

%%  初始化隶属度矩阵并归一
U = rand(data_n,cluster_n); 

S = sum(U,2);
for i = 1:data_n
    for j = 1:cluster_n
        U(i,j) = U(i,j) / S(i,1);
    end
end

%% 计算核函数
opeo = zeros(M,1);
for i = 1:M
    oper = zeros(data_n,data_n);
    for j = 1:data_n
        for p = 1:data_n
        oper(j,p) = (-1 * ((data(j,i) - data(p,i)).^2) / log(minr));
        end 
    end
    oprt = reshape(oper,1,data_n * data_n);
    oprt(oprt==0) = inf; 
    opeo(i,1) = min(oprt);
end

for i = 1:M
    for p = 1:data_n
        for q = 1:data_n
        K(i,p,q) = exp((-1 *(data(p,i) - data(q,i)).^2)/ opeo(i,1));
        end
    end
end

%% 
for o = 1:max_iter    
    tic    
    %%%%%标准化隶属度矩阵
    sum = 0;
    for i = 1:data_n
        for j = 1:cluster_n
            sum = sum + U(i,j).^expo;
        end
    end
    for i = 1:data_n
        for j = 1:cluster_n
            Uv(i,j) = U(i,j).^expo / sum;
        end
    end
    
    %%%%%%%求参数a  
    tmp = zeros(M,data_n,cluster_n);
    for i = 1:data_n
        for c = 1:cluster_n
            for k =1:M
                for j = 1:data_n
                    tmp(k,i,c) = tmp(k,i,c) + Uv(j,c) * K(k,i,j);
                end
            end
        end
    end
    
    for i = 1:data_n
        for c = 1:cluster_n
            for k =1:M
                tp = 0;
                for j = 1:data_n
                    for q = data_n
                        tp = tp + Uv(j,c)* Uv(q,c) * K(k,j,q);
                    end
                end
                A(k,i,c) = K(k,i,i) - 2 * tmp(k,i,c) + tp;
            end
        end
    end
    
    %%%%%求系数B
    for k = 1:M
        for i = data_n
            for j = 1:cluster_n
               B(k,1) = B(k,1) + U(i,j).^expo * A(k,i,j);
            end
        end
    end
    
    %%%%%%求权重w
    sum = 0;
    for k = 1:M
        sum = sum + B(k,1).^(-1);
    end
    for k = 1:M
        W(k,1) = B(k,1).^(-1) / sum;
    end
    
    %%%%%%求欧式距离
    for i = 1:data_n
        for c = 1:cluster_n
            for k = 1:M
                D(i,c) = D(i,c) + A(k,i,c) * W(k,1).^2;
            end
        end
    end
    
    %%%%%%更新隶属度矩阵
    for i = 1:data_n
        for c = 1:cluster_n
            sum = 0;
            for q = 1:cluster_n
                sum = sum + (D(i,c) / D(i,q)).^(1 / (expo - 1));
            end
            membership(i,c) = 1 / sum;
        end
    end
    
    sum = zeros(data_n,1);
    for i =1:data_n
        for j = 1:cluster_n
            sum(i,1) = sum(i,1) + membership(i,j);
        end
    end
    
    for i = 1:data
        for j = 1:cluster_n
            membership(i,j) = membership(i,j) / sum(i,1);
        end
    end    

    %%%%%%得到目标函数
    temp = 0;
    for i = 1:data_n
        for c = 1:cluster_n
            temp = temp + membership(i,c).^expo;
            trmp = 0;
            for k = 1:M
                trmp = trmp + A(k,i,c) * W(k,1).^2 ;
            end
        end
    end
    obj_fcn(o) = temp * trmp;
    
    if display,
        fprintf('Iteration count = %d, obj_fcn = %f\n', o, obj_fcn(o));
%         fprintf('Iteration count = %d\n', o);
    end
    
    if o > 1
        if max(abs(membership - U)) < min_impro;
            break
        else
            U  = membership;
        end
    end
    toc
   
end

B = ones(data_n,1);
for i = 1:data_n
    for j = 1:cluster_n
        if membership(i,B(i,1)) < membership(i,j)
            B(i,1) = j;
        end
    end
end
B = reshape(B,1,data_n);
A = wine(:,1);
A = reshape(A,1,data_n);
t = nmi(A,B);

plot(membership(:,1),'-ro');
grid on
hold on
plot(membership(:,2),'-g*');
plot(membership(:,3),'-b+');

ylabel('Membership degrees')
legend('MKFC1','MKFC2','MKFC3','location','northeast');

