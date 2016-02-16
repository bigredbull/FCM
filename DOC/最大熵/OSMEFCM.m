%% 最大熵ＦＣＭ
clear
tic

% mu1=[0,0];
% sigma1=[1,0.1;0.1,1];
% mu2=[10,0];
% sigma2=[1,0.1;0.1,1];
% mu3=[0,10];
% sigma3=[1,0.1;0.1,1];
% mu4=[10,10];
% sigma4=[1,0.1;0.1,1];
% 
% data= mvnrnd(mu1,sigma1,50);
% dlabel=ones(50,1);
% data=[data;mvnrnd(mu2,sigma2,50);mvnrnd(mu3,sigma3,50);mvnrnd(mu4,sigma4,50)];
% dlabel=[dlabel;ones(50,1)*2;ones(50,1)*3;ones(50,1)*4];
% 
% gscatter(data(:,1),data(:,2),dlabel)
% axis equal

load iris.dat
data=iris(:,1:4) ;


data_n = size(data, 1); %数据多少
in_n = size(data, 2);% 数据维数


alpha=0.1;
beta=0.01;
cluster_n=3;


default_options = [2;	% exponent for the partition matrix U
		300;	% max. number of iteration
		1e-5;	% min. amount of improvement
		1];	% info display during iteration 


options = default_options;


expo = options(1);		% Exponent for U 隶属度函数的幂次方
max_iter = options(2);		% Max. iteration 最大迭代次数
min_impro = options(3);		% Min. improvement  最小进化步长
display = options(4);		% Display info or not 显示信息与否

obj_fcn = zeros(max_iter, 1);	% Array for objective function
%初始化隶属度矩阵
% U = rand(cluster_n, data_n);  
% col_sum = sum(U);
% U = U./col_sum(ones(cluster_n, 1), :);
%  U0=U;
%  save U0 U0
load U0
U=U0;

tempalpha=[];

hk=zeros(cluster_n,size(data,1));

% Main loop  data 每一行为一个数据的
for i = 1:max_iter,%迭代次数控制
	%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%
    mf = U;      % MF matrix after exponential modification
%     center = mf*data./((ones(size(data, 2), 1)*sum(mf'))'); % new center
%      center = (alpha*mf*data+(1-alpha)*hk*data)./((ones(size(data, 2), 1)*alpha*sum(mf'))'+(ones(size(data, 2), 1)*(1-alpha)*sum(hk'))') ;
      center = (alpha*mf*data+(1-alpha)*hk*data)./((ones(size(data, 2), 1)*alpha*sum(mf'))'+(ones(size(data, 2), 1)*(1-alpha)*sum(hk'))') ;
    %%%%%%%%%%%%%%%%%%
       out = zeros(size(center, 1), size(data, 1));  %每个点到每个中心的距离，行数为中心数
       % fill the output matrix
        if size(center, 2) > 1,
            for k = 1:size(center, 1),
 	             out(k, :) = sqrt(sum(((data-ones(size(data, 1), 1)*center(k, :)).^2)'));
            end
       else	% 1-D data
            for k = 1:size(center, 1),
 	             out(k, :) = abs(center(k)-data)';
            end
        end
     dist =out;
     %%%%%%%%%%%%%%%%%%%%
   for ppp=1:size(data,1)
    [ttt,sss]=min(dist(:,ppp));   
    hk(:,ppp)=0;
    hk(sss,ppp)=1;
   end
     %%%%%%%%%%%%%%%%%%%%
   
    obj_fcn(i) =sum(sum(alpha*(dist.^2).*mf))+alpha*(beta.^-1)*sum(sum(mf.*log(mf)))+sum(sum((dist.^2).*hk*(1-alpha)));  % objective function
    tmp = exp(-beta*(dist.^2));
    U= tmp./(ones(cluster_n, 1)*sum(tmp));  % 新的隶属度矩阵
    %%%%%%%%%%%%%%%%%%%%%%%%%%%

  % tempalpha=[tempalpha,alpha];
 %%%%%%%%%%%抑制式修改
%     U1=U.*alpha;
%     for j=1:data_n
%         [max_data,max_loca]=max(U(:,j));
%         tempU=alpha.*U(max_loca,j)+(1-alpha);
%         U1(max_loca,j)=tempU;
%     end
%     U=U1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %sum(tmp) 本来是一列，即某一个点到各个中心的距离之和，
 %ones(cluster_n, 1)*sum(tmp)  再扩展成n行，方便矩阵预算
 
%   Note that the situation of "singularity" (one of the data points is
%   exactly the same as one of the cluster centers) is not checked.
%   However, it hardly occurs in practice.
%EXPO: exponent (> 1) for the partition matrix.
%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%
	if display, 
		fprintf('Iteration count = %d, obj. fcn = %f\n', i, obj_fcn(i));
	end
	% check termination condition
	if i > 1,  %进化步长控制
		if abs(obj_fcn(i) - obj_fcn(i-1)) < min_impro, break; end,
	end
end

iter_n = i;	% Actual number of iterations 
obj_fcn(iter_n+1:max_iter) = [];

%绘图 坐标图像
% figure
%       plot(data(:,1), data(:,2),'o');
%       hold on;
%       maxU = max(U);
%       % Find the data points with highest grade of membership in cluster 1
%       index1 = find(U(1,:) == maxU);
%       % Find the data points with highest grade of membership in cluster 2
%       index2 = find(U(2,:) == maxU);
%       line(data(index1,1),data(index1,2),'marker','*','color','g');
%       line(data(index2,1),data(index2,2),'marker','*','color','r');
%       % Plot the cluster centers
%       plot([center([1 2],1)],[center([1 2],2)],'*','color','k')
%       title(iter_n)
%       hold off;
for i=1:data_n
    [ttt,sss]=max(U(:,i));
    for j=1:cluster_n
        U(j,i)=0;
    end
    U(sss,i)=1;
end

% [~,I]=max(U);
% 
% figure
% gscatter(data(:,1),data(:,2),I);
% axis equal    



%显示图像
% [maxU,data1]=max(U);
% image1=reshape(data1,row,column);
% figure,imagesc(image1),colormap(gray)
%maxU=max(U)
%figure;
%plot(tempalpha,'-ro');
toc