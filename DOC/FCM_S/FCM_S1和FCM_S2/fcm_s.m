%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 标准的FCM_S算法

a=imread('picture.bmp');
[height,width,c]=size(a);
if c~=1
    a=rgb2gray(a);
end
a=double(a);
[row,column]=size(a);  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 几个参数
cluster_num=4;
threshold=0.00001;
m=1.75;
iter_num=100;
alpha=2.0;
%%%%%%%%%%%%%%%%%%%%%%%
% 初始化隶属度

for(i=1:height)
    for(j=1:width)
        member_sum=0;
        for(k=1:cluster_num)
            membership(i,j,k)=rand();
            member_sum=member_sum+membership(i,j,k);
        end        
        for(k=1:cluster_num)
            membership(i,j,k)=membership(i,j,k)/member_sum;
        end
    end
end

%计算初始中心
for(k=1:cluster_num)
    sum_down=0;
    sum_up=0;
    for(i=1:height)
        for(j=1:width)
            sum_down=sum_down+(1+alpha)*membership(i,j,k)^m;  
            num_neigh=0;
            sum_gray=0;
            up=i-1;down=i+1;left=j-1;right=j+1;
            if(up<1) up=1; end
            if (down>height) down=height; end
            if (left<1) left=1; end
            if (right>width) right=width; end
            
            for(p=up:down)
                for(q=left:right)
                    sum_gray=sum_gray+a(p,q);
                    num_neigh=num_neigh+1;
                end
            end
            num_neigh=num_neigh-1;
            sum_gray=sum_gray-a(i,j);
            sum_up=sum_up+membership(i,j,k)^m*(a(i,j)+alpha*sum_gray/num_neigh);              
        end
    end
    center(k)=sum_up/sum_down;
end

for(it=1:iter_num)
    costfunction(it)=0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%计算目标函数

for(i=1:height)
    for(j=1:width)
        for(k=1:cluster_num)
            costfunction(1)=costfunction(1)+membership(i,j,k)^m*(a(i,j)-center(k))^2;
            % 加上邻域项
            num_neigh=0;
            temp_sum=0;
            up=i-1;down=i+1;left=j-1;right=j+1;
            if(up<1) up=1; end
            if (down>height) down=height; end
            if (left<1) left=1; end
            if (right>width) right=width; end
            
            for(p=up:down)
                for(q=left:right)
                    temp_sum=temp_sum+(a(p,q)-center(k))^2;
                    num_neigh=num_neigh+1;
                end
            end
            temp_sum=temp_sum-(a(i,j)-center(k))^2;
            num_neigh=num_neigh-1;
            costfunction(1)=costfunction(1)+membership(i,j,k)^m*alpha*temp_sum/num_neigh;            
        end
    end
end
fprintf('iter.count = %d, obj. fcn = %f\n',1, costfunction(1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 开始迭代

for(it=2:iter_num)
    %计算隶属度
    for(i=1:height)
        for(j=1:width)
            sum_down=0;       
            
            for(k=1:cluster_num)%求分母，它是对所有类求
                num_neigh=0;
                temp_sum=0;
                up=i-1;down=i+1;left=j-1;right=j+1;
                if(up<1) up=1; end
                if (down>height) down=height; end
                if (left<1) left=1; end
                if (right>width) right=width; end
                
                for(p=up:down)
                    for(q=left:right)
                        temp_sum=temp_sum+(a(p,q)-center(k))^2;
                        num_neigh=num_neigh+1;
                    end
                end
                temp_sum=temp_sum-(a(i,j)-center(k))^2;
                num_neigh=num_neigh-1;
                
                sum_down=sum_down+((a(i,j)-center(k))^2+alpha*temp_sum/num_neigh)^(-1/(m-1));                         
            end
            
            for(k=1:cluster_num)%membership(i,j,k)
                num_neigh=0;
                temp_sum=0;
                up=i-1;down=i+1;left=j-1;right=j+1;
                if(up<1) up=1; end
                if (down>height) down=height; end
                if (left<1) left=1; end
                if (right>width) right=width; end
                
                for(p=up:down)
                    for(q=left:right)
                        temp_sum=temp_sum+(a(p,q)-center(k))^2;
                        num_neigh=num_neigh+1;
                    end
                end
                temp_sum=temp_sum-(a(i,j)-center(k))^2;
                num_neigh=num_neigh-1;
                
                sum_up=(alpha*temp_sum/num_neigh+(a(i,j)-center(k))^2)^(-1/(m-1));
                membership(i,j,k)=sum_up/sum_down;
            end           
        end        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %求聚类中心
    for(k=1:cluster_num)
        sum_down=0;sum_up=0;
        for(i=1:height)
            for(j=1:width)
                sum_down=sum_down+(1+alpha)*membership(i,j,k)^m;            
                num_neigh=0;
                sum_gray=0;
                up=i-1;down=i+1;left=j-1;right=j+1;
                if(up<1) up=1; end
                if (down>height) down=height; end
                if (left<1) left=1; end
                if (right>width) right=width; end
                
                for(p=up:down)
                    for(q=left:right)
                        sum_gray=sum_gray+a(p,q);
                        num_neigh=num_neigh+1;
                    end
                end
                num_neigh=num_neigh-1;
                sum_gray=sum_gray-a(i,j);
                sum_up=sum_up+membership(i,j,k)^m*(a(i,j)+alpha*sum_gray/num_neigh);  
            end
        end
        center(k)=sum_up/sum_down;
    end      
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %求目标函数
    
    for(i=1:height)
        for(j=1:width)
            for(k=1:cluster_num)
                costfunction(it)=costfunction(it)+membership(i,j,k)^m*(a(i,j)-center(k))^2;
                % 加上邻域项
                num_neigh=0;
                temp_sum=0;
                up=i-1;down=i+1;left=j-1;right=j+1;
                if(up<1) up=1; end
                if (down>height) down=height; end
                if (left<1) left=1; end
                if (right>width) right=width; end
                
                for(p=up:down)
                    for(q=left:right)
                        temp_sum=temp_sum+(a(p,q)-center(k))^2;
                        num_neigh=num_neigh+1;
                    end
                end
                temp_sum=temp_sum-(a(i,j)-center(k))^2;
                num_neigh=num_neigh-1;
                costfunction(it)=costfunction(it)+membership(i,j,k)^m*alpha*temp_sum/num_neigh;            
            end
        end
    end
    fprintf('iter.count = %d, obj. fcn = %f\n',it, costfunction(it));
    
    if abs(costfunction(it)-costfunction(it-1))<threshold
        break;
    end   
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%输出图像
for(i=1:row)
    for(j=1:column)         
        maxmembership=membership(i,j,1);
        index=1;
        for(k=2:cluster_num)            
            if(membership(i,j,k)>maxmembership)
                maxmembership=membership(i,j,k);
                index=k;
            end
        end
        newimg(i,j)=round(255*(1-(index-1)/(cluster_num-1)));
    end
end

imshow(uint8(newimg));