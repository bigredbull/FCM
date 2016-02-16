%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ��׼��EnFCM�㷨

a=imread('picture.bmp');
[row,column,c]=size(a);
if c~=1
    a=rgb2gray(a);
end

a=double(a);
[row,column]=size(a);
for(i=1:256)
    histogram(i)=0;
end
 
for(i=1:row)
    for(j=1:column)
        grayvalue=a(i,j)+1;
        histogram(grayvalue)=histogram(grayvalue)+1;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ��������
cluster_num=4;
threshold=0.000001;
m=1.75;
iter_num=100;
%%%%%%%%%%%%%%%%%%%%%%%
% ��ʼ��������
for(i=1:256)    
    memsum=0;
    for(k=1:cluster_num)
        membership(i,k)=rand();
        memsum=memsum+membership(i,k);
    end
    for(k=1:cluster_num)
        membership(i,k)=membership(i,k)/memsum;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ��ʼ����������

for(k=1:cluster_num)
    center(k)=0;
    memsum=0;
    for(i=1:256)
        center(k)=center(k)+membership(i,k)^m*i;
        memsum=memsum+membership(i,k)^m;
    end
    center(k)=center(k)/memsum;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ���ʼ����
for (i=1:256)
    for(k=1:cluster_num)
        dist(i,k)=abs(i-center(k));
    end
end

for(i=1:iter_num)
    costfunction(i)=0.0;
end

for(i=1:256)
    for(k=1:cluster_num)
        costfunction(1)=costfunction(1)+membership(i,k)^m* dist(i,k)^2*histogram(i);
    end    
end

fprintf('iter.count = %d, obj. fcn = %f\n',1, costfunction(1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ��ʼ����
for(it=2:iter_num)
    % ����������
    costfunction(it)=0.0;
    for(i=1:256)
        for(k=1:cluster_num)
            membership(i,k)=0;
            for(kk=1:cluster_num)
                membership(i,k)=membership(i,k)+(dist(i,k)/dist(i,kk))^(2/(m-1));
            end
            membership(i,k)=1/membership(i,k);
        end
    end     
             
    % ���¾�������
    for(k=1:cluster_num)
        center(k)=0;
        memsum=0;
        for(i=1:256)
            center(k)=center(k)+membership(i,k)^m*i;
            memsum=memsum+membership(i,k)^m;
        end        
        center(k)=center(k)/memsum;
    end
    
    % ���¾���
    for (i=1:256)
        for(k=1:cluster_num)
            dist(i,k)=abs(i-center(k));
        end        
    end
    
    % Ŀ�꺯��
    for(i=1:256)
        for(k=1:cluster_num)
            costfunction(it)=costfunction(it)+membership(i,k)^m*dist(i,k)^2;
        end
    end
   

    fprintf('iter.count = %d, obj. fcn = %f\n',it, costfunction(it));
    
    if abs(costfunction(it)-costfunction(it-1))<threshold
        break;
    end   
end

for(i=1:row)
    for(j=1:column)         
        maxmembership=membership(a(i,j)+1,1);
        index=1;
        for(k=2:cluster_num)            
            if(membership(a(i,j)+1,k)>maxmembership)
                maxmembership=membership(a(i,j)+1,k);
                index=k;
            end
        end
        newimg(i,j)=round(255*(1-(index-1)/(cluster_num-1)));
    end
end

imshow(uint8(newimg));