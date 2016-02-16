function [C]=calculateC(I1,w)
[m,n]=size(I1);
C=cell(m,n);
for i=1:m
    for j=1:n
        C{i,j}=zeros();
    end
end

for i=1:m
    for j=1:n
%         %%%%计算中心点的方差均值比%%%%
        up=i-1;down=i+1;left=j-1;right=j+1;
            if (up<1) up=1; end
            if (down>m) down=m; end
            if (left<1) left=1; end
            if (right>n) right=n; end
        q=1;A=zeros();%sum=0;
        for x=up:down
            for y=left:right
                A(q)=I1(x,y);
%                 sum=sum+(I1(x,y)-I1(i,j))^2;
                q=q+1;
            end
        end
%         C_var=sqrt(sum/(q-2));        
        C_means=mean(A);
        C_var=var(A);
        if (C_means~=0)
            C_center=C_var/C_means^2;
        else
            C_center=0;
        end
        %%%%计算邻域点的方差均值比%%%%
        k=1;C_sum=0;C1=zeros();
        for x=up:down
            for y=left:right
                if (x~=i||y~=j)
                   up1=x-1;down1=x+1;left1=y-1;right1=y+1;%邻域点的3x3的小窗口
                   if(up1<1) up1=1; end
                   if (down1>m) down1=m; end
                   if (left1<1) left1=1; end
                   if (right1>n) right1=n; end
                   q=1;A=zeros();%sum=0;
                   for x1=up1:down1
                       for y1=left1:right1
%                            if (x1~=x||y1~=y)
                               A(q)=I1(x1,y1);
%                                sum=sum+(I1(x1,y1)-I1(x,y))^2;
                               q=q+1;
%                            end
                       end
                   end
%                   C_var=sqrt(sum/(q-1));
                  C_var=var(A);
                  C_means=mean(A);
                  if (C_means~=0)
                      C1(k)=C_var/C_means^2;
                  else
                      C1(k)=0;
                  end
                  C_sum=C_sum+C1(k);
                  k=k+1;
                end
            end
        end
%         C_mean=C_sum/k-1;
        C_mean=(C_sum+C_center)/k;
        sum=0;
        for q=1:k-1
                C{i,j}(1,q)=exp((-1)*(C1(q)-C_mean));
                sum=sum+C{i,j}(1,q);
        end
        sum=sum+exp((-1)*(C_center-C_mean));
         for q=1:k-1
             if(C1(q)<C_mean)
                C{i,j}(1,q)=w+C{i,j}(1,q)/sum;
             else
                C{i,j}(1,q)=w-C{i,j}(1,q)/sum;
             end
         end
        
    end
end