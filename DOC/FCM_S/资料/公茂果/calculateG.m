function [G1,U1]=calculateG(U1,D1,I1,m1,c,C)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%���룺��������V1�������Ⱦ���U1����ԭͼ�����غ��I1
%�����ÿ�����Ӧ��G
%���ܣ�����ÿ�����Ӧÿ���������ĵ�G
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alpha=0.5;
[m,n]=size(I1);
G1=cell(m,n);
U9=cell(m,n);

for i=1:m
    for j=1:n
        G1{i,j}=zeros(c,1);
    end
end

for i=1:m
    for j=1:n
        for k=1:c
           temp=0;
           for x=i-1:i+1
                for y=j-1:j+1 
                    q=1;
                    if x>0&&x<m&&y>0&&y<n&&(x~=i||y~=j)                         
                        d=sqrt((i-x)^2+(j-y)^2); 
%                         temp=temp+1/(1+d)*(1-U1{x,y}(k,1))^m1*(I1(x,y)-V2(k))^2;
%                         temp=temp+C{i,j}(1,q)*(1-U1{x,y}(k,1))^m1*(I1(x,y)-V2(k))^2;
                        temp=temp+1/(1+d)*C{i,j}(1,q)*(1-U1{x,y}(k,1))^m1*(1-D1{x,y}(k,1));%���ٳͷ���
                        q=q+1;
                    end
                end
            end
            G1{i,j}(k,1)=temp;
            G(k)=temp;
        end
        for k=1:c
            temp=0;
            for z=1:c
%                 temp=temp+(G(k)/G(z))^(1/(m1-1));
                 a=1-D1{i,j}(z,1)+G(z);
                 if a~=0     
                    temp=temp+((1-D1{i,j}(k,1)+G(k))/a)^(1/(m1-1));
                 else
                    temp=temp+((1-D1{i,j}(k,1)+G(k))/0.01)^(1/(m1-1)); 
                 end
            end
            if temp~=0
               U1{i,j}(k,1)=1/temp;
            else
               U1{i,j}(k,1)=1;
            end
      
        end
    end
end
%%%%%%%%%%%%%����ʽ�޸�
%        for i=1:m
%          for j=1:n
%              for k=1:c    
%                 U9{i,j}(:,1)=U1{i,j}(:,1).*alpha;
%          
%                 [~,max_loca]=max(U1{i,j}(:,1));
%                 tempU=alpha.*U1{i,j}(max_loca,1)+(1-alpha);
%                 U9{i,j}(max_loca,1)=tempU;
%           
%                 
%              end
%          end
%        end
%        U1=U9;
%%%%%%%%%%%%%%%%%%%

   