%2012 12 16 by lichao
%2012 12 24 修改，克服了黑色网格化现象。
%传统相机成像
%距透镜d的物体，在透镜后v处成的像
%适用与近轴光学，即物体离透镜的距离远大于焦距
%可分别令v=16,17,0,17.094,15.84

%obj:           场景
%d:             场景距离主透镜的距离
%D:             主透镜直径
%F:             主透镜焦距
%v:             成像位置，即传感器位置
%N_line:        主透镜采样率
%sen_d:         传感器直径
%sen_N:         传感器个数
%nx:            场景x坐标
%ny:            场景y坐标

%%
clc
clear all
close all

addpath(genpath('.\image'));
addpath(genpath('.\sub'));

%%  参数 单位mm
D=4/1.4;                                                                   %直径 1/F = 1.4 = 焦距/直径
F=4;                                                                       %焦距
d=2000;                                                                    %物距
%v=F*d/(d-F);  %像距
v=4.01;                                                                       %在焦距上
N_line=5;                                                                  %对主透镜的离散化，每个点在主透镜有多少光线
sen_N=10000;                                                                 %传感器总个数
sen_d=D/sen_N;                                                             %传感器直径

%选择图像 暂时不需要通过选择的方式实现，而且这个名字也比较奇怪
% obj=select_obj();                                                        %选择物体
% if obj==0
%     clearall
%     disp('请重新选择参数！');
%     break;
% end

%% 读取图像，并将图像转为灰度图
obj_imread=imread('PeppersRGB.tif'); %843*1268
figure
imshow(obj_imread,[]);
obj_r = obj_imread(:,:,1)
obj_g = obj_imread(:,:,2)
obj_b = obj_imread(:,:,3)
figure
imshow(obj_r,[]);
imshow(obj_g,[]);
imshow(obj_b,[]);
for i_rgb =1:1:1

    
    
obj_gray=obj_imread(:,:,i_rgb);
% obj_gray=rgb2gray(obj_imread);
obj_gray_double=double(obj_gray);

obj = obj_gray_double;

%% 模拟成像过程
tic
%[im,error]=LF_sim(obj,d,D,F,v,N_line,sen_d,sen_N);
im=zeros(sen_N);                                                            %清空sen_N 传感器
%%  物体信息
Kx=size(obj,1);                                                             %返回行数
Ky=size(obj,2);                                                             %返回列数
d_D=d/D;                                                                    %场景是主透镜的多少倍
d_m=1;                                                                      %每个传感器上对应场景中的点的个数
%[nx,ny,error]=select_nx_ny(Kx,Ky,F,d,sen_N,D,d_D,d_m);                      %d_D,d_m为可选参数，不输入时为d_D=8,d_m=2
dd=(D/sen_N)*(d-F)/F; 
% dd=10; 
%表示每个场景间的最大间隔
nx=linspace(-dd*Kx/(2*d_m),dd*Kx/(2*d_m),Kx);
ny=linspace(-dd*Ky/(2*d_m),dd*Ky/(2*d_m),Ky);
%     disp('不呈现黑色网格现象：');
    disp(['每个传感器对应的场景像素个数为：d_m=',num2str(d_m)]);
% if error==1
%     return;
% else
    fprintf('\n开始相机传输过程：\n');
% end
%% 相机成像过程 薄透镜
M_len=[1 0;-1/F 1];                                                         %主透镜折射
%这个是什么问题？ 怎么传输，是主透镜到sensor之间
T=[1 v;0 1];                                                                %主透镜与微透镜之间传播 

sen_N_1=1/sen_d;                                                            %单位距离传感器个数

para_x=zeros(N_line,1);
para_rx=zeros(N_line,1);                                                    %光线在传感器上坐标
para_th=zeros(N_line,1);                                                    %光线与微透镜的夹角

para_y=zeros(N_line,1);
para_ry=zeros(N_line,1);
para_be=zeros(N_line,1);

n=linspace(-fix(D/2),fix(D/2),N_line);%透镜离散化

k_num=0.1;
for ix=1:Kx
    for jy=1:Ky
        if obj(ix,jy)~=0
            dx=nx(ix);                                                      %选择入射光在物体处的x坐标
            dy=ny(jy);                                                      %选择入射光在物体处的y坐标
           %% 每条光线的光强  此处有近似
            In=obj(ix,jy)/(N_line ^2);  %%不清楚光强和物体有什么关系？可以理解为灰度值，就是光强，每一个点上的光束都均匀的射到镜面上，
                                        %%后续就按照镜面为入射面进行计算就可以额
                                        %%镜面上物[1,0;0,1]可以不乘，就得到T*M*[R,theta]；
                                        %%这样基本上能说的通了
          
            for i=1:N_line
                x=n(i);                                                     %透镜上离散坐标
                theta=atan((x-dx)/d);                                       %入射光与透镜的角度
              %  pq1_x=T*M_len*[x;theta];      %[x;theta]表示入射光线向量，和z高度，以及和z水平夹角             
              %透镜的转换及移动距离v
                pq1_x=T*M_len*[x;theta];
                x1=pq1_x(1);
                para_th(i)=pq1_x(2);
                para_x(i)=sen_N_1*x1+fix(sen_N/2)+1 ;                       %调整后的x2  1--4001
                para_rx(i)=round(para_x(i));                                %传感器坐标为整数
               
            end
            for j=1:N_line
                y=n(j);
                beta=atan((y-dy)/d);
                pq1_y=T*M_len*[y;beta];                                     %透镜的转换及移动距离v
                y1=pq1_y(1);
                para_be(j)=pq1_y(2);
                para_y(j)=sen_N_1*y1+fix(sen_N/2)+1 ;                       %调整后的y2
                para_ry(j)=round(para_y(j));                                %传感器坐标为整数
            end
            
            %% 计算光强
            for i=1:N_line
                for j=1:N_line
                    if ((n(i)^2+n(j)^2)<=(D/2)^2)&&(para_rx(i)<=sen_N)&&(para_ry(j)<=sen_N)&&(para_rx(i)>0)&&(para_ry(j)>0)
                        %((n(i)^2+n(j)^2)<=(D/2)^2)                                     主透镜圆形
                        %(para_rx(i)<=sen_N_total)&&(para_ry(j)<=sen_N_total)&&(para_rx(i)>0)&&(para_ry(j)>0) 在传感器范围内
                        rx=para_rx(i);
                        ry=para_ry(j);
                        midpara=im(rx,ry);%% 变换光强？
                        In1=midpara+In*sqrt(cos(para_th(i))^2+cos(para_be(j))^2);  %%理解为x方面和y方向的光强 ，最终合并为z方向的光强
                        im(rx,ry)=In1;
                    end
                end
            end
        end
    end
    if ix==fix(k_num*Kx)
        disp(['已完成',num2str(k_num*100),'%']);
        k_num=k_num+0.1;
    end
end

t=toc;
disp(['模拟成像过程花费时间为t=',num2str(t),'s']);
% 
% if error==1
%     clearall
%     disp('请重新选择参数！');
%     break;
% end
%%  画出图像
figure
subplot(1,2,1),imshow(obj,[]);title('原始图像');
subplot(1,2,2),imshow(im,[]);title('传感器图像im');
disp('已画出传感器图像！');

im_revi=sub_revise_im(im);
%figure
%imshow(im_revi,[]);title('消除0行0列的图像')

im_reve=sub_reversal_im(im_revi);
figure
imshow(im_reve,[]);title('反转后的图像')
if i_rgb ==1
    result1 = im_reve;
end 
if i_rgb ==2
    result2 = im_reve;
end
if i_rgb ==3
    result3 = im_reve;
end 
end
%%以上代码主要是显示灰度图，需要按照各个光线的灰度再处理一次。



%%完成灰度图处理之后，可以可以实现bayer数据以及YUV数据的转换

%%YUV 数据完成后可以开始编码 Mjpeg，尝试将多个图像合并为一个视屏来编码H.264

%%编码完成后，可以开始打包，并发送码流