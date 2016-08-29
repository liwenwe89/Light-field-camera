%2012 12 16 by lichao
%2012 12 24 修改，克服了黑色网格化现象。
%传统相机成像
%距透镜d的物体，在透镜后v处成的像
%适用与近轴光学，即物体离透镜的距离远大于焦距
%可分别令v=16,17,0,17.094,15.84


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
v=4;                                                                       %在焦距上
N_line=5;                                                                  %对主透镜的离散化，每个点在主透镜有多少光线
sen_N=500;                                                                 %传感器总个数
sen_d=D/sen_N;                                                             %传感器直径

%选择图像 暂时不需要通过选择的方式实现，而且这个名字也比较奇怪
% obj=select_obj();                                                        %选择物体
% if obj==0
%     clearall
%     disp('请重新选择参数！');
%     break;
% end

%% 读取图像，并将图像转为灰度图
obj_imread=imread('17.jpg'); %843*1268
obj_gray=rgb2gray(obj_imread);
obj_gray_double=double(obj_gray);

obj = obj_gray_double;

%% 模拟成像过程
tic
[im,error]=LF_sim(obj,d,D,F,v,N_line,sen_d,sen_N);
t=toc;
disp(['模拟成像过程花费时间为t=',num2str(t),'s']);

if error==1
    clearall
    disp('请重新选择参数！');
    break;
end
%%  画出图像
figure
subplot(1,2,1),imshow(obj,[]);title('原始图像');
subplot(1,2,2),imshow(im,[]);title('传感器图像im');
disp('已画出传感器图像！');

im_revi=sub_revise_im(im);
figure
imshow(im_revi,[]);title('消除0行0列的图像')

im_reve=sub_reversal_im(im_revi);
figure
imshow(im_reve,[]);title('反转后的图像')

