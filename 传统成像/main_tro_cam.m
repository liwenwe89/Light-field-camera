%2012 12 16 by lichao
%2012 12 24 �޸ģ��˷��˺�ɫ��������
%��ͳ�������
%��͸��d�����壬��͸����v���ɵ���
%����������ѧ����������͸���ľ���Զ���ڽ���
%�ɷֱ���v=16,17,0,17.094,15.84


%%
clc
clear all
close all

addpath(genpath('.\image'));
addpath(genpath('.\sub'));

%%  ���� ��λmm
D=4/1.4;                                                                   %ֱ�� 1/F = 1.4 = ����/ֱ��
F=4;                                                                       %����
d=2000;                                                                    %���
%v=F*d/(d-F);  %���
v=4;                                                                       %�ڽ�����
N_line=5;                                                                  %����͸������ɢ����ÿ��������͸���ж��ٹ���
sen_N=500;                                                                 %�������ܸ���
sen_d=D/sen_N;                                                             %������ֱ��

%ѡ��ͼ�� ��ʱ����Ҫͨ��ѡ��ķ�ʽʵ�֣������������Ҳ�Ƚ����
% obj=select_obj();                                                        %ѡ������
% if obj==0
%     clearall
%     disp('������ѡ�������');
%     break;
% end

%% ��ȡͼ�񣬲���ͼ��תΪ�Ҷ�ͼ
obj_imread=imread('17.jpg'); %843*1268
obj_gray=rgb2gray(obj_imread);
obj_gray_double=double(obj_gray);

obj = obj_gray_double;

%% ģ��������
tic
[im,error]=LF_sim(obj,d,D,F,v,N_line,sen_d,sen_N);
t=toc;
disp(['ģ�������̻���ʱ��Ϊt=',num2str(t),'s']);

if error==1
    clearall
    disp('������ѡ�������');
    break;
end
%%  ����ͼ��
figure
subplot(1,2,1),imshow(obj,[]);title('ԭʼͼ��');
subplot(1,2,2),imshow(im,[]);title('������ͼ��im');
disp('�ѻ���������ͼ��');

im_revi=sub_revise_im(im);
figure
imshow(im_revi,[]);title('����0��0�е�ͼ��')

im_reve=sub_reversal_im(im_revi);
figure
imshow(im_reve,[]);title('��ת���ͼ��')

