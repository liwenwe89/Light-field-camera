%2012 12 16 by lichao
%2012 12 24 �޸ģ��˷��˺�ɫ��������
%��ͳ�������
%��͸��d�����壬��͸����v���ɵ���
%����������ѧ����������͸���ľ���Զ���ڽ���
%�ɷֱ���v=16,17,0,17.094,15.84

%obj:           ����
%d:             ����������͸���ľ���
%D:             ��͸��ֱ��
%F:             ��͸������
%v:             ����λ�ã���������λ��
%N_line:        ��͸��������
%sen_d:         ������ֱ��
%sen_N:         ����������
%nx:            ����x����
%ny:            ����y����

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
v=4.01;                                                                       %�ڽ�����
N_line=5;                                                                  %����͸������ɢ����ÿ��������͸���ж��ٹ���
sen_N=10000;                                                                 %�������ܸ���
sen_d=D/sen_N;                                                             %������ֱ��

%ѡ��ͼ�� ��ʱ����Ҫͨ��ѡ��ķ�ʽʵ�֣������������Ҳ�Ƚ����
% obj=select_obj();                                                        %ѡ������
% if obj==0
%     clearall
%     disp('������ѡ�������');
%     break;
% end

%% ��ȡͼ�񣬲���ͼ��תΪ�Ҷ�ͼ
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

%% ģ��������
tic
%[im,error]=LF_sim(obj,d,D,F,v,N_line,sen_d,sen_N);
im=zeros(sen_N);                                                            %���sen_N ������
%%  ������Ϣ
Kx=size(obj,1);                                                             %��������
Ky=size(obj,2);                                                             %��������
d_D=d/D;                                                                    %��������͸���Ķ��ٱ�
d_m=1;                                                                      %ÿ���������϶�Ӧ�����еĵ�ĸ���
%[nx,ny,error]=select_nx_ny(Kx,Ky,F,d,sen_N,D,d_D,d_m);                      %d_D,d_mΪ��ѡ������������ʱΪd_D=8,d_m=2
dd=(D/sen_N)*(d-F)/F; 
% dd=10; 
%��ʾÿ��������������
nx=linspace(-dd*Kx/(2*d_m),dd*Kx/(2*d_m),Kx);
ny=linspace(-dd*Ky/(2*d_m),dd*Ky/(2*d_m),Ky);
%     disp('�����ֺ�ɫ��������');
    disp(['ÿ����������Ӧ�ĳ������ظ���Ϊ��d_m=',num2str(d_m)]);
% if error==1
%     return;
% else
    fprintf('\n��ʼ���������̣�\n');
% end
%% ���������� ��͸��
M_len=[1 0;-1/F 1];                                                         %��͸������
%�����ʲô���⣿ ��ô���䣬����͸����sensor֮��
T=[1 v;0 1];                                                                %��͸����΢͸��֮�䴫�� 

sen_N_1=1/sen_d;                                                            %��λ���봫��������

para_x=zeros(N_line,1);
para_rx=zeros(N_line,1);                                                    %�����ڴ�����������
para_th=zeros(N_line,1);                                                    %������΢͸���ļн�

para_y=zeros(N_line,1);
para_ry=zeros(N_line,1);
para_be=zeros(N_line,1);

n=linspace(-fix(D/2),fix(D/2),N_line);%͸����ɢ��

k_num=0.1;
for ix=1:Kx
    for jy=1:Ky
        if obj(ix,jy)~=0
            dx=nx(ix);                                                      %ѡ������������崦��x����
            dy=ny(jy);                                                      %ѡ������������崦��y����
           %% ÿ�����ߵĹ�ǿ  �˴��н���
            In=obj(ix,jy)/(N_line ^2);  %%�������ǿ��������ʲô��ϵ���������Ϊ�Ҷ�ֵ�����ǹ�ǿ��ÿһ�����ϵĹ��������ȵ��䵽�����ϣ�
                                        %%�����Ͱ��վ���Ϊ��������м���Ϳ��Զ�
                                        %%��������[1,0;0,1]���Բ��ˣ��͵õ�T*M*[R,theta]��
                                        %%������������˵��ͨ��
          
            for i=1:N_line
                x=n(i);                                                     %͸������ɢ����
                theta=atan((x-dx)/d);                                       %�������͸���ĽǶ�
              %  pq1_x=T*M_len*[x;theta];      %[x;theta]��ʾ���������������z�߶ȣ��Լ���zˮƽ�н�             
              %͸����ת�����ƶ�����v
                pq1_x=T*M_len*[x;theta];
                x1=pq1_x(1);
                para_th(i)=pq1_x(2);
                para_x(i)=sen_N_1*x1+fix(sen_N/2)+1 ;                       %�������x2  1--4001
                para_rx(i)=round(para_x(i));                                %����������Ϊ����
               
            end
            for j=1:N_line
                y=n(j);
                beta=atan((y-dy)/d);
                pq1_y=T*M_len*[y;beta];                                     %͸����ת�����ƶ�����v
                y1=pq1_y(1);
                para_be(j)=pq1_y(2);
                para_y(j)=sen_N_1*y1+fix(sen_N/2)+1 ;                       %�������y2
                para_ry(j)=round(para_y(j));                                %����������Ϊ����
            end
            
            %% �����ǿ
            for i=1:N_line
                for j=1:N_line
                    if ((n(i)^2+n(j)^2)<=(D/2)^2)&&(para_rx(i)<=sen_N)&&(para_ry(j)<=sen_N)&&(para_rx(i)>0)&&(para_ry(j)>0)
                        %((n(i)^2+n(j)^2)<=(D/2)^2)                                     ��͸��Բ��
                        %(para_rx(i)<=sen_N_total)&&(para_ry(j)<=sen_N_total)&&(para_rx(i)>0)&&(para_ry(j)>0) �ڴ�������Χ��
                        rx=para_rx(i);
                        ry=para_ry(j);
                        midpara=im(rx,ry);%% �任��ǿ��
                        In1=midpara+In*sqrt(cos(para_th(i))^2+cos(para_be(j))^2);  %%���Ϊx�����y����Ĺ�ǿ �����պϲ�Ϊz����Ĺ�ǿ
                        im(rx,ry)=In1;
                    end
                end
            end
        end
    end
    if ix==fix(k_num*Kx)
        disp(['�����',num2str(k_num*100),'%']);
        k_num=k_num+0.1;
    end
end

t=toc;
disp(['ģ�������̻���ʱ��Ϊt=',num2str(t),'s']);
% 
% if error==1
%     clearall
%     disp('������ѡ�������');
%     break;
% end
%%  ����ͼ��
figure
subplot(1,2,1),imshow(obj,[]);title('ԭʼͼ��');
subplot(1,2,2),imshow(im,[]);title('������ͼ��im');
disp('�ѻ���������ͼ��');

im_revi=sub_revise_im(im);
%figure
%imshow(im_revi,[]);title('����0��0�е�ͼ��')

im_reve=sub_reversal_im(im_revi);
figure
imshow(im_reve,[]);title('��ת���ͼ��')
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
%%���ϴ�����Ҫ����ʾ�Ҷ�ͼ����Ҫ���ո������ߵĻҶ��ٴ���һ�Ρ�



%%��ɻҶ�ͼ����֮�󣬿��Կ���ʵ��bayer�����Լ�YUV���ݵ�ת��

%%YUV ������ɺ���Կ�ʼ���� Mjpeg�����Խ����ͼ��ϲ�Ϊһ������������H.264

%%������ɺ󣬿��Կ�ʼ���������������