function [nx,ny,error]=select_nx_ny(Kx,Ky,F,d,sen_N,D,d_D,d_m)
%% �Ƿ���ֺ�ɫ��������
error=0;
nx=zeros(1,Kx);
ny=zeros(1,Ky);
grad_num=menu('�Ƿ���ֺ�ɫ��������','��','��','����');

if grad_num==1
    %% �к�ɫ��������
    % ��Ϊÿ���������϶�Ӧ�����ֻ������ĵ�ĸ�����Ϊ����
    if (exist('d_D','var')==0)||(d_D==0)
        d_D=200;                                                             %��������͸���Ķ��ٱ�
    end
    nx=linspace(-d_D*D/2,d_D*D/2,Kx);
    ny=linspace(-d_D*D/2,d_D*D/2,Ky);
    disp('���ֺ�ɫ��������');
    disp(['��������͸��D��d_D����d_D=',num2str(d_D)]);
elseif grad_num==2
    %% �˷��˺�ɫ��������
    % ʹÿ���������϶�Ӧ�ĳ����е����ֻ��ĵ�ĸ���Ϊ����
    %v_re=F*d/(d-F)
    %sen_d/v_re=dd/d
    %dd=d*sen_d/v_re=d*(D/sen_N)*((d-F)/(F*d))=(D/sen_N)*(F-d)/F=(4/500)*(250-16)/16=0.117
    if (exist('d_m','var')==0)||(d_m==0)
        d_m=200;                                                           %ÿ���������϶�Ӧ�����еĵ�ĸ���
    end
    dd=(D/sen_N)*(d-F)/F;                                                  %��ʾÿ��������������
    nx=linspace(-dd*Kx/(2*d_m),dd*Kx/(2*d_m),Kx);
    ny=linspace(-dd*Ky/(2*d_m),dd*Ky/(2*d_m),Ky);
    disp('�����ֺ�ɫ��������');
    disp(['ÿ����������Ӧ�ĳ������ظ���Ϊ��d_m=',num2str(d_m)]);
else
    error=1;
    return;
end