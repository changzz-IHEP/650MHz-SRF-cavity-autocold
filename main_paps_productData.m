clc;clear;
format long

%% ����������
%   �����б�������������������
% �����������mg����λkg/s��
% ����¶�Tin����λK��
% ���ѹ��Pin����λPa��
% ��ֱ��Dout����λm��
% ��ֱ��Din����λm��
% ����L����λm��
% �����¶ȣ�Te����λK����Ϊ�ⲿ©���ǻ����¶Ⱥͱ����¶ȵĺ����������ǳ�����
% ������ʼ�¶ȣ�Tw,��ֵ����һά���飬��λK��
% �ռ䲽����step_L,���Բ��䣬Ĭ��0.01m��
% ʱ�䲽����step_t,���Բ��䣬Ĭ��1s��
% 
% ����б���������������������
% �����¶�Tout��һά���飬����ĳ���Ϊn��
% ����ѹ��Pout��һά���飬���鳤��Ϊn��
% һ��ʱ�䲽��ı����¶�Tw_new��һά���顣

Pin=166325;%pa
Dout=0.0269;%m
Din=0.02268;%m
Din2=0.0236;
L1=4.7; %m
L_b=33;

step_L=1;%m
step_t=20; %s

%%
%����¶ȣ���������������¶�,�����¶�4��������4��ѭ�����Ӷ���������ʼѭ��
file1 = fopen('CalData_BTCM.txt','a');
a=20;b=7;c=16;d=10; %���ɱ߽�����
% a=2;b=2;c=2;d=2; %���Դ��룬16��߽�
for i=1:a
    mg=0.0003*i;
    for j=1:b
        Te=268+5*j;
        for k=1:c
            TwS=Te-(Te-100)/15*(k-1);
            for m=1:d
                Tin=TwS-6*m; 
       
%% 
n=360; %��ʱ��Ϊn*step_t
To1=zeros(1,n);
x=1:n;T1d=zeros(1,n);
Tw=TwS;Tw_b=TwS; %����ֵ���ñ����¶Ȼָ�
for t=1:n
[Tout,Pout,Tw_new]=pipe1d(mg,Tin,Pin,Dout,Din,L1,Te,Tw,step_L,step_t);
Tw=Tw_new;

Tin_b=Tout(end);Pin_b=Pout(end);
[Tout_b,Pout_b,Tw_new_b]=pipe1d(mg,Tin_b,Pin_b,Dout,Din2,L_b,Te,Tw_b,step_L,step_t);
Tw_b=Tw_new_b;

To1(t)=Tout_b(end);
i,j,k,m
t
To1(t)

    if (t>1)
    T1d(t)=(To1(t)-To1(t-1))/step_t;
    end
end
%%
%����ı�
mg1(1:n)=mg;Tin1(1:n)=Tin;Te1(1:n)=Te;Tw1(1:n)=TwS; %�������������¶ȣ������¶�,�����¶ȡ�
S=[step_t*x;mg1;Te1;Tw1;Tin1;To1;T1d]; %�������Ϣ����һ�������ʱ�䣬�¶ȵȲ���

fprintf(file1,'%8.1f %9.6f %9.6f %9.6f %9.6f %9.6f %12.8f\n',S);
% fclose(file1);

            end
        end
    end 
end
fclose(file1);