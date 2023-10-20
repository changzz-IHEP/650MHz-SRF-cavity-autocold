clc;clear;
format long

%% 参数输入区
%   输入列表——————————
% 入口质量流量mg，单位kg/s；
% 入口温度Tin，单位K；
% 入口压力Pin，单位Pa；
% 外直径Dout，单位m；
% 内直径Din，单位m；
% 长度L，单位m；
% 环境温度：Te，单位K；认为外部漏热是环境温度和壁面温度的函数，而不是常数；
% 外壁面初始温度：Tw,单值或者一维数组，单位K；
% 空间步长：step_L,可以不输，默认0.01m；
% 时间步长：step_t,可以不输，默认1s；
% 
% 输出列表———————————
% 出口温度Tout：一维数组，数组的长度为n。
% 出口压力Pout：一维数组，数组长度为n。
% 一个时间步后的壁面温度Tw_new：一维数组。

Pin=166325;%pa
Dout=0.0269;%m
Din=0.02268;%m
Din2=0.0236;
L1=4.7; %m
L_b=33;

step_L=1;%m
step_t=20; %s

%%
%入口温度，入口流量，环境温度,壁面温度4个变量跑4重循环，从独立变量开始循环
file1 = fopen('CalData_BTCM.txt','a');
a=20;b=7;c=16;d=10; %生成边界条件
% a=2;b=2;c=2;d=2; %测试代码，16组边界
for i=1:a
    mg=0.0003*i;
    for j=1:b
        Te=268+5*j;
        for k=1:c
            TwS=Te-(Te-100)/15*(k-1);
            for m=1:d
                Tin=TwS-6*m; 
       
%% 
n=360; %总时长为n*step_t
To1=zeros(1,n);
x=1:n;T1d=zeros(1,n);
Tw=TwS;Tw_b=TwS; %给初值，让壁面温度恢复
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
%输出文本
mg1(1:n)=mg;Tin1(1:n)=Tin;Te1(1:n)=Te;Tw1(1:n)=TwS; %入口流量，入口温度，环境温度,壁面温度。
S=[step_t*x;mg1;Te1;Tw1;Tin1;To1;T1d]; %把相关信息放入一个矩阵里，时间，温度等参数

fprintf(file1,'%8.1f %9.6f %9.6f %9.6f %9.6f %9.6f %12.8f\n',S);
% fclose(file1);

            end
        end
    end 
end
fclose(file1);