clc;clear;close all;
format long

%% 参数输入�?
%   输入列表—�?��?��?��?��?��?��?��?��??
% 入口质量流量mg，单位kg/s�?
% 入口温度Tin，单位K�?
% 入口压力Pin，单位Pa�?
% 外直径Dout，单位m�?
% 内直径Din，单位m�?
% 长度L，单位m�?
% L_tank_max,代表高温tank罐长度；
% L_tank_min,代表低温tank罐长度；
% Do_tank_max,高温tank罐外直径�?
% Do_tank_min,低温tank罐外直径�?
% e_tank_max,高温tank罐厚度；
% e_tank_min,低温tank罐厚度；
% 环境温度：Te，单位K；认为外部漏热是环境温度和壁面温度的函数，�?�不是常数；
% 外壁面初始温度：Tw,单�?�或者一维数组，单位K�?
% 空间步长：step_L,可以不输，默�?0.01m�?
% 时间步长：step_t,可以不输，默�?1s�?
% 
% 输出列表—�?��?��?��?��?��?��?��?��?��??
% 出口温度Tout：一维数组，数组的长度为n�?
% 出口压力Pout：一维数组，数组长度为n�?
% �?个时间步后的壁面温度Tw_new：一维数组�??

Pin=156325;%pa
Dout=0.0269;%m
Din=0.02268;%m
Din2=0.0236;
L1=4.7; %m
L_b=40; %到相分离器入�?33，到超导腔按40�?

L_tank_max=0.4515;
L_tank_min=0.4399;

Do_tank_max=0.3;
Do_tank_min=0.15;
e_tank_max=0.09;
e_tank_min=0.04;
% Do_tank_max=0.1475;
% Do_tank_min=0.101;
% e_tank_max=0.04532;
% e_tank_min=0.01684;

%相分离器尺寸
L_tank_separator=0.76;
Do_tank_separator=0.41;
e_tank_separator=0.05; %0.005

step_L=1;%m
step_t=10; %s

%%
mg=0.004;
TinS=[260];
% TinS=[271.9
% 268.9
% 254.8
% 244.2
% 239.5
% 232.8
% 227.7
% 225.7
% 222.4
% 216.2
% 211.1
% 206.6
% 201.5
% 196.9
% 191.3
% 185.6
% 179.1
% 174.9
% 169.4
% 161.7
% 158
% 155
% 150.1
% 144.7
% ];
TwS=290;
Te=300;
k=length(TinS);
Tw=TwS;Tw_b=TwS; TwS_max=TwS;TwS_min=TwS;TwS_separator=TwS;%给初值，让壁面温度恢�?
for i=1:k
    Tin=TinS(i);
    
%% 
n=1440; %总时长为n*step_t
To1=zeros(1,n);
Tw_max=zeros(1,n);Tw_min=zeros(1,n);Tw_separator=zeros(1,n);
% x=1:n;

for t=1:n
[Tout,Pout,Tw_new]=pipe1d(mg,Tin,Pin,Dout,Din,L1,Te,Tw,step_L,step_t);
Tw=Tw_new;

Tin_b=Tout(end);Pin_b=Pout(end);
[Tout_b,Pout_b,Tw_new_b]=pipe1d(mg,Tin_b,Pin_b,Dout,Din2,L_b,Te,Tw_b,step_L,step_t);
Tw_b=Tw_new_b;

To1(t)=Tout_b(end);
Tin_separator=Tout_b(end);Pin_separator=Pout_b(end);
[Tout_separator,Tw_separator(t)]=tank0d316(mg,Tin_separator,Pin_separator,Do_tank_separator,e_tank_separator,L_tank_separator,TwS_separator,step_t);
TwS_separator=Tw_separator(t);

Tin_max=Tout_separator;Tin_min=Tout_separator;
Pin_max=Pout_b(end);Pin_min=Pout_b(end);

[Tout_max,Tw_max(t)]=tank0dNb(mg/2,Tin_max,Pin_max,Do_tank_max,e_tank_max,L_tank_max,TwS_max,step_t);
TwS_max=Tw_max(t);
[Tout_min,Tw_min(t)]=tank0dNb(mg/2,Tin_min,Pin_min,Do_tank_min,e_tank_min,L_tank_min,TwS_min,step_t);
TwS_min=Tw_min(t);

i,t
Tw_max(t)
Tw_min(t)
end
To_2(i,:)=To1;
Tw_separator_2(i,:)=Tw_separator;
Tw_max_2(i,:)=Tw_max;
Tw_min_2(i,:)=Tw_min;
end
%%
To_3=reshape(To_2',k*n,1);
Tw_max_3=reshape(Tw_max_2',k*n,1);
Tw_min_3=reshape(Tw_min_2',k*n,1);
Tw_separator_3=reshape(Tw_separator_2',k*n,1);
x=1:n*k;
figure(1)
plot(x*step_t,Tw_max_3,'b-',x*step_t,Tw_min_3,'r-',x*step_t,Tw_separator_3,'g-',x*step_t,To_3,'--')
delta_T=Tw_max_3-Tw_min_3;
figure(2)
plot(x*step_t,delta_T)