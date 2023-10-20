clc;clear;close all;
format long
%%
Pin=156325;%pa
Dout=0.0269;%m
Din=0.02268;%m
Din2=0.0236;
L1=4.7; %m
L_b=40; %到相分离器入�?33，到超导腔按40�?

% L_tank_max=0.4515;
% L_tank_min=0.4399;
% Do_tank_max=0.1475;
% Do_tank_min=0.101;
% e_tank_max=0.04532;
% e_tank_min=0.01684;
%相分离器尺寸
% L_tank_separator=0.76;
% Do_tank_separator=0.41;
% e_tank_separator=0.01; %0.005

L_tank_max=0.4515;
L_tank_min=0.4399;
Do_tank_max=0.3;
Do_tank_min=0.15;
e_tank_max=0.09;
e_tank_min=0.04;
L_tank_separator=0.76;
Do_tank_separator=0.41;
e_tank_separator=0.05; %0.005

step_L=1;%m
step_t=10; %s
mg=0.004;
Te=300;
%入口温度，壁面温�?2个变量跑2重循环，流量固定�?4g
file1 = fopen('CalData_BTCM_deltaT.txt','a');
a=18;b=12; %生成边界条件
% a=2;b=2; %测试代码�?4组边�?

 for i=1:a
     TwS=300-10*i;
     for j=1:b
        Tin=TwS-5*j; 
              
        Tw=TwS;Tw_b=TwS; TwS_max=TwS;TwS_min=TwS;TwS_separator=TwS;%给初值，让壁面温度恢�?
        DdeltaT=1;deltaTmax=0;
        t=1; %时间�?
        while DdeltaT>=0 || deltaTmax<0.001 
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
            deltaT(t)=Tw_max(t)-Tw_min(t);
            if t>1
            DdeltaT=deltaT(t)-deltaT(t-1);
            end

            deltaTmax=deltaT(t);
            t=t+1
            i,j
            DdeltaT,deltaTmax
            
            if t==14400
                break
            end
        end
        deltaTmax2=deltaT(t-1); %入口温度，壁面温度�??
        S=[TwS,Tin,deltaTmax2,step_t*t]; %把相关信息放入一个向量里
        fprintf(file1,'%4.1f %4.1f %9.6f %9.1f\n',S);            
     end
 end
fclose(file1);