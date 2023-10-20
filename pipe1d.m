function [Tout,Pout,Tw_new] = pipe1d(mg,Tin,Pin,Dout,Din,L,Te,Tw,step_L,step_t) 
%pipe1d 用于计算一维圆截面直管道的非稳态流动传热过程
%   输入列表——————————
% 入口质量流量mg，单位kg/s；
% 入口温度Tin，单位K；
% 入口压力Pin，单位Pa；
% 外直径Dout，单位m；
% 内直径Din，单位m；
% 长度L，单位m；
% 环境温度：Te，单位K；认为外部漏热是环境温度和壁面温度的函数，而不是常数；
% 外壁面初始温度：Tw,单值或者一维数组，单位K；
% 空间步长：step_L,默认0.01m；
% 时间步长：step_t,默认1s；
% 
% 输出列表———————————
% 出口温度Tout：一维数组，数组的长度为n。
% 出口压力Pout：一维数组，数组长度为n。单位kPa；
% 一个时间步后的壁面温度Tw_new：一维数组。
%%调用refprop.m等3个文件,steel304.m,G10.m两个物性文件,被main函数调用
if nargin<8
    error("输入变量不能少于8个！")
end

if nargin<9|isempty(step_L)
    step_L=0.01;
end

if nargin<10|isempty(step_t)
    step_t=1; 
end

n=ceil(L/step_L); %总分段数，向上取整。
L_last=L-(n-1)*step_L; %最后一段的长度。防止除不尽的情况。

if length(Tw)==1
    Tw(1:n)=Tw;  %如果只有一个初值，则给出全局外壁温度
end 

q_leak(1:n)=G10(Tw(1:n)/2).*(Te-Tw(1:n)); %改写成G10热导率的函数了。
% q_leak(1:n)=0.5*(Te-Tw(1:n));  %给出漏热热流密度，这里的常数是一个系数，单位和对流换热系数一致（W/m2*K）

Tin(1:n)=Tin;
Pin(1:n)=Pin/1000; %预分配内存，同时初始化,同时把压强单位变成KPa，方便refprop调用
Tout(1:n)=0;Pout(1:n)=0;
Tw_inner(1:n)=Tw(1:n);
QQ(1:n)=0; %每段在一个时间步长内的吸热量

for i=1:n-1
    %计算气体物性，这里得改一下refpropm.m里的路径，把c盘路径挪到本地来。
    %另外为了节省计算时间，不用出入口平均值算物性了，直接用入口来算，结果差不多，省去迭代过程。
    fluid1='air.ppf';
    fluid2='water';
    fluid3='helium';
   
    Cp=refpropm('C','T',Tin(i),'P',Pin(i),fluid3); %Cp [J/(kg K)]
    Den=refpropm('D','T',Tin(i),'P',Pin(i),fluid3); %Density [kg/m^3]
    Cond=refpropm('L','T',Tin(i),'P',Pin(i),fluid3); %Thermal conductivity [W/(m K)]
    Vis=refpropm('V','T',Tin(i),'P',Pin(i),fluid3); %Dynamic viscosity [Pa*s]
    
    Vel=(mg/Den)/(pi*(Din/2)*(Din/2)); % 本段流动速度，m/s
    Re=Vel*Den*Din/Vis; %本段雷诺数
    Pr=Cp*Vis/Cond; %本段普朗特数
    
%%%%%%%%%%%%%%计算压降%%%%%%%%%%%%%%%%%%%%%
    Rug=5e-5; %普通管内粗糙度，抛光管取5e-7.
    XA=(2.457*log(1/((7/Re)^0.9+0.27*Rug/Din)))^16;%系数A
    XB=(37530/Re)^16 ; %系数B
    XF=8*((8/Re)^12+1/(XA+XB)^(3/2))^(1/12); %摩擦系数f
    Pout(i)=(Pin(i)*1000-0.5*XF*step_L/Din*Den*Vel^2)/1000; %出口压力同样化成kpa单位，方便下一轮调用。  
    Pin(i+1)=Pout(i);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if (Re>0)&&(Re<=2300)
        Nu=4;
    elseif Re>4000
        Nu=0.023*Re^0.8*Pr^0.4;
    elseif (Re>2300)&&(Re<=4000)
        Nu=(4^16+(0.023*Re^0.8*Pr^0.4)^16)^(1/16);
    else
        error("雷诺数不能小于0，请检查输入数据！")
    end
    
    h=Nu*Cond/Din; %本段的对流换热系数
    [Cond_S,Cp_S,Den_S]=steel304(Tw(i)); %给出304金属管的热导率，热容,密度。
%     A11=h*Din+2*Cond_S/log(Dout/Din);
%     A12=-h/2*Din;
%     b1=-q_leak(i)*Dout+Tw(i)*2*Cond_S/log(Dout/Din)+h*Tin(i)*Din/2;
%     A21=h*pi*Din*step_L;
%     A22=-h*pi*Din*step_L/2-mg*Cp;
%     b2=(h*pi*Din*step_L/2-mg*Cp)*Tin(i);

    A11=h*Din+2*Cond_S/log(Dout/Din);
    A12=-h/2*Din;
    b1=-q_leak(i)*Dout+Tw(i)*2*Cond_S/log(Dout/Din)+h*Tin(i)*Din/2;
    A21=h*pi*Din*step_L;
    A22=-h*pi*Din*step_L/2-mg*Cp;
    b2=(h*pi*Din*step_L/2-mg*Cp)*Tin(i)-q_leak(i)*Dout*pi*step_L;

    AA=[A11,A12;A21,A22];   %列了两个方程，用矩阵除法求解线性方程组，得到管壁内壁温和本段出口温度。
    bb=[b1;b2];
    YY=AA\bb;
    Tw_inner(i)=YY(1);Tout(i)=YY(2);
    Tin(i+1)=Tout(i);
    QQ(i)=mg*Cp*(Tout(i)-Tin(i))-q_leak(i)*Dout*pi*step_L;
    
    %%%%%%%%%%%%%%计算一个时间步长之后的外壁面温度%%%%%%%%%%%%
    Tw_new(i)=Tw(i)-QQ(i)*step_t/((Dout^2-Din^2)/4*pi*step_L*Cp_S*Den_S);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end


   %算完前面n-1段，算最后一段，最后一段是为了防止总长度不能整除分段长度而预留的一段。     
    fluid1='air.ppf';
    fluid2='water';
    fluid3='helium';
    
    Cp=refpropm('C','T',Tin(i),'P',Pin(n),fluid3);
    Den=refpropm('D','T',Tin(i),'P',Pin(n),fluid3);
    Cond=refpropm('L','T',Tin(i),'P',Pin(n),fluid3);
    Vis=refpropm('V','T',Tin(i),'P',Pin(n),fluid3);
    
    Vel=(mg/Den)/(pi*(Din/2)*(Din/2)); 
    Re=Vel*Den*Din/Vis; 
    Pr=Cp*Vis/Cond; 
%%%%%%%%%%%%%%计算压降%%%%%%%%%%%%%%%%%%%%%
    Rug=5e-5; %普通管内粗糙度，抛光管取5e-7.
    XA=(2.457*log(1/((7/Re)^0.9+0.27*Rug/Din)))^16;
    XB=(37530/Re)^16 ; 
    XF=8*((8/Re)^12+1/(XA+XB)^(3/2))^(1/12); 
    Pout(n)=(Pin(n)*1000-0.5*XF*L_last/Din*Den*Vel^2);   %这里就不用除以1000了，出来的单位直接是Pa。
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    if (Re>0)&&(Re<=2300)
        Nu=4;
    elseif Re>4000
        Nu=0.023*Re^0.8*Pr^0.4;
    elseif (Re>2300)&&(Re<=4000)
        Nu=(4^16+(0.023*Re^0.8*Pr^0.4)^16)^(1/16);
    else
        error("雷诺数不能小于0，请检查输入数据！")
    end
    
    h=Nu*Cond/Din; %本段的对流换热系数
    [Cond_S,Cp_S]=steel304(Tw(n)); %给出304金属管的热导率，热容。
%     A11=h*Din+2*Cond_S/log(Dout/Din);
%     A12=-h/2*Din;
%     b1=-q_leak(n)*Dout+Tw(n)*2*Cond_S/log(Dout/Din)+h*Tin(n)*Din/2;
%     A21=h*pi*Din*L_last;
%     A22=-h*pi*Din*L_last/2-mg*Cp;
%     b2=(h*pi*Din*L_last/2-mg*Cp)*Tin(n);
    
    A11=h*Din+2*Cond_S/log(Dout/Din);
    A12=-h/2*Din;
    b1=-q_leak(n)*Dout+Tw(n)*2*Cond_S/log(Dout/Din)+h*Tin(n)*Din/2;
    A21=h*pi*Din*L_last;
    A22=-h*pi*Din*L_last/2-mg*Cp;
    b2=(h*pi*Din*L_last/2-mg*Cp)*Tin(n)-q_leak(n)*Dout*pi*L_last;
    
    AA=[A11,A12;A21,A22];   %列了两个方程，用矩阵除法求解线性方程组，得到管壁内壁温和本段出口温度。
    bb=[b1;b2];
    YY=AA\bb;
    Tw_inner(n)=YY(1);Tout(n)=YY(2);
    QQ(n)=mg*Cp*(Tout(n)-Tin(n))-q_leak(n)*Dout*pi*L_last;
    Tw_new(n)=Tw(n)-QQ(n)*step_t/((Dout^2-Din^2)/4*pi*L_last*Cp_S*Den_S);
end

