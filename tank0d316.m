function [Tout,Tw_new] = tank0d316(mg,Tin,Pin,do,e,L,Tw,step_t) 
%pipe1d 用于计算0维圆截面tank罐的非稳态流动传热过程
%   输入列表——————————
% 入口质量流量mg，单位kg/s；
% 入口温度Tin，单位K；
% 入口压力Pin，单位Pa；
% 外直径Dout，单位m；
% 厚度，单位m；
% 长度L，单位m；
% 外壁面初始温度：Tw,单位K；
% 时间步长：step_t,默认1s；
% 
% 输出列表———————————
% 一个时间步后的壁面温度Tw_new，K。
%%调用refprop.m等3个文件,被main函数调用
if nargin<7
    error("输入变量不能少于7个！")
end

if nargin<8|isempty(step_t)
    step_t=1; 
end

di=do-2*e;
Qin=mg;
Pin=Pin/1000; %Pa化kpa，方便refprop调用

TT=log10(Tw);
aa=-1879.464;
bb=3643.198;
cc=76.70125;
dd=-6176.028;
ee=7437.6247;
ff=-4305.7217;
gg=1382.4627;
hh=-237.22704;
ii=17.05262;
Cp_316=10^(aa+bb*TT+cc*TT^2+dd*TT^3+ee*TT^4+ff*TT^5+gg*TT^6+hh*TT^7+ii*TT^8);     %Cp/(J/kg/K) 316L制造的tank罐 300K-50K

rou_316=7980;%rou(kg/m^3)
V=pi/4*(do^2-di^2)*L;
M=rou_316*V;

char refpropm;
fluid='helium';
lamda=refpropm('L','T',Tin,'P',Pin,fluid); %%L代表Thermal conductivity [W/(m K)]
rou=refpropm('D','T',Tin,'P',Pin,fluid);   %D代表density [kg/m3]
miu=refpropm('V','T',Tin,'P',Pin,fluid);   %V代表viscosity [Pa*s]
Cp=refpropm('C','T',Tin,'P',Pin,fluid);   %C代表Cp [J/kg*K]
Pr=refpropm('^','T',Tin,'P',Pin,fluid);   %^   Prandtl number [-]
u=Qin/(rou*(pi*di^2/4));
Re=u*di*rou/miu;
Nu_tur=0.023*Re^0.8*Pr^0.4;
Nu_lam=4;
Nu_jet=2.2*(Qin/(pi*di*miu)*Pr)^0.632;
Nu=(Nu_lam^16+Nu_tur^16+Nu_jet^16)^(1/16);

h=lamda*Nu/di;
A=pi*di*L;
Tout=(2*Cp*Qin*Tin+2*h*A*Tw-h*A*Tin)/(2*Cp*Qin+h*A);
Tw_new=Tw-((Cp*Qin*(Tout-Tin)*step_t)/(Cp_316*M));
end