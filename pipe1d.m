function [Tout,Pout,Tw_new] = pipe1d(mg,Tin,Pin,Dout,Din,L,Te,Tw,step_L,step_t) 
%pipe1d ���ڼ���һάԲ����ֱ�ܵ��ķ���̬�������ȹ���
%   �����б�������������������
% �����������mg����λkg/s��
% ����¶�Tin����λK��
% ���ѹ��Pin����λPa��
% ��ֱ��Dout����λm��
% ��ֱ��Din����λm��
% ����L����λm��
% �����¶ȣ�Te����λK����Ϊ�ⲿ©���ǻ����¶Ⱥͱ����¶ȵĺ����������ǳ�����
% ������ʼ�¶ȣ�Tw,��ֵ����һά���飬��λK��
% �ռ䲽����step_L,Ĭ��0.01m��
% ʱ�䲽����step_t,Ĭ��1s��
% 
% ����б���������������������
% �����¶�Tout��һά���飬����ĳ���Ϊn��
% ����ѹ��Pout��һά���飬���鳤��Ϊn����λkPa��
% һ��ʱ�䲽��ı����¶�Tw_new��һά���顣
%%����refprop.m��3���ļ�,steel304.m,G10.m���������ļ�,��main��������
if nargin<8
    error("���������������8����")
end

if nargin<9|isempty(step_L)
    step_L=0.01;
end

if nargin<10|isempty(step_t)
    step_t=1; 
end

n=ceil(L/step_L); %�ֶܷ���������ȡ����
L_last=L-(n-1)*step_L; %���һ�εĳ��ȡ���ֹ�������������

if length(Tw)==1
    Tw(1:n)=Tw;  %���ֻ��һ����ֵ�������ȫ������¶�
end 

q_leak(1:n)=G10(Tw(1:n)/2).*(Te-Tw(1:n)); %��д��G10�ȵ��ʵĺ����ˡ�
% q_leak(1:n)=0.5*(Te-Tw(1:n));  %����©�������ܶȣ�����ĳ�����һ��ϵ������λ�Ͷ�������ϵ��һ�£�W/m2*K��

Tin(1:n)=Tin;
Pin(1:n)=Pin/1000; %Ԥ�����ڴ棬ͬʱ��ʼ��,ͬʱ��ѹǿ��λ���KPa������refprop����
Tout(1:n)=0;Pout(1:n)=0;
Tw_inner(1:n)=Tw(1:n);
QQ(1:n)=0; %ÿ����һ��ʱ�䲽���ڵ�������

for i=1:n-1
    %�����������ԣ�����ø�һ��refpropm.m���·������c��·��Ų����������
    %����Ϊ�˽�ʡ����ʱ�䣬���ó����ƽ��ֵ�������ˣ�ֱ����������㣬�����࣬ʡȥ�������̡�
    fluid1='air.ppf';
    fluid2='water';
    fluid3='helium';
   
    Cp=refpropm('C','T',Tin(i),'P',Pin(i),fluid3); %Cp [J/(kg K)]
    Den=refpropm('D','T',Tin(i),'P',Pin(i),fluid3); %Density [kg/m^3]
    Cond=refpropm('L','T',Tin(i),'P',Pin(i),fluid3); %Thermal conductivity [W/(m K)]
    Vis=refpropm('V','T',Tin(i),'P',Pin(i),fluid3); %Dynamic viscosity [Pa*s]
    
    Vel=(mg/Den)/(pi*(Din/2)*(Din/2)); % ���������ٶȣ�m/s
    Re=Vel*Den*Din/Vis; %������ŵ��
    Pr=Cp*Vis/Cond; %������������
    
%%%%%%%%%%%%%%����ѹ��%%%%%%%%%%%%%%%%%%%%%
    Rug=5e-5; %��ͨ���ڴֲڶȣ��׹��ȡ5e-7.
    XA=(2.457*log(1/((7/Re)^0.9+0.27*Rug/Din)))^16;%ϵ��A
    XB=(37530/Re)^16 ; %ϵ��B
    XF=8*((8/Re)^12+1/(XA+XB)^(3/2))^(1/12); %Ħ��ϵ��f
    Pout(i)=(Pin(i)*1000-0.5*XF*step_L/Din*Den*Vel^2)/1000; %����ѹ��ͬ������kpa��λ��������һ�ֵ��á�  
    Pin(i+1)=Pout(i);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if (Re>0)&&(Re<=2300)
        Nu=4;
    elseif Re>4000
        Nu=0.023*Re^0.8*Pr^0.4;
    elseif (Re>2300)&&(Re<=4000)
        Nu=(4^16+(0.023*Re^0.8*Pr^0.4)^16)^(1/16);
    else
        error("��ŵ������С��0�������������ݣ�")
    end
    
    h=Nu*Cond/Din; %���εĶ�������ϵ��
    [Cond_S,Cp_S,Den_S]=steel304(Tw(i)); %����304�����ܵ��ȵ��ʣ�����,�ܶȡ�
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

    AA=[A11,A12;A21,A22];   %�����������̣��þ������������Է����飬�õ��ܱ��ڱ��ºͱ��γ����¶ȡ�
    bb=[b1;b2];
    YY=AA\bb;
    Tw_inner(i)=YY(1);Tout(i)=YY(2);
    Tin(i+1)=Tout(i);
    QQ(i)=mg*Cp*(Tout(i)-Tin(i))-q_leak(i)*Dout*pi*step_L;
    
    %%%%%%%%%%%%%%����һ��ʱ�䲽��֮���������¶�%%%%%%%%%%%%
    Tw_new(i)=Tw(i)-QQ(i)*step_t/((Dout^2-Din^2)/4*pi*step_L*Cp_S*Den_S);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end


   %����ǰ��n-1�Σ������һ�Σ����һ����Ϊ�˷�ֹ�ܳ��Ȳ��������ֶγ��ȶ�Ԥ����һ�Ρ�     
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
%%%%%%%%%%%%%%����ѹ��%%%%%%%%%%%%%%%%%%%%%
    Rug=5e-5; %��ͨ���ڴֲڶȣ��׹��ȡ5e-7.
    XA=(2.457*log(1/((7/Re)^0.9+0.27*Rug/Din)))^16;
    XB=(37530/Re)^16 ; 
    XF=8*((8/Re)^12+1/(XA+XB)^(3/2))^(1/12); 
    Pout(n)=(Pin(n)*1000-0.5*XF*L_last/Din*Den*Vel^2);   %����Ͳ��ó���1000�ˣ������ĵ�λֱ����Pa��
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    if (Re>0)&&(Re<=2300)
        Nu=4;
    elseif Re>4000
        Nu=0.023*Re^0.8*Pr^0.4;
    elseif (Re>2300)&&(Re<=4000)
        Nu=(4^16+(0.023*Re^0.8*Pr^0.4)^16)^(1/16);
    else
        error("��ŵ������С��0�������������ݣ�")
    end
    
    h=Nu*Cond/Din; %���εĶ�������ϵ��
    [Cond_S,Cp_S]=steel304(Tw(n)); %����304�����ܵ��ȵ��ʣ����ݡ�
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
    
    AA=[A11,A12;A21,A22];   %�����������̣��þ������������Է����飬�õ��ܱ��ڱ��ºͱ��γ����¶ȡ�
    bb=[b1;b2];
    YY=AA\bb;
    Tw_inner(n)=YY(1);Tout(n)=YY(2);
    QQ(n)=mg*Cp*(Tout(n)-Tin(n))-q_leak(n)*Dout*pi*L_last;
    Tw_new(n)=Tw(n)-QQ(n)*step_t/((Dout^2-Din^2)/4*pi*L_last*Cp_S*Den_S);
end

