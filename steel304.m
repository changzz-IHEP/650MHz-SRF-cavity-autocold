function [lambda,Cp,Den] = steel304(Temperature)
%steel304 ���ڼ���304L����������¶ȵ����ݱ仯���ȵ��ʱ仯��
%   �����¶ȣ�����ȵ��ʺͱ�����,ֻ�ܽ��ձ�����
if Temperature<=300&Temperature>6
    lambda=-0.71544+0.16077*Temperature-(7.02864e-4)*Temperature^2+8.61764e-7*Temperature^3+3.31512e-9*Temperature^4-8.04863e-12*Temperature^5;
elseif  Temperature>4&Temperature<=6
    lambda=0.1;
elseif  Temperature>300
    lambda=22.9378-0.06563*Temperature+(1.56186e-4)*Temperature^2-(1e-7)*Temperature^3; %ʵ��ֻ���õ�700K������700K�Ͳ�׼�ˡ�
else
    error('�����¶Ȳ���С��4��')
end

if Temperature<=300&Temperature>4
    TT=log10(Temperature);
    Cp=10^(22.0061-127.5528*TT+303.647*TT^2-381.0098*TT^3+274.0328*TT^4-112.9212*TT^5+24.7593*TT^6-2.239153*TT^7);
elseif Temperature>300
    Cp=210.26181+1.52402*Temperature-0.00266*Temperature^2+1.66667e-6*Temperature^3; %ʵ��ֻ���õ�700K�����������Χ�û�����Ϲ�ʽ��
else
    error('�����¶Ȳ���С��4��')
end

Den=7900; %304���ܶ�
end