function [lambda] = G10(Temperature)
%%G10�ĵ���ϵ����������������pipe1d���©��ϵ��
TT=log10(Temperature);
a=-4.1236;
b=13.788;
c=-26.068;
d=26.272;
e=-14.663;
f=4.4954;
g=-0.6905;
h=0.0397;
i=0;
lambda=10.^(a+b*TT+c*TT.^2+d*TT.^3+e*TT.^4+f*TT.^5+g*TT.^6+h*TT.^7+i*TT.^8);
end