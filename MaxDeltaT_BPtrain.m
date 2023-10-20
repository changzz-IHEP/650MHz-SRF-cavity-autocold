clc;clear;close all;
format long
data=importdata('CalData_BTCM_deltaT.txt');
%%
%生成随机�?
LL=size(data,1); %数据的�?�个�?
k=rand(1,LL);
[m,n]=sort(k);%n的大小与k相同，描述了k的元素沿排序维数排列成m的情况�?�即返回排序后元素的索引
%%
input(:,1)=data(:,1);
input(:,2)=data(:,3);
output(:,1)=data(:,2); %入口温度是输�?
%%
%随机提取80%个样本为训练样本�?20%个样本为预测样本
LL_train=ceil(LL*0.8);
LL_test=floor(LL*0.2);
input_train=input(n(1:LL_train),:)';
output_train=output(n(1:LL_train),:)';
input_test=input(n((LL_train+1):LL),:)';
output_test=output(n((LL_train+1):LL),:)';

%%
%归一�?
[inputn,inputps]=mapminmax(input_train);
[outputn,outputps]=mapminmax(output_train);
%% BP网络训练
% 初始化网络结�?
net=newff(inputn,outputn,5);

net.trainParam.epochs=100;
net.trainParam.lr=0.1;
net.trainParam.goal=0.00004;

%% 网络训练
net=train(net,inputn,outputn);

%% BP网络预测
% 预测数据归一�?
inputn_test=mapminmax('apply',input_test,inputps); %apply表示按照训练数据的归�?化方式进行归�?�?
 
%% 网络预测输出
an=sim(net,inputn_test);
 
%% 网络输出反归�?�?
BPoutput=mapminmax('reverse',an,outputps); %按照训练数据的归�?化方式进行反归一化�??

%%
% 单个测试
input_test2=[290;12]; %单个数据
inputn_test2=mapminmax('apply',input_test2,inputps); 
test2=sim(net,inputn_test2);
output_test2=mapminmax('reverse',test2,outputps)

%%
% BP结果分析
figure(1)
plot(BPoutput,':og')
hold on
plot(output_test,'-*');
legend('预测输出','期望输出')
title('BP网络预测输出','fontsize',12)
ylabel('函数输出','fontsize',12)
xlabel('样本','fontsize',12)
%%
% 预测误差
error=(BPoutput-output_test)./output_test;
figure(2)
plot(error,'.')
ylabel('预测温度的相对误�?','fontsize',12)
% set(gca,'yticklabel',{'-50%', '-25%','0','25%', '50%' ,'75%','100%'});
xlabel('样本','fontsize',12)
%%
save ('BTCM_BP_3D_net1.mat','net');
save ('BTCM_BP_3D_net1_ioPS.mat','inputps','outputps');%输出对应的inputps和outputps