%%testblackbody.m
%%for B-v
%%for matlab


clear;clc; clf;
c1=3.741832; % 第一辐射常数
c2=14387.86; % 第二辐射常数
h=6.63*10^(-34);
c=3*10^8;
k=1.38*10^(-23);
figure();
for t=-4:1:0 % 设置辐射温度(K)
T=10^t;
f=logspace(7,13); % 设置频率范围范围及计算步长
M=(2*h*f.^3/c/c)./(exp(h*f./k/T)-1); % 计算指定温度光谱辐出度
loglog(f,M,'-b','LineWidth',1.4) % 绘制光谱辐出度曲线
maxM=max(M) % 找出指定温度最大光谱辐出度
i=find(maxM==M); % 找峰值波长点
hold on % 在指定位置按给定方式标记对应温度
end

xlabel('f') % 横坐标名称及单位
ylabel('I') % 纵坐标名称及单位
hold off

%%for B-l
clear;clc; clf;
c1=3.741832; % 第一辐射常数
c2=14387.86; % 第二辐射常数
figure();
for T=800:200:1500 % 设置辐射温度(K)
%T=10^t;
l=0.0001:0.02:15; % 设置波长范围及计算步长
M=1e+4.*(c1./(l.^5)./(exp(c2./(l.*T))-1)); % 计算指定温度光谱辐出度
plot(l,M,'-b','LineWidth',1.4) % 绘制光谱辐出度曲线
maxM=max(M) % 找出指定温度最大光谱辐出度
i=find(maxM==M); % 找峰值波长点
text(l(i+20),M(i+20),[num2str(T),'K'],'VerticalAlignment',...
'baseline','HorizontalAlignment','left','fontsize',7)
hold on % 在指定位置按给定方式标记对应温度
end
set(gca,'XTick',[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]) %设置横坐标点
xlabel('\lambda / \mum') % 横坐标名称及单位
ylabel('M_{b\lambda} / W\cdotcm^{-2}\cdot\mum^{-1}') % 纵坐标名称及单位
hold off
