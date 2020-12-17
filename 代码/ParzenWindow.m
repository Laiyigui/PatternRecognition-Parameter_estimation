% 该代码用于绘制不同窗函数下的估计概率函数曲线图
%% 给出待估计概率密度函数值，并绘制相应曲线
clear all;clc;
% 函数p(x)=0.2*N(-1,1)+0.8*N(1,1)，并绘制图像，以供比较
x1 = normrnd(-1,1,1,2000);
x2 = normrnd(1,1,1,8000);
xi = [x1 x2];
yi = 0.2*1/sqrt(2*pi)*exp(-1*(xi+1).^2/2)+0.8*1/sqrt(2*pi)*exp(-1*(xi-1).^2/2);
figure(1);
plot(xi,yi,'b.');
ylim([0 0.4]);
xlabel('xi');
xlabel('yi');
grid on; 
title('待估计概率函数图形');
saveas(gcf,'待估计概率函数图形.jpg')

%% 参数初始化
%给出样本量N，可调节参量h1的取值
N_set = [5,10,50,100,1000,10000];
h1_set = [0.25,0.5,1,2,4];
%subplot子图的横、纵向个数
image_col=size(h1_set,2);%纵向
image_row=size(N_set,2);%横向
%待估计函数pN(x)的自变量取值以及在该自变量下函数取值，将该值用于绘制图线
x = linspace(-5,5);
y = 0.2*1/sqrt(2*pi)*exp(-1*(x+1).^2/2)+0.8*1/sqrt(2*pi)*exp(-1*(x-1).^2/2);

%% 方窗
figure(2);
for i = 1:length(N_set)
    N = N_set(i);
    for j = 1:length(h1_set)
        hN = h1_set(j);
        px = zeros(length(x),1);
        index = randperm(10000,N);%10000为给定概率密度函数的自变量个数
        for k = 1:length(x)
            px(k) = sum(SquareWindow((x(k)-xi(index))/hN),1)/N;
        end
        subplot(image_row,image_col, image_col*(i-1)+j);
        plot(x,px,'r');
        hold on
        plot(x,y,'b');
        if i == 1
            line = [ 'h1=',num2str(hN)];
            title(line);
        end
        if j == 1
            yline = [ 'N=',num2str(N)];
            ylabel(yline);
        end
        grid on;       
    end
end
suptitle('方窗估计')
saveas(gcf,'方窗估计概率函数图形.jpg')

%% 高斯窗估计
figure(3);
for i = 1:length(N_set)
    N = N_set(i);
    for j = 1:length(h1_set)
        h1 = h1_set(j);
        hN = h1/sqrt(N);
        px = zeros(length(x),1);
        index = randperm(10000,N);%10000为给定概率密度函数的自变量个数
        for k = 1:length(x)
            px(k) = sum(GaussianWindow((x(k)-xi(index))/hN),2)/(N*hN);
        end
        subplot(image_row,image_col, image_col*(i-1)+j);
        plot(x,px,'r');
        hold on
        plot(x,y,'b');
        if i == 1
            line = [ 'h1=',num2str(h1)];
            title(line);
        end
        if j == 1
            yline = [ 'N=',num2str(N)];
            ylabel(yline);
        end
        grid on;       
    end
end
suptitle('高斯窗估计')
saveas(gcf,'高斯窗估计概率函数图形.jpg')

%% 指数窗估计
figure(4);
for i = 1:length(N_set)
    N = N_set(i);
    for j = 1:length(h1_set)
        h1 = h1_set(j);
        hN = h1/sqrt(N);
        px = zeros(length(x),1);
        index = randperm(10000,N);%10000为给定概率密度函数的自变量个数
        for k = 1:length(x)
            px(k) = sum(ExponentWindow((x(k)-xi(index))/hN),2)/(N*hN);
        end
        subplot(image_row,image_col, image_col*(i-1)+j);
        plot(x,px,'r');
        hold on
        plot(x,y,'b');
        if i == 1
            line = [ 'h1=',num2str(h1)];
            title(line);
        end
        if j == 1
            yline = [ 'N=',num2str(N)];
            ylabel(yline);
        end
        grid on;       
    end
end
suptitle('指数窗估计')
saveas(gcf,'指数窗估计概率函数图形.jpg')

%% 三角窗估计
figure(5);
for i = 1:length(N_set)
    N = N_set(i);
    for j = 1:length(h1_set)
        hN = h1_set(j);
        px = zeros(length(x),1);
        index = randperm(10000,N);%10000为给定概率密度函数的自变量个数
        for k = 1:length(x)
            px(k) = sum(TriangleWindow((x(k)-xi(index))/hN),1)/N;
        end
        subplot(image_row,image_col, image_col*(i-1)+j);
        plot(x,px,'r');
        hold on
        plot(x,y,'b');
        if i == 1
            line = [ 'h1=',num2str(hN)];
            title(line);
        end
        if j == 1
            yline = [ 'N=',num2str(N)];
            ylabel(yline);
        end
        grid on;       
    end
end
suptitle('三角窗估计')
saveas(gcf,'三角窗估计概率函数图形.jpg')