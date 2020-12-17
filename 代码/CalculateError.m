% 该代码用于计算估计概率密度函数值与标准概率密度函数值之间的差值、方差、期望值
%% 给出待估计概率密度函数值
clear all;clc;
% 函数p(x)=0.2*N(-1,1)+0.8*N(1,1)，并绘制图像，以供比较
x1 = normrnd(-1,1,1,2000);
x2 = normrnd(1,1,1,8000);
xi = [x1 x2];
yi = 0.2*1/sqrt(2*pi)*exp(-1*(xi+1).^2/2)+0.8*1/sqrt(2*pi)*exp(-1*(xi-1).^2/2);

%% 参数初始化
%给出样本量N，可调节参量h1的取值
N_set = [5,10,50,100,1000,10000];
h1_set = [0.25,0.5,1,2,4];

%待估计函数pN(x)的自变量取值以及在该自变量下函数取值，将该值用于绘制图线
dx = 0.02;
x = -5:dx:5;
y = 0.2*1/sqrt(2*pi)*exp(-1*(x+1).^2/2)+0.8*1/sqrt(2*pi)*exp(-1*(x-1).^2/2);

%用于多次计算差值的次数
times = 20;

%% 方窗
%用于储存e的值
e = zeros(times,1);
e_Square_E = zeros(length(N_set),length(h1_set));
e_Square_Var = zeros(length(N_set),length(h1_set));
for i = 1:length(N_set)
    N = N_set(i);
    for j = 1:length(h1_set)
        hN = h1_set(j);
        px = zeros(length(x),1);
        index = randperm(10000,N);%10000为给定概率密度函数的自变量个数
        for t = 1:times
            for k = 1:length(x)
                px(k) = sum(SquareWindow((x(k)-xi(index))/hN),1)/N;
            end
            e(t) = sum((px-y').^2,1)*dx;
        end
        e_Square_E(i,j) = mean(e);
        e_Square_Var(i,j) = var(e);
    end
end

%% 高斯窗
%用于储存e的值
e = zeros(times,1);
e_Gaussian_E = zeros(length(N_set),length(h1_set));
e_Gaussian_Var = zeros(length(N_set),length(h1_set));
for i = 1:length(N_set)
    N = N_set(i);
    for j = 1:length(h1_set)
        h1 = h1_set(j);
        hN = h1/sqrt(N);
        px = zeros(length(x),1);
        index = randperm(10000,N);%10000为给定概率密度函数的自变量个数
        for t = 1:times
            for k = 1:length(x)
                px(k) = sum(GaussianWindow((x(k)-xi(index))/hN),2)/(N*hN);
            end
            e(t) = sum((px-y').^2,1)*dx;
        end
        e_Gaussian_E(i,j) = mean(e);
        e_Gaussian_Var(i,j) = var(e);
    end
end

%% 指数窗
%用于储存e的值
e = zeros(times,1);
e_Exponent_E = zeros(length(N_set),length(h1_set));
e_Exponent_Var = zeros(length(N_set),length(h1_set));
for i = 1:length(N_set)
    N = N_set(i);
    for j = 1:length(h1_set)
        h1 = h1_set(j);
        hN = h1/sqrt(N);
        px = zeros(length(x),1);
        index = randperm(10000,N);%10000为给定概率密度函数的自变量个数
        for t = 1:times
            for k = 1:length(x)
                px(k) = sum(ExponentWindow((x(k)-xi(index))/hN),2)/(N*hN);
            end
            e(t) = sum((px-y').^2,1)*dx;
        end
        e_Exponent_E(i,j) = mean(e);
        e_Exponent_Var(i,j) = var(e);
    end
end

%% 三角窗
%用于储存e的值
e = zeros(times,1);
e_Triangle_E = zeros(length(N_set),length(h1_set));
e_Triangle_Var = zeros(length(N_set),length(h1_set));
for i = 1:length(N_set)
    N = N_set(i);
    for j = 1:length(h1_set)
        hN = h1_set(j);
        px = zeros(length(x),1);
        index = randperm(10000,N);%10000为给定概率密度函数的自变量个数
        for t = 1:times
            for k = 1:length(x)
                px(k) = sum(TriangleWindow((x(k)-xi(index))/hN),1)/N;
            end
            e(t) = sum((px-y').^2,1)*dx;
        end
        e_Triangle_E(i,j) = mean(e);
        e_Triangle_Var(i,j) = var(e);
    end
end