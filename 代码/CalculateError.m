% �ô������ڼ�����Ƹ����ܶȺ���ֵ���׼�����ܶȺ���ֵ֮��Ĳ�ֵ���������ֵ
%% ���������Ƹ����ܶȺ���ֵ
clear all;clc;
% ����p(x)=0.2*N(-1,1)+0.8*N(1,1)��������ͼ���Թ��Ƚ�
x1 = normrnd(-1,1,1,2000);
x2 = normrnd(1,1,1,8000);
xi = [x1 x2];
yi = 0.2*1/sqrt(2*pi)*exp(-1*(xi+1).^2/2)+0.8*1/sqrt(2*pi)*exp(-1*(xi-1).^2/2);

%% ������ʼ��
%����������N���ɵ��ڲ���h1��ȡֵ
N_set = [5,10,50,100,1000,10000];
h1_set = [0.25,0.5,1,2,4];

%�����ƺ���pN(x)���Ա���ȡֵ�Լ��ڸ��Ա����º���ȡֵ������ֵ���ڻ���ͼ��
dx = 0.02;
x = -5:dx:5;
y = 0.2*1/sqrt(2*pi)*exp(-1*(x+1).^2/2)+0.8*1/sqrt(2*pi)*exp(-1*(x-1).^2/2);

%���ڶ�μ����ֵ�Ĵ���
times = 20;

%% ����
%���ڴ���e��ֵ
e = zeros(times,1);
e_Square_E = zeros(length(N_set),length(h1_set));
e_Square_Var = zeros(length(N_set),length(h1_set));
for i = 1:length(N_set)
    N = N_set(i);
    for j = 1:length(h1_set)
        hN = h1_set(j);
        px = zeros(length(x),1);
        index = randperm(10000,N);%10000Ϊ���������ܶȺ������Ա�������
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

%% ��˹��
%���ڴ���e��ֵ
e = zeros(times,1);
e_Gaussian_E = zeros(length(N_set),length(h1_set));
e_Gaussian_Var = zeros(length(N_set),length(h1_set));
for i = 1:length(N_set)
    N = N_set(i);
    for j = 1:length(h1_set)
        h1 = h1_set(j);
        hN = h1/sqrt(N);
        px = zeros(length(x),1);
        index = randperm(10000,N);%10000Ϊ���������ܶȺ������Ա�������
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

%% ָ����
%���ڴ���e��ֵ
e = zeros(times,1);
e_Exponent_E = zeros(length(N_set),length(h1_set));
e_Exponent_Var = zeros(length(N_set),length(h1_set));
for i = 1:length(N_set)
    N = N_set(i);
    for j = 1:length(h1_set)
        h1 = h1_set(j);
        hN = h1/sqrt(N);
        px = zeros(length(x),1);
        index = randperm(10000,N);%10000Ϊ���������ܶȺ������Ա�������
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

%% ���Ǵ�
%���ڴ���e��ֵ
e = zeros(times,1);
e_Triangle_E = zeros(length(N_set),length(h1_set));
e_Triangle_Var = zeros(length(N_set),length(h1_set));
for i = 1:length(N_set)
    N = N_set(i);
    for j = 1:length(h1_set)
        hN = h1_set(j);
        px = zeros(length(x),1);
        index = randperm(10000,N);%10000Ϊ���������ܶȺ������Ա�������
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