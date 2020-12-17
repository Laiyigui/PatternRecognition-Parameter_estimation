% �ô������ڻ��Ʋ�ͬ�������µĹ��Ƹ��ʺ�������ͼ
%% ���������Ƹ����ܶȺ���ֵ����������Ӧ����
clear all;clc;
% ����p(x)=0.2*N(-1,1)+0.8*N(1,1)��������ͼ���Թ��Ƚ�
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
title('�����Ƹ��ʺ���ͼ��');
saveas(gcf,'�����Ƹ��ʺ���ͼ��.jpg')

%% ������ʼ��
%����������N���ɵ��ڲ���h1��ȡֵ
N_set = [5,10,50,100,1000,10000];
h1_set = [0.25,0.5,1,2,4];
%subplot��ͼ�ĺᡢ�������
image_col=size(h1_set,2);%����
image_row=size(N_set,2);%����
%�����ƺ���pN(x)���Ա���ȡֵ�Լ��ڸ��Ա����º���ȡֵ������ֵ���ڻ���ͼ��
x = linspace(-5,5);
y = 0.2*1/sqrt(2*pi)*exp(-1*(x+1).^2/2)+0.8*1/sqrt(2*pi)*exp(-1*(x-1).^2/2);

%% ����
figure(2);
for i = 1:length(N_set)
    N = N_set(i);
    for j = 1:length(h1_set)
        hN = h1_set(j);
        px = zeros(length(x),1);
        index = randperm(10000,N);%10000Ϊ���������ܶȺ������Ա�������
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
suptitle('��������')
saveas(gcf,'�������Ƹ��ʺ���ͼ��.jpg')

%% ��˹������
figure(3);
for i = 1:length(N_set)
    N = N_set(i);
    for j = 1:length(h1_set)
        h1 = h1_set(j);
        hN = h1/sqrt(N);
        px = zeros(length(x),1);
        index = randperm(10000,N);%10000Ϊ���������ܶȺ������Ա�������
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
suptitle('��˹������')
saveas(gcf,'��˹�����Ƹ��ʺ���ͼ��.jpg')

%% ָ��������
figure(4);
for i = 1:length(N_set)
    N = N_set(i);
    for j = 1:length(h1_set)
        h1 = h1_set(j);
        hN = h1/sqrt(N);
        px = zeros(length(x),1);
        index = randperm(10000,N);%10000Ϊ���������ܶȺ������Ա�������
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
suptitle('ָ��������')
saveas(gcf,'ָ�������Ƹ��ʺ���ͼ��.jpg')

%% ���Ǵ�����
figure(5);
for i = 1:length(N_set)
    N = N_set(i);
    for j = 1:length(h1_set)
        hN = h1_set(j);
        px = zeros(length(x),1);
        index = randperm(10000,N);%10000Ϊ���������ܶȺ������Ա�������
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
suptitle('���Ǵ�����')
saveas(gcf,'���Ǵ����Ƹ��ʺ���ͼ��.jpg')