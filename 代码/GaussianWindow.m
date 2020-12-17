% 构建高斯窗函数
% 输入：自变量x
% 输入：因变量y
function y = GaussianWindow(x)
    y = exp(-1/2*x.^2)/sqrt(2*pi);
end

