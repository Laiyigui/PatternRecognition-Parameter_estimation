% 构建方窗函数
% 输入：自变量x
% 输入：因变量y
function y= SquareWindow(x)
    % 构建与自变量同等维度的y
    y = zeros(length(x),1);
    % 定义方窗函数
    index = find(x >= -1/2 & x <= 1/2);
    y(index) = 1;
end

