% 构建三角窗函数
% 输入：自变量x
% 输入：因变量y
function y= TriangleWindow(x)
    % 构建与自变量同等维度的y
    y = zeros(length(x),1);
    % 定义三角窗函数
    index = find( x >= -1 & x <= 1);
    y(index) = 1 - abs(x(index));
end

