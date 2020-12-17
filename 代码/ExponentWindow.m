% 构建指数窗函数
% 输入：自变量x
% 输入：因变量y
function y = ExponentWindow(x)
    y = exp(-abs(x))/2;
end
