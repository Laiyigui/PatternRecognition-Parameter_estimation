% �������Ǵ�����
% ���룺�Ա���x
% ���룺�����y
function y= TriangleWindow(x)
    % �������Ա���ͬ��ά�ȵ�y
    y = zeros(length(x),1);
    % �������Ǵ�����
    index = find( x >= -1 & x <= 1);
    y(index) = 1 - abs(x(index));
end

