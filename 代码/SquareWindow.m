% ������������
% ���룺�Ա���x
% ���룺�����y
function y= SquareWindow(x)
    % �������Ա���ͬ��ά�ȵ�y
    y = zeros(length(x),1);
    % ���巽������
    index = find(x >= -1/2 & x <= 1/2);
    y(index) = 1;
end

