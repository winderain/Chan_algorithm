% ���� myChan2 �Ķ�λЧ��

% ��վ��Ŀ
BSN = 6;

% ������վ��λ��
BS = [0, 100, 150, 50, 200, 100; 
      -100,      0, 300, 190,   0,  300   ]; %��վλ��
% BS = [0, sqrt(3), 0.5*sqrt(3), -0.5*sqrt(3), -sqrt(3), -0.5*sqrt(3), 0.5*sqrt(3); 
%       0,      0,         1.5,          1.5,        0,         -1.5,        -1.5]; 
% BS = BS(:,1:BSN);
% BS = BS .* 50;

% MS��ʵ��λ��
MS = [1, 1];

% R0i�Ǹ���BS��MS��ʵ�ʾ��룬������
for i = 1: BSN
    R0(i) = sqrt((BS(1,i) - MS(1))^2 + (BS(2,i) - MS(2))^2); 
end

% ��������
Noise = 1;
numEpochs = 1000;
t = linspace(0,1,BSN);
errs=zeros(numEpochs,1);
for ep = 1:numEpochs
    R0_meas = R0 + Noise* randn(1,BSN);
%     R0_meas = R0;
    % R=R_{i,1},�Ǽ�����������BSi��BS1��MS�ľ�����ʵ��ʹ����Ӧ���� TDOA * c���
    for i = 1: BSN-1
        R(i) = R0_meas(i+1) - R0_meas(1) ; 
    end
    X = hchan(BSN, BS, R);
    err = sqrt(sum((X.' - MS).^2));
    errs(ep) = err;
end
RMSerr = sqrt(mean(errs.^2));
X


