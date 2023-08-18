% 测试 myChan2 的定位效果

% 基站数目
BSN = 6;

% 各个基站的位置
BS = [0, 100, 150, 50, 200, 100; 
      -100,      0, 300, 190,   0,  300   ]; %基站位置
% BS = [0, sqrt(3), 0.5*sqrt(3), -0.5*sqrt(3), -sqrt(3), -0.5*sqrt(3), 0.5*sqrt(3); 
%       0,      0,         1.5,          1.5,        0,         -1.5,        -1.5]; 
% BS = BS(:,1:BSN);
% BS = BS .* 50;

% MS的实际位置
MS = [1, 1];

% R0i是各个BS与MS的实际距离，无噪声
for i = 1: BSN
    R0(i) = sqrt((BS(1,i) - MS(1))^2 + (BS(2,i) - MS(2))^2); 
end

% 噪声方差
Noise = 1;
numEpochs = 1000;
t = linspace(0,1,BSN);
errs=zeros(numEpochs,1);
for ep = 1:numEpochs
    R0_meas = R0 + Noise* randn(1,BSN);
%     R0_meas = R0;
    % R=R_{i,1},是加上了噪声后，BSi与BS1到MS的距离差，在实际使用中应该由 TDOA * c算得
    for i = 1: BSN-1
        R(i) = R0_meas(i+1) - R0_meas(1) ; 
    end
    X = hchan(BSN, BS, R);
    err = sqrt(sum((X.' - MS).^2));
    errs(ep) = err;
end
RMSerr = sqrt(mean(errs.^2));
X


