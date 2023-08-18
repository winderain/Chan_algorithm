function [X] = hchan(BSN,BS,R)

% 第一次WLS

 %k=X^2+Y^2
for i = 1:BSN                      %BSN为基站个数
    k(i) = BS(1,i)^2 + BS(2,i)^2;  %BS为基站坐标
end

%h = 1/2(Ri^2-ki+k1)
for i =1:BSN-1
    h(i) = 0.5*(R(i)^2 - k(i+1) + k(1));  %注意k(i+1)
end


%Ga = [Xi,Yi,Ri]
for i = 1:BSN-1
    Ga(i,1) = -BS(1,i+1)+BS(1,1);
    Ga(i,2) = -BS(2,i+1)+BS(2,1);
    Ga(i,3) = -R(i);
end

%Q为TDOA系统的协方差矩阵
Q = cov(R);

%MS与BS距离较远时
za = inv(Ga' * inv(Q) * Ga) * Ga' * inv(Q) * h';

%第二次WLS
%h'
X1 = BS(1,1);
Y1 = BS(2,1);
h2 = [
    (za(1,1) - X1)^2;
    (za(2,1) - Y1)^2;
     za(3,1)^2
      ];

%Ga'
Ga2 = [1,0;0,1;1,1];

%B'
B2 = [
      za(1,1)-X1,0,0;
      0,za(2,1)-Y1,0;
      0,0,za(3,1)
      ];
  
%za',距离较远时
% za2 = inv( Ga2' * inv(B2) * Ga' * inv(Q) * Ga * inv(B2) * Ga2) * (Ga2' * inv(B2) * Ga' * inv(Q) * Ga * inv(B2)) * h2;
za3 =  (   ( Ga2' / B2)   * Ga' / Q * Ga / B2 * Ga2) \ (  (Ga2' / B2 * Ga' / Q * Ga / B2 )  * h2); %更改了inv，变成/和\

za4 = zeros(2,1);
za4(1,1) = abs(sqrt(za3(1,1)));
za4(2,1) = abs(sqrt(za3(2,1)));

zp(1,1) = za4(1,1) + X1;
zp(2,1) = za4(2,1) + Y1;


X = zp;