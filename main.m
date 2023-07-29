clc,clear;
close all;

%% parameters

N = 10;%ULA 阵元数目
M = 40;%信号样本数目
lambda = 1;
d = lambda/2;
SNR=[10;20];sigma=1;% T =1;


Am=sqrt(2*sigma^2*10.^(SNR/10));


%信号方案一
x1=Am(1)*ones(1,M);
x2=Am(2)*exp(1i*2*pi*rand(1,M));

Theta = [-10 ,30]/180*pi;
n = 0:N-1;
A = exp(-1i*2*pi*d/lambda*n'*sin(Theta));

s = A*[x1;x2];

x = awgn(s,10,'measured');%加入高斯白噪声

%计算协方差矩阵

R = zeros(N,N);
for idx = 1:length(x1)
R = R+x(:,idx)*x(:,idx)';
end
R = R/length(x1);
[V,D] = eig(R);
eig_value = real(diag(D));
[B,I] = sort(eig_value,'descend');

Vspace = V(:,I(3:end));%构造噪声子空间


theta = -pi/2:pi/400:pi/2;
P = zeros(1,length(theta));
for i=1:length(theta)
   a= exp(-1i*2*pi*d/lambda*n'*sin(theta(i)));
   P(i) = 1/(a'*(Vspace*Vspace')*a);

end
mag = db(abs(P)/max(abs(P)));
ymin = min(mag);
ymax = max(mag);
hold on

plot(theta*180/pi,mag,'linewidth',1.2);
plot([Theta(1)*180/pi,Theta(1)*180/pi],[ ymin ymax ],'r','linewidth',1);%角度1
plot([Theta(2)*180/pi,Theta(2)*180/pi],[ymin ymax],'r','linewidth',1);%角度2
grid minor
xlabel('\theta');
ylabel('amplitude (dB)');




%%%%%%% RMSE VS SNR %%%%%%
clc,clear;
close all;


N = 10;%ULA 阵元数目
Num_snap=100;%snapshot
w=pi/4;%角频率
lambda=2*pi*3e8/w;%波长lamda
d=0.5*lambda;%阵元间距

SNR=-10:1:20;%SNR(dB)
RMSE1=zeros(1,length(SNR));

Theta = -10/180*pi;
n = 0:N-1;
A = exp(-1i*2*pi*d/lambda*n'*sin(Theta));
source_number=1;

theta1=[-90:90]*pi/180;
M=1000;%蒙特卡洛仿真次数

RMSE=zeros(1,length(SNR));
CRb=zeros(1,length(SNR));

for cnt=1:length(SNR)
    SErr=0;

    for ii=1:M
        x1=exp(1i*w*[0:Num_snap-1]);
        s = A*x1;
% x =s;
        x = awgn(s,SNR(cnt),'measured');

        R=x*x'/Num_snap;
        [V,D] = eig(R);
        D=diag(D);
        Vspace=V(:,1:N-source_number);%噪声部分所对应的特征向量构成的噪声子空间


        i = 1;
        for kk=1:length(theta1)
           a= exp(-1i*2*pi*d/lambda*n'*sin(theta1(kk)));
           P(i) = 1/abs(a'*(Vspace*Vspace')*a);
           i = i+1;
        end
        [~,ind]=max(abs(P));
        theta_hat=theta1(ind);
        SErr=SErr+(theta_hat-Theta).^2;


    end
    RMSE(cnt)=sqrt(SErr/M);
    cnt

end
figure
plot(SNR,RMSE*180/pi,'b-o','linewidth',1);

xlabel('SNR (dB)');
ylabel('RMSE (degree)');
grid on;
