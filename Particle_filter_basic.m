%% Passive target tracking by particle filter

clc;
clear all;
close all;

T=1; % Time interval
F=[1 T 0 0;0 1 0 0;0 0 1 T;0 0 0 1]; % State transition matrix
C=[T^2/2 0;T 0;0 T^2/2;0 T]; % Process noise distribution matrix
H=[1 0 0 0;0 0 1 0]; % Measurement transition matrix

x=20; % Initial position on x
y=20; % Initial position on y
vx=2; % Initial velocity on x
vy=1; % Initial position on y

q1=1; % Process noise on x
q2=1; % Process noise on y

r_angle=1; % Measurement noise ( 1бу = 1*pi/180 rad )

len = 200; % Total tracking steps

% Sitimulate the real target track on x-y
for k=1:len
    ZS(1,k)=x+vx*k*T;
    ZS(2,k)=vx;
    ZS(3,k)=y+vy*k*T;
    ZS(4,k)=vy;
    RS(k)=sqrt(ZS(1,k)^2+ZS(3,k));
    THS(k)=atand(ZS(3,k)/ZS(1,k));
end

% Parameters for particle filter
N = 100;   % Total particle number

% Tracking world size
WorldSizex = 250;   % x range
WorldSizey = 250;   % y range
X = zeros(2, len);  % States
P = zeros(2, N);    % Particle sets

PCenter = zeros(1, len);  % Center of the particles
w = zeros(N, 1);          % Particle weights
err = zeros(1,len);       % Errors
RTH=zeros(1,len);         % Measurements (bearing)

X(:,1)=[x;y];             % Real initial statement  
XTH(1)=atand(X(2,1)./X(1,1));
W=normrnd(0,r_angle,1,1);
RTH(1)=atand(ZS(3,1)/ZS(1,1))+W;

time=100; % Monte Carlo simulating times

for j=1:time
    % Initialization
    for i = 1 : N
        P(:, i) = [WorldSizex*rand;WorldSizey*rand];
        P1(i)=atand( P(2, i)/ P(1, i));
        dist = abs(P1(i)-RTH(1));   
        w(i) = (1 / sqrt(r_angle) / sqrt(2 * pi)) * exp(-(dist)^2 / 2 / r_angle);  
    end
    PCenter(1) = sum(P1) / N;     
    err(j,1)=abs(PCenter(1)-RTH(1));
    
    % Start to track
    for k = 2:len
        v1=normrnd(0,q1,1,1);               
        v2=normrnd(0,q2,1,1);               
        X(:,k)=X(:,k-1)+T*[vx;vy]+[v1;v2];  
        XTH(k)=atand(X(2,k)/X(1,k));        
        XR(k)=sqrt(X(2,k)^2+X(1,k)^2);       
        W=normrnd(0,r_angle,1,1);
        RTH(k)=XTH(k)+W;                   

        for i = 1:N
            v1=normrnd(0,q1,1,1);               
            v2=normrnd(0,q2,1,1);               
            P(:, i)=P(:, i)+T*[vx;vy]+[v1;v2];
            P1(i)=atand(P(2, i)/P(1, i));
            dist = abs(P1(i)-RTH(k));
            w(i) = (1 / sqrt(r_angle) / sqrt(2 * pi)) * exp(-(dist)^2 / 2 / r_angle);   % Calculate weights
        end
        
        % Normalization weights
        wsum = sum(w);
        for i = 1 : N
            w(i) = w(i) / wsum;
        end
        
        % Resampling
        for i = 1 : N
            wmax = 2 * max(w) * rand;  
            index = randi(N, 1); 
            while(wmax > w(index))
                wmax = wmax - w(index);
                index = index + 1;
                if index > N
                    index = 1;
                end
            end
            P(:, i) = P(:, index);   
            
        end
        
        PCenter(k) = sum(P1) / N;   % Calculte the central of all particles as the estimating result
                
        %     figure(10);
        %     set(gca,'FontSize',12);
        %     clf;
        %     hold on
        %     plot(RTH(1, k), RTH(2, k), 'r.', 'markersize',50);  
        %     axis([0 WorldSizex 0 WorldSizey]);
        %     grid on
        %     plot(P(1, :), P(2, :), 'k.', 'markersize',5); 
        %     plot(PCenter(1, k), PCenter(2, k), 'b.', 'markersize',25);
        %     legend('True State', 'Particle', 'The Center of Particles');
        %     hold off
        %     pause(0.1);
    end
    
    % Error calculation
    err(j,:)=abs(PCenter-XTH);
    if j==1
        figure;  % Cartesian coordinates
        set(gca,'FontSize',12);
        plot(1:len, XTH, 'r*-',1:len, RTH , 'g.-',1:len,PCenter , 'b-');
        legend('True States', 'Measurements', 'Particle filter results');
        xlabel('t/s'); ylabel('ж╚/бу');
        grid on
        
    end
end
err1=mean(err);
%%
figure;  % Error plotting
set(gca,'FontSize',12);
plot(err1,'.-');
xlabel('t/s');
ylabel('RMSE/m');
legend('Particle Filter')
title('The err');
grid on