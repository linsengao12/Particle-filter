%% Passive target tracking by PF (particle filter) - PDA (probability density assiociation)

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

r_angle=5; % Measurement noise ( 1бу = 1*pi/180 rad )

lamda=0.004; % Density of spurious measurement
Pd=0.9;      % Decetion probability

len=300;     % Total tracking steps
ZS(1,1)=x;
ZS(3,1)=y;
ZS(2,:)=vx;
ZS(4,:)=vy;

% Sitimulate the real target track on x-y
for k=2:len
    ZS(1,k)=ZS(1,k-1)+vx*T;
    ZS(3,k)=ZS(3,k-1)+vy*T;   
    RS(k)=sqrt(ZS(1,k)^2+ZS(3,k));
    THS(k)=atand(ZS(3,k)/ZS(1,k));
end

% Parameters for particle filter
N = 100;   % Total particle number

% Tracking world size
WorldSizex=180; 
WorldSizey=180; 
X=zeros(2,len); % States
P=zeros(2,N);   % Particle sets
PCenter=zeros(1,len);  % Center of the particles
w=zeros(N,1);          % Particle weights
err=zeros(1,len);      % Errors
nc=5;                  % Measurements number
RTH=zeros(len,nc);     % Measurements (bearing)
X(:,1)=[x;y];          % Real initial statement
XTH(1)=atand(X(2,1)./X(1,1));
W=normrnd(0,r_angle,1,1);
RTH(1,1:5)=XTH(1)+W;      

time=100;              % Monte Carlo simulating times

for m=1:time
    % Initialization
    for i=1:N
        P(:,i)=[WorldSizex*rand;WorldSizey*rand];
        P1(i)=atand(P(2,i)/P(1,i));
        w(1,i)=1/N;         
    end
    PCenter(1)=mean(P1);     
    % Start to track
    for k=2:len
        v1=normrnd(0,q1,1,1);             
        v2=normrnd(0,q2,1,1);               
        X(:,k)=X(:,k-1)+T*[vx;vy]+[v1;v2]; 
        XTH(k)=atand(X(2,k)/X(1,k));        
        XR(k)=sqrt(X(2,k)^2+X(1,k)^2); 
        % Create spurious measurements
        for i=1:nc
            W=normrnd(0,r_angle,1,1);
            RTH(k,i)=XTH(k)+W;              
        end
        for i=1:N
            v1=normrnd(0,q1,1,1);              
            v2=normrnd(0,q2,1,1);               
            P(:,i)=P(:,i)+T*[vx;vy]+[v1;v2];  
            P1(i)=atand(P(2,i)/P(1,i));       
            for j=1:nc
                v(i,j)=RTH(k,j)-P1(i);
                e(i,j)=exp(-0.5*(v(i,j))'*(r_angle)^(-1)*(v(i,j)));
            end
            b=(1-Pd)/Pd*lamda*sqrt((2*pi)^nc*r_angle);
            l(i)=sum(e(i,:))+b;      
            wh(i)=w(k-1,i)*l(i);     
        end
        
        % Calculate and normalize the weights
        for i=1:N
            w(k,i)=wh(i)/sum(wh);
            PCenter1(i)=w(k,i)*P1(i);
        end
        PCenter(k)=sum(PCenter1);
        
        % Resampling
        for i=1:N 
            wmax=rand;
            wss=0;
            for j=1:N
                wss=wss+w(k,j);
                if wss>=wmax
                    P(i)=P(j);
                    break;
                end
            end
        end
        
    end
    
    % Calculte error
    err(m,:)=abs(PCenter-XTH);
    if m==1
        RTH1=mean(RTH,2);   
        figure; 
        set(gca,'FontSize',12);
        plot(1:len, XTH, 'r*-',1:len, RTH1 , 'g.-',1:len,PCenter , 'bo-');
        legend('True State', 'Measurements central position', 'PF-PDA result');
        xlabel('t/s'); ylabel('ж╚/бу');
        grid on
        title('Particle Filter')       
    end
end
err1=mean(err);
%%
figure;  % Error plotting
set(gca,'FontSize',12);
plot(err1,'.-');
xlabel('t/s');
ylabel('RMSE/бу');
% legend('Particle Filter')
title('The err');
grid on