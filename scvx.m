clc 
clear all

A=[0 0 1 0; 0 0 0 1; 0 0 0 0; 0 0 0 0];
B=[0 0; 0 0; 1 0;0 1];
C=eye(4,4);
D=zeros(4,2);

sys=ss(A,B,C,D);

Ts=0.1;
sysd=c2d(sys,Ts);
Ad=sysd.A;
Bd=sysd.B;
T=50;
%% intial trajectory
count=1;
X(:,1)=[0;0;0;0];
for t=0:Ts:5
    if t<2
        u(:,count)=[0.1;0.5];
    else
        u(:,count)=[0.5;0.1];
    end
    X(:,count+1)=Ad*X(:,count)+Bd*u(:,count);
    count=count+1;
end

XX=reshape(X(:,1:end-1),[length((X(:,1)))*length(X(1,:))-4,1]);
uu=reshape(u,[length((u(:,1)))*length(u(1,:)),1]);
z=stack_vec(X,u);
plot(X(1,:),X(2,:),'.')
hold on

%%


cvx_begin

    variable uu(102,1)
    variable XX(204,1)
    minimize(norm(uu,2))

    XX(1:4)==[0;0;0;0];
    uu(1:2)==[0;0];
    for i=2:51
        XX((i-1)*4+1:(i-1)*4+4)==Ad*XX((i-1-1)*4+1:(i-1-1)*4+4)+Bd*uu((i-1-1)*2+1:(i-1-1)*2+2);
        -1<=uu((i-1-1)*2+1)<=1;
        -1<=uu((i-1-1)*2+1)<=1;
    end


    XX(201:204)==[3.2225;4.5805;0;0];

cvx_end

[X,u]=unstack_vec(XX,uu);

plot(X(1,:),X(2,:),'.')







%% 
syms x y xobs yobs
P=[x;y]
P_obs=[xobs;yobs]
D=norm(P-P_obs)
pretty(diff(D,x))
pretty(diff(D,y))



%% utility funcitons


function XX=stack_vec(X,u)
    XX=reshape(X(:,1:end-1),[length((X(:,1)))*length(X(1,:))-4,1]);
    uu=reshape(u,[length((u(:,1)))*length(u(1,:)),1]);
    XX=[XX;uu];
end


function [X,u]=unstack_vec(XX,uu)
    X=reshape(XX,[4,length(XX)/4]);
    u=reshape(uu,[2,length(uu)/2]);
end