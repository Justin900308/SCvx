clc 
clear all
clf

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
        u(:,count)=1*[1;01];
    else
        u(:,count)=1*[1;1];
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
N=length(X(1,:));
obs_center=[5;5];
R=3;
alpha=1;

r_default=0.8;

lambda=10000;
rho0=0.01;
rho1=0.2;
rho2=0.9;
i=1;
figure(1)
tol=0.001;





theta=linspace(0,2*pi,201);
x_theta=R*cos(theta);
y_theta=R*sin(theta);
plot(obs_center(1)+x_theta,obs_center(2)+y_theta)
for k=1:60
    



    cvx_begin quiet
        
        variable w(2,N-1)

        variable v(4,N-1)
        variable d(4,N)
        variable U(1,N-1)
        %variable s(N)
        %variable s_pos(N)
        minimize (  500*sum(U*Ts) + lambda*sum(sum(abs(v)))   )
        subject to
        E=eye(4);
        d(:,1)==[0;0;0;0];
        
        for i=1:N-1
            X(:,i+1)+d(:,i+1)==(Ad*X(:,i)+Ad*d(:,i))+(Bd*u(:,i)+Bd*w(:,i))+E*v(:,i);
            U(i)==norm(u(:,i),2);
            
            -r_default<=w(1,i)<=r_default;
            -r_default<=w(2,i)<=r_default;
            dh=2*(X(1:2,i)-obs_center);
            h=norm(X(1:2,i)-obs_center,2);
            %-0.1<=d(3:4,i)<=0.1;


            R-norm(X(1:2,i)-obs_center,2)-(X(1:2,i)-obs_center)'*(X(1:2,i)+d(1:2,i)-obs_center)/norm(X(1:2,i)-obs_center,2)<=0;
            %R-norm(X(1:2,i)-obs_center,2)-(X(1:2,i)-obs_center)'*(X(1:2,i)-obs_center)/norm(X(1:2,i)-obs_center,2)==s(i);
        end
        
        X(:,N)+d(:,N)==[10;10;0;0];
        
    cvx_end

    % 
    for i=1:N-1

        s(i)=R-norm(X(1:2,i)-obs_center,2);

    end

    w=full(w);
    v=full(v);
    d=full(d);
    X=X+d;
    u=u+w;
    hold on
    plot(X(1,:),X(2,:),'.')
    pause(0.01)
    if max(s)<0 && k>4
        break;
    end
end



figure(2)
hold on
plot(X(1,:),X(2,:),'.')
theta=linspace(0,2*pi,201);
x_theta=R*cos(theta);
y_theta=R*sin(theta);
plot(obs_center(1)+x_theta,obs_center(2)+y_theta)


%%

        for i=1:N-1

            R-norm(X(1:2,i)-obs_center,2)-(X(1:2,i)-obs_center)'*(X(1:2,i)+d(1:2,i)-obs_center)/norm(X(1:2,i)-obs_center,2)

        end
        i=29
%%

%X=[5.1;5;0;0]
% R-norm(X(1:2)-obs_center,2)-(X(1:2)-obs_center)'*(X(1:2)-obs_center)/norm(X(1:2)-obs_center,2)
% 
% lambda*sum(s(s>=0))
%(X(1:2,i)-obs_center)'*(X(1:2,i)-obs_center)*d(:,i)

% %%
% obs_center=[5;5];
% R=1;
% alpha=1;
% 
% 
% 
% cvx_begin
% 
%     variable uu(102,1)
%     variable XX(204,1)
% 
% 
%     minimize(norm(uu,2))
% 
%     XX(1:4)==[0;0;0;0];
%     uu(1:2)==[0;0];
%     for i=2:51
%         XX((i-1)*4+1:(i-1)*4+4)==Ad*XX((i-1-1)*4+1:(i-1-1)*4+4)+Bd*uu((i-1-1)*2+1:(i-1-1)*2+2);
%         h=norm(XX((i-1-1)*4+1:(i-1-1)*4+2)-obs_center,2)-R;
%         dh=2*(XX((i-1-1)*4+1:(i-1-1)*4+2)-obs_center)';
% 
%     end
% 
% 
%     XX(201:204)==[10;10;0;0];
% 
% cvx_end
% 
% [X,u]=unstack_vec(XX,uu);
% 

% 
% 
% 
% 
% 
% 
% 
% %% 
% syms x y xobs yobs
% P=[x;y]
% P_obs=[xobs;yobs]
% D=norm(P-P_obs)
% pretty(diff(D,x))
% pretty(diff(D,y))



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