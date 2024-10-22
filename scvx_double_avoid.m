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


count=1;
X2(:,1)=[10;10;0;0];
for t=0:Ts:5
    if t<2
        u2(:,count)=-1*[1;01];
    else
        u2(:,count)=-1*[1;1];
    end
    X2(:,count+1)=Ad*X2(:,count)+Bd*1*u(:,count);
    count=count+1;   
end


plot(X2(1,:),X2(2,:),'.')
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
        variable s(1,N-1)

        variable w2(2,N-1)
        variable v2(4,N-1)
        variable d2(4,N)
        variable U2(1,N-1)
        variable s2(1,N-1)
        %variable s(N)
        %variable s_pos(N)
        minimize (  500*sum(U*Ts) + lambda*sum(sum(abs(v))) + 500*sum(U2*Ts) + lambda*sum(sum(abs(v2))) ...
            + lambda*sum(max(s,0)) + lambda*sum(max(s2,0)))
        subject to
        E=eye(4);
        d(:,1)==[0;0;0;0];
        d2(:,1)==[0;0;0;0];
        for i=1:N-1
            X(:,i+1)+d(:,i+1)==(Ad*X(:,i)+Ad*d(:,i))+(Bd*u(:,i)+Bd*w(:,i))+E*v(:,i);
            U(i)==norm(u(:,i),2);
            
            -r_default<=w(1,i)<=r_default;
            -r_default<=w(2,i)<=r_default;
            c1=w(1,i);
            c2=w(2,i);
            
            %0.001<=c1<=r_default
            %0.001<=c1<=r_default
            R-norm(X(1:2,i)-X2(1:2,i),2)-(X(1:2,i)-X2(1:2,i))'*(X(1:2,i)+d(1:2,i)-X2(1:2,i))/norm(X(1:2,i)-X2(1:2,i),2)<=s(i);

            %
            X2(:,i+1)+d2(:,i+1)==(Ad*X2(:,i)+Ad*d2(:,i))+(Bd*u2(:,i)+Bd*w2(:,i))+E*v2(:,i);
            U2(i)==norm(u2(:,i),2);
            
            -r_default<=w2(1,i)<=r_default;
            -r_default<=w2(2,i)<=r_default;
            R-norm(X2(1:2,i)-X(1:2,i),2)-(X2(1:2,i)-X(1:2,i))'*(X2(1:2,i)+d2(1:2,i)-X(1:2,i))/norm(X2(1:2,i)-X(1:2,i),2)<=s2(i);
        end
        
        X(:,N)+d(:,N)==[10;10;0;0];
        X2(:,N)+d2(:,N)==[0;0;0;0];
    cvx_end

    % 

    w=full(w);
    v=full(v);
    d=full(d);
    X=X+d;
    u=u+w;

    w2=full(w2);
    v2=full(v2);
    d2=full(d2);
    X2=X2+d2;
    u2=u2+w2;
    hold on
    for i=1:N-1

        ss(i)=R-norm(X(1:2,i)-X2(1:2,i),2);
        ss2(i)=R-norm(X2(1:2,i)-X(1:2,i),2);
    end
    
    plot(X(1,:),X(2,:),'.')
    plot(X2(1,:),X2(2,:),'.')
    pause(0.01)
    if max(ss)<0 && k>4 && max(ss2)<0 
        break;
    end
end
R=R/2;
%%
clf
for i=1:N

hold on
plot(X(1,1:i),X(2,1:i),'r.')
plot(X2(1,1:i),X2(2,1:i),'b.')
theta=linspace(0,2*pi,201);
x_theta=R*cos(theta);
y_theta=R*sin(theta);
plot(X(1,i)+x_theta,X(2,i)+y_theta,'r')
plot(X2(1,i)+x_theta,X2(2,i)+y_theta,'b')

xlim([0 15])
ylim([0 15])
pause(0.01)


end

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