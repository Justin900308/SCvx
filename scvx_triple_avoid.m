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
T=100;
%% intial trajectory
count=1;
X(:,1)=[0;0;0;0];
for t=0:Ts:5
    if t<2
        u(:,count)=0.*[1;01];
    else
        u(:,count)=0.*[1;1];
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
        u2(:,count)=0.*[1;01];
    else
        u2(:,count)=0.*[1;1];
    end
    X2(:,count+1)=Ad*X2(:,count)+Bd*1*u(:,count);
    count=count+1;   
end


plot(X2(1,:),X2(2,:),'.')
hold on

count=1;
X3(:,1)=[-5;5;0;0];
for t=0:Ts:5
    if t<2
        u3(:,count)=0.*[1;01];
    else
        u3(:,count)=0.*[1;1];
    end
    X3(:,count+1)=Ad*X3(:,count)+Bd*1*u3(:,count);
    count=count+1;   
end


plot(X3(1,:),X3(2,:),'.')
hold on
%%
N=length(X(1,:));
obs_center=[5.1;5];
R=4;
alpha=1;

r_default=1.8;

lambda=10000;
rho0=0.01;
rho1=0.2;
rho2=0.9;
i=1;
figure(1)
tol=0.001;





theta=linspace(0,2*pi,201);
x_theta=1.5*R*cos(theta);
y_theta=1.5*R*sin(theta);
plot(obs_center(1)+x_theta,obs_center(2)+y_theta)
final_center=[20;20];
for k=1:60
    



    cvx_begin 
        
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

        variable w3(2,N-1)
        variable v3(4,N-1)
        variable d3(4,N)
        variable U3(1,N-1)
        variable s3(1,N-1)

        variable P(3,N-1)
        %variable s(N)
        %variable s_pos(N)
        minimize (  500*sum(U*Ts) + lambda*sum(sum(abs(v))) + 500*sum(U2*Ts) + lambda*sum(sum(abs(v2))) ...
            + lambda*sum(max(s,0)) + lambda*sum(max(s2,0)) + ...
            + 500*sum(U3*Ts) + lambda*sum(sum(abs(v3))) ...
            + lambda*sum(max(s3,0)) )
        subject to
        E=eye(4);
        d(:,1)==[0;0;0;0];
        d2(:,1)==[0;0;0;0];
        d3(:,1)==[0;0;0;0];
        for i=1:N-1
            X(:,i+1)+d(:,i+1)==(Ad*X(:,i)+Ad*d(:,i))+(Bd*u(:,i)+Bd*w(:,i))+E*v(:,i);
            U(i)==norm(u(:,i),2);
            
            -r_default<=w(1,i)<=r_default;
            -r_default<=w(2,i)<=r_default;
            c1=w(1,i);
            c2=w(2,i);
            
            %0.001<=c1<=r_default
            %0.001<=c1<=r_default
            1.5*R-norm(X(1:2,i)-obs_center,2)-(X(1:2,i)-obs_center)'*(X(1:2,i)+d(1:2,i)-obs_center)/norm(X(1:2,i)-obs_center,2)<=s(i);

            %
            X2(:,i+1)+d2(:,i+1)==(Ad*X2(:,i)+Ad*d2(:,i))+(Bd*u2(:,i)+Bd*w2(:,i))+E*v2(:,i);
            U2(i)==norm(u2(:,i),2);
            
            -r_default<=w2(1,i)<=r_default;
            -r_default<=w2(2,i)<=r_default;
            1.5*R-norm(X2(1:2,i)-obs_center,2)-(X2(1:2,i)-obs_center)'*(X2(1:2,i)+d2(1:2,i)-obs_center)/norm(X2(1:2,i)-obs_center,2)<=s2(i);

            %
            X3(:,i+1)+d3(:,i+1)==(Ad*X3(:,i)+Ad*d3(:,i))+(Bd*u3(:,i)+Bd*w3(:,i))+E*v3(:,i);
            U3(i)==norm(u3(:,i),2);
            
            -r_default<=w3(1,i)<=r_default;
            -r_default<=w3(2,i)<=r_default;
            1.5*R-norm(X3(1:2,i)-obs_center,2)-(X3(1:2,i)-obs_center)'*(X3(1:2,i)+d3(1:2,i)-obs_center)/norm(X3(1:2,i)-obs_center,2)<=s3(i);
        
            
            P(1,i)==5-norm(X(1:2,i)-X2(1:2,i),2);
            P(2,i)==5-norm(X2(1:2,i)-X3(1:2,i),2);
            P(3,i)==5-norm(X3(1:2,i)-X(1:2,i),2);
            % 
            % if i>110
            % X2(:,i)+d2(:,i)<=[X(1,i)+d(1,i)+5;X(2,i)+d(2,i);0;0];
            % X3(:,i)+d3(:,i)<=[X(1,i)+d(1,i)+2.5;X(2,i)+d(2,i)+5/2*sqrt(3);0;0];
            % end
        end

        
        X(:,N)+d(:,N)==[final_center(1);final_center(2);0;0];
        X2(:,N)+d2(:,N)==[X(1,N)+d(1,N)+5;X(2,N)+d(2,N);0;0];
        X3(:,N)+d3(:,N)==[X(1,N)+d(1,N)+2.5;X(2,N)+d(2,N)+5/2*sqrt(3);0;0];
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

    w3=full(w3);
    v3=full(v3);
    d3=full(d3);
    X3=X3+d3;
    u3=u3+w3;

    hold on
    for i=1:N-1

        ss(i)=1.5*R-norm(X(1:2,i)-obs_center,2);
        ss2(i)=1.5*R-norm(X2(1:2,i)-obs_center,2);
    end
    
    plot(X(1,:),X(2,:),'.')
    plot(X2(1,:),X2(2,:),'.')
    plot(X3(1,:),X3(2,:),'.')
    pause(0.01)
    if max(ss)<0 && k>3 && max(ss2)<0 
        break;
    end
end

%%
clf
R=4;
x_theta=R*cos(theta);
y_theta=R*sin(theta);
plot(5+x_theta,5+y_theta,'g')
R=R/2;
for i=1:N

hold on
plot(X(1,1:i),X(2,1:i),'r.')
plot(X2(1,1:i),X2(2,1:i),'b.')
plot(X3(1,1:i),X3(2,1:i),'g.')
theta=linspace(0,2*pi,201);
x_theta=R*cos(theta);
y_theta=R*sin(theta);
plot(X(1,i)+x_theta,X(2,i)+y_theta,'r')
plot(X2(1,i)+x_theta,X2(2,i)+y_theta,'b')
plot(X3(1,i)+x_theta,X3(2,i)+y_theta,'g')

xlim([-10 45])
ylim([-10 45])
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