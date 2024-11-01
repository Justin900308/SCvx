clc 
clear all
close all
clf

A=[0 0;0 0];
B=[1 0;0 1];
C=eye(2,2);
D=zeros(2,2);

sys=ss(A,B,C,D);

Ts=0.1;
sysd=c2d(sys,Ts);
Ad=sysd.A;
Bd=sysd.B;
T=50;
%% intial trajectory
count=1;
X(:,1)=[0;0];
for t=0:Ts:5
    if t<2
        u(:,count)=1*0*[1;1];
    else
        u(:,count)=1*0*[1;1];
    end
    X(:,count+1)=Ad*X(:,count)+Bd*u(:,count);
    count=count+1;   
end


plot(X(1,:),X(2,:),'.')
hold on

%%


P_des=[10;10];


N=length(X(1,:));
obs_center=[5;4];
R=3.5;
alpha=1;

r_default=0.3;

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
for k=1:100
    

    cvx_solver SDPT3
    cvx_precision best
    %cvx_solver sedumi
    cvx_begin 
        
        variable w(2,N-1)

        variable v(2,N-1)
        variable d(2,N)
        %variable U(2,N-1)
        variable s(N-1)
        %s=zeros(N,1);
        %variable s_pos(N)
        minimize (  0.1*sum(sum(abs((u+w)*Ts))) + lambda*sum(sum(abs(v)))  + 1*lambda*sum(max(s,0)) )
        %minimize (  01*sum(sum(abs((U)*Ts))) + lambda*sum(sum(abs(v)))  + 1*lambda*sum(max(s)) )
  
        subject to
        E=eye(2);
        
        X(:,1)+d(:,1)==[0;0];
        for i=1:N-1
            %s(i)<=-R*0.5*0;
            
            
            
            X(:,i+1)+d(:,i+1)==(Ad*X(:,i)+Ad*d(:,i))+(Bd*u(:,i)+Bd*w(:,i))+E*v(:,i);
            %U(:,i)  ==     0.5*(P_des-X(:,i))     -    (u(:,i)+w(:,i));  % desired velocity - actual velocity


            %v_des=([10;10]-(X(1:2,i)+d(:,i)));
            %U(i)==norm([00;0]-(X(1:2,i)))  -  (u(:,i)+w(:,i));

            
            -r_default<=w(1,i)<=r_default;
            -r_default<=w(2,i)<=r_default;
            dh=2*(X(1:2,i)-obs_center);
            h=norm(X(1:2,i)-obs_center,2);
            %-0.1<=d(3:4,i)<=0.1;


            %R-norm(X(1:2,i)-obs_center,2)-(X(1:2,i)-obs_center)'*(X(1:2,i)+d(1:2,i)-obs_center)/norm(X(1:2,i)-obs_center,2)<=s(i);
            R-norm(X(1:2,i)-obs_center,2)-(X(1:2,i)-obs_center)'*(X(1:2,i)+d(1:2,i)-obs_center)/norm(X(1:2,i)-obs_center,2)<=s(i);

            %R-(X(1:2,i)+d(1:2,i)-obs_center)/norm(X(1:2,i)-obs_center,2)*(X(1:2,i)-obs_center)'<=s(i);
            %R-norm(X(1:2,i)-obs_center,2)-(X(1:2,i)-obs_center)'*(X(1:2,i)-obs_center)/norm(X(1:2,i)-obs_center,2)
        end
        
        X(:,N)+d(:,N)==P_des;
        
    cvx_end

    % 


    w=full(w);
    v=full(v);
    d=full(d);
    X=X+d;
    u=u+w;
    hold on
    plot(X(1,:),X(2,:),'.')
    for i=1:N-1

        ss(i)=R-norm(X(1:2,i)-obs_center,2);

    end
    pause(0.01)
    if max(ss)<0 && k>10
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





