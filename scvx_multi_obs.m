clc 
clear all
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
        u(:,count)=5*[1;1];
    else
        u(:,count)=5*[1;1];
    end
    X(:,count+1)=Ad*X(:,count)+Bd*u(:,count);
    count=count+1;   
end

X0=X;
plot(X(1,:),X(2,:),'.')
hold on

%%
N=length(X(1,:));
obs_center=[7,9; ...
    14.1,11.1];
R=[2,2.5];

% obs_center=[7,9];
% R=[2];
% obs_center=[5,5;12.1,12.1];
% R=[4,2];
alpha=1;
obs_num=length(R);
r_default=1;

lambda=10000;
rho0=0.01;
rho1=0.2;
rho2=0.9;
i=1;
figure(1)
tol=0.001;



hold on

theta=linspace(0,2*pi,201);
for j=1:obs_num
    x_theta=R(j)*cos(theta);
    y_theta=R(j)*sin(theta);
    plot(obs_center(j,1)+x_theta,obs_center(j,2)+y_theta)
end

for k=1:60
    



    cvx_begin 
        
        variable w(2,N-1)

        variable v(2,N-1)
        variable d(2,N)
        variable U(1,N-1)
        variable s(N,obs_num)
        variable X_diff(N-1,1)
        %variable s_pos(N)
        %minimize (  500*sum(U*Ts) + lambda*sum(sum(abs(v)))  + lambda*(sum(max(s(:,1),0))   +   1*sum(max(s(:,2),0))) )
        
        minimize (  50000*sum(X_diff) + lambda*sum(sum(abs(v)))  + lambda*(   sum(sum(max(s,0)))   )) 
        %minimize (  50000*sum(X_diff) + lambda*sum(sum(abs(v)))) 
        %minimize (  50000*sum(X_diff) + lambda*sum(sum(abs(v)))  + lambda^2*(sum(max(s(:,1),0))   +   1*sum(max(s(:,2),0))) )
        subject to
        E=eye(2);
        d(:,1)==[0;0];

        for i=1:N-1
            X(:,i+1)+d(:,i+1)==(Ad*X(:,i)+Ad*d(:,i))+(Bd*u(:,i)+Bd*w(:,i))+E*v(:,i);
            U(i)==norm(u(:,i),2);

            X_diff(i)==norm(X(1:2,i+1)-X(1:2,i),2);

            -r_default<=w(1,i)<=r_default;
            -r_default<=w(2,i)<=r_default;
            
            dh=2*(X(1:2,i)-obs_center);
            h=norm(X(1:2,i)-obs_center,2);



            %-dh*u(:,i)-alpha*h<=0;


            % 
            for j=1:obs_num
                R(j)-norm(X(1:2,i)-obs_center(j,:)',2)- ...
                    (X(1:2,i)-obs_center(j,:)')'*(X(1:2,i)+d(1:2,i)-obs_center(j,:)') ...
                    /norm(X(1:2,i)-obs_center(j,:)',2)<=s(i,j);
            end
            

        end

        X(:,N)+d(:,N)==[20;20];
        
    cvx_end

    % 


    w=full(w);
    v=full(v);
    d=full(d);
    X=X+d;
    u=u+w;
    hold on
    plot(X(1,:),X(2,:),'.')
    ss=0;

    for j=1:obs_num
        for i=1:N-1
            ss(i,j)=R(j)-norm(X(1:2,i)-obs_center(j,:)',2);
        end
    end

    
    %ss=reshape(ss,[1,(N-1)*obs_num]);
    if max(max(ss))<0 && k>4
        break;
    end
    pause(0.01)
end


%%
figure(2)
hold on
plot(X(1,:),X(2,:),'.')
plot(X0(1,:),X0(2,:),'.')
theta=linspace(0,2*pi,201);
xlim([0 20])
ylim([0 20])
for j=1:obs_num
    x_theta=R(j)*cos(theta);
    y_theta=R(j)*sin(theta);
    plot(obs_center(j,1)+x_theta,obs_center(j,2)+y_theta)
end






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