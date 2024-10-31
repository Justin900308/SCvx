clc 
clear all
clf

A=[0 0;0 0];
B=[1 0;0 1];
C=eye(2,2);
D=zeros(2,2);

sys=ss(A,B,C,D);


Num_agen=4;


Ts=0.1;
sysd=c2d(sys,Ts);
Ad=sysd.A;
Bd=sysd.B;
T=50;
%% intial trajectory
count=1;
X(:,1)=[0;0 ;  20;0 ; 20;20   ;  0;20 ];
for t=0:Ts:5
    if t<2
        u(:,count)=0*ones(2*Num_agen,1);
    else
        u(:,count)=0*ones(2*Num_agen,1);
    end
    for i=1:Num_agen
        X((2*i-1):(2*i),count+1)=Ad*X((2*i-1):(2*i),count)+Bd*u((2*i-1):(2*i),count);
    end
    count=count+1;   
end

X0=X;
for i=1:Num_agen
    hold on
    plot(X((2*i-1),:),X((2*i),:),'.')
end


hold on

%%
R_obs=2;
R_agent=0.5;
N=length(X(1,:));
obs_center=[
    X(1:2,1)'; ...
    X(3:4,1)'; ...
    X(5:6,1)'; ...
    X(7:8,1)'];
R_plot=[R_agent,R_agent,R_agent,R_agent];
R=2*R_plot;


alpha=1;
obs_num=length(R);
r_default=ones(Num_agen,1);

lambda=10000;

i=1;
figure(1)
tol=0.001;



hold on

theta=linspace(0,2*pi,201);
for j=1:obs_num
    x_theta=R_plot(j)*cos(theta);
    y_theta=R_plot(j)*sin(theta);
    plot(obs_center(j,1)+x_theta,obs_center(j,2)+y_theta)
end

for iteration=1:60
    




    cvx_begin 
        
        variable w(2*Num_agen,N-1)

        variable v(2*Num_agen,N-1)
        variable d(2*Num_agen,N)
        variable U(1,N-1)
        variable s(N*(obs_num-1),Num_agen)

        variable X_diff(N-1,Num_agen)
        %variable s_pos(N)
        %minimize (  500*sum(U*Ts) + lambda*sum(sum(abs(v)))  + lambda*(sum(max(s(:,1),0))   +   1*sum(max(s(:,2),0))) )
        
        minimize (  500*sum(sum(X_diff)) + lambda*sum(sum(abs(v)))  + lambda*0.001*(   sum(sum(max(s,0)))   )) 

        subject to
        cvx_precision best
        %cvx_precision low
        E=eye(2);
        d(:,1)==zeros(2*Num_agen,1);




        
        for i=1:N-1

            obs_center=[
                X(1:2,i)'; ...
                X(3:4,i)'; ...
                X(5:6,i)'; ...
                X(7:8,i)'];

            for j=1:Num_agen


                -r_default(j,1)<=w((2*j-1):(2*j),i)<=r_default(j,1);
  


                X((2*j-1):(2*j),i+1)+d((2*j-1):(2*j),i+1)== ...
                    (Ad*X((2*j-1):(2*j),i)+Ad*d((2*j-1):(2*j),i))+ ...
                    (Bd*u((2*j-1):(2*j),i)+Bd*w((2*j-1):(2*j),i))+E*v((2*j-1):(2*j),i);
                U(i)==norm(u(:,i),2);

                X_diff(i,j)==norm(X((2*j-1):(2*j),i+1)-X((2*j-1):(2*j),i),2);




                countk=1;

                for k=1:obs_num
                    if k==j
                        continue
                    end
                    R(k)-norm(X((2*j-1):(2*j),i)-obs_center(k,:)',2)- ...
                        (X((2*j-1):(2*j),i)-obs_center(k,:)')'*(X((2*j-1):(2*j),i)+d((2*j-1):(2*j),i)-obs_center(k,:)') ...
                        /norm(X((2*j-1):(2*j),i)-obs_center(k,:)',2)<=s((countk-1)*(N)+i,j);
                    countk=countk+1;
                end
            end

        end

        X(:,N)+d(:,N)==[20;20  ;  0;20   ;   0;0   ;  20;0];
        
    cvx_end

    % 


    w=full(w);
    v=full(v);
    d=full(d);
    X=X+d;
    u=u+w;
    hold on
    for i=1:Num_agen
        hold on
        plot(X((2*i-1),:),X((2*i),:),'.')
    end
    ss=0;



    for i=1:N
        
        obs_center=[
            X(1:2,i)'; ...
            X(3:4,i)'; ...
            X(5:6,i)'; ...
            X(7:8,i)'];
        for j=1:Num_agen
            countk=1;
            for k=1:obs_num
                if k==j
                    continue
                end
                ss((countk-1)*(N)+i,j)=R(k)-norm(X((2*j-1):(2*j),i)-obs_center(k,:)',2);
                countk=countk+1;
            end
        end
        
    end

    ss_max=max(ss);

    for i=1:Num_agen
        if ss_max(i)<0
            r_default(i)=0;
        else
            r_default(i)=1;
        end
    end

    
    
    if max(max(ss))<0 && iteration>4
        break;
    end
    pause(0.01)
end


%%
figure(2)
hold on
for i=1:N
    for j=1:Num_agen
        hold on
        plot(X((2*j-1),i),X((2*j),i),'g.')
        %plot(X0((2*i-1),:),X0((2*i),:),'.')
    end

    obs_center=[
        X(1:2,i)'; ...
        X(3:4,i)'; ...
        X(5:6,i)'; ...
        X(7:8,i)'];


    theta=linspace(0,2*pi,201);
    xlim([0 25])
    ylim([0 25])
    for j=1:obs_num
        hold on
        x_theta=R_plot(j)*cos(theta);
        y_theta=R_plot(j)*sin(theta);
        plot(obs_center(j,1)+x_theta,obs_center(j,2)+y_theta,'b')
    end

pause(0.05)
%clf


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