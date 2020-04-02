%Modeling for Repression
clear all
%First, I want to calculate the Mean and Noise of (level of)transcription,
%for different concentration of repressors.

%Define variables and parameters

%R=[0.001,0.01,0.1,1,10,100]; %Rep. conc. divided by KR (dissociation constant for repressor)
R=logspace(-3,2);
%A=[0.001,0.01,0.1,1,10,100]; %Act. conc. divided by KA (dissociation constant for activator)
A=logspace(-3,2);
w=1; %Interaction between activator and repressor

r= 1; %rate of transcription initiation with activator (molecules/min?)
rR= 0.2; %rate of transcription initiation with activator and repressor

for i=1:length(R)
    for j=1:length(A)
        Mean(i,j)=(r*A(j)+rR*A(j)*R(i)*w)/(1+A(j)+R(i)+A(j)*R(i)*w);
    end
end

for i=1:length(R)
    for j=1:length(A)
        Fold(i,j)=((r*A(j)+rR*A(j)*R(i)*w)/(1+A(j)+R(i)+A(j)*R(i)*w))/((r*A(j))/(1+A(j)));
    end
end

%% plot for repressor conc.
kStart=1;
kEnd=length(R);
kTotal=kEnd-kStart;

 colormap(jet(256));
 cmap=colormap;
 Color=cmap(round(((kStart:kEnd)-kStart)/(kEnd-kStart)*255)+1,:);


for k=1:length(A)
    semilogx(R,Mean(:,k),'color',Color(k-kStart+1,:))
    hold on
    title('Gene Expression for repressor concentration')
    xlabel('Repressor Concentration (R/K_R)')
    ylabel('Rate of Gene Expression')
    set(gca,'Fontsize',30)
    %ylim([0.2 1])
    pause
    (A(k))
    
end

%hold off
colorbar
%% Rate for w (interaction between Act. and Rep.)
A(40);
r=100;
rR=50;

w=[0:0.5:1];

kStart=1;
kEnd=length(w);
kTotal=kEnd-kStart;

colormap(jet(256));
cmap=colormap;
Color=cmap(round(((kStart:kEnd)-kStart)/(kEnd-kStart)*255)+1,:);
for k=1:length(w)
    for i=1:length(R)
        Mean(i,k)=(r*A(40)+rR*A(40)*R(i)*w)/(1+A(40)+R(i)+A(40)*R(i)*w);
    
    end
end

hold on

for k=1:length(w)
    semilogx(R,Mean(:,k),'color',Color(k-kStart+1,:))
    %hold on
    title('Fold-change for Repressor')
    xlabel('Repressor Concentration (R/K_R)')
    ylabel('Fold-change')
    %ylim([0 1])
    pause
    set(gca,'Fontsize',30)
end
legend('w=0','w=0.5','w=1')

%% 
for i=1:length(R)
    Mean(i,1)=(r*A(40))/(1+A(40)+R(i));
end

semilogx(R,Mean(:,1))
hold on
for i=1:length(R)
    Mean(i,2)=(r*A(40)+rR*A(40)*R(i)*1)/(1+A(40)+R(i)+A(40)*R(i)*1);
end

semilogx(R,Mean(:,2))

title('Rate of Transcription Initiation')
xlabel('Repressor Concentration (R/K_R)')
ylabel('Rate(AU)')
legend('w=0','w=1')
set(gca,'Fontsize',30)

%% Prediction for hbP2- using Hill function

na=6;

rate=1;
rR=0.1;

A=[0.001,0.01,0.1,1,3,10,100];
R=[0.001,0.01,0.1,1,10,100];

AP=0:0.01:1;
a=5*exp(-3*AP);

%Repressor conc.
r=-R(5)*1*(AP-0.5).^2+4;
% for i=1:length(AP)
%     if AP(i)<0.5
%         r(i)=R(4)*AP(i);
%     else
%         r(i)=R(4)-R(4)*AP(i);
%     end
% end
w=0 ;
for nr=1:3
    for j=1:length(AP)
        Rate(nr,j)=(rate*(a(j))^(na)+(rR)*w*(a(j))^(na)*(r(j))^nr)/(1+(a(j))^(na)+(r(j))^nr+w*(a(j))^(na)*(r(j))^nr);
    end
end

for j=1:length(AP)
        Rate(4,j)=(rate*(a(j))^na)/(1+(a(j))^na);
end


figure(1)
plot(AP,a,'b')
hold on
plot(AP,r,'r')
title('Input transcription factors')
xlabel('AP')
ylabel('Concentration/K_D (AU)')
legend('Activator','Repressor')
set(gca,'fontsize',30)
hold off

figure(2)
hold on
plot(AP,Rate(4,:),'k') %nr=0

plot(AP,Rate(1,:),'b')
plot(AP,Rate(2,:),'g')
plot(AP,Rate(3,:),'r')


title('Prediction for the rate of transcription')
xlabel('AP')
ylabel('Rate (AU)')
legend('0','1','2','3')
set(gca,'fontsize',30)
