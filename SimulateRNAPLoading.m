function SimulateRNAPLoading

% AUTHOR : Yang Joon Kim (yjkim90@berkeley.edu)
% DESCRIPTION :
% This script is for simulating MS2 traces based on specific RNAP loading
% rate profile (time trace). Now, I will assume that it's nc 13 (~ 15 min), 
% and the elongation rate as Garcia, 2013, (~1.5kb/min),
% and the gene length to be ~4.5kb, which gives us elongation time of 3
% min. So, the number of RNAPs loaded at time T should fall off at (T + 3)
% min. We assumed that the Txn starts at T=0, and ends in T=12 min.

% Define the rates
r1 = 10; % #/min
r2 = 3; % #/min
% Define the Time space
Time = linspace(0,15,100);
Tstop = 12; % min
Telongate = 3; % min
dt = Time(2)-Time(1);
deltaT = ceil(Telongate/dt);
N = zeros((size(Time)));
N2 = zeros((size(Time)));

% 1) linear rate
Rate = 10*(ones(size(Time)))-0.5*Time;

% 2) exponetial rate
% Amp = 10;
% Decay = 10;

Amp2 = 10;
Decay2 = 8;

% 3) Gaussian rate (or Dirac-Delta function-like?)

% Rate = Amp*exp(-(Time)/Decay);
Rate2 = Amp2*exp(-(Time)/Decay2);

% 4) parabola rate
%Rate = (15*(ones(size(Time)))-Time).^2 /22;

% Calculate the number of RNAP on the gene
% In here, we consider the falling off rate to be same as the loading rate
% in the past (T-T_elongation)
for i=2:length(Time)
    if Time(i)<3
        N(i) = N(i-1) + Rate(i)*dt;
    elseif Time(i)>=3 && Time(i) < Tstop
        N(i) = N(i-1) + Rate(i)*dt - Rate(i-deltaT)*dt;
    else
        N(i) = N(i-1) - Rate(i-deltaT)*dt;
    end
end

for i=2:length(Time)
    if Time(i)<3
        N2(i) = N2(i-1) + Rate2(i)*dt;
    elseif Time(i)>=3 && Time(i) < Tstop
        N2(i) = N2(i-1) + Rate2(i)*dt - Rate2(i-deltaT)*dt;
    else
        N2(i) = N2(i-1) - Rate2(i-deltaT)*dt;
    end
end

figure(1) 
hold on
plot(Time,Rate)
plot(Time,Rate2)

title('RNAP loading rate')
xlabel('Time (min)')
ylabel('RNAP loading Rate (AU)')
legend('linear','exp')
%% Plot the Number of RNAP loaded on the gene
figure(2)
hold on
plot(Time,N)
plot(Time,N2)

title('Loaded RNAP on the gene (MS2 fluo)')
xlabel('Time (min)')
ylabel('Number of RNAP (AU)')
legend('linear','exp')

%% Test different scenarios to break the model
% First, I will think about how adding more Repressor binding sites would
% change the r(t) form
N0 = zeros((size(Time)));
N = zeros((size(Time)));
N2 = zeros((size(Time)));

% 1) linear rate
%Rate = 10*(ones(size(Time)))-0.5*Time;

% 2) exponetial rate
Amp = 10;
Decay = 5;

Amp2 = 10;
Decay2 = 8;

Rate = Amp*exp(-(Time)/Decay);
Rate2 = Amp2*exp(-(Time)/Decay2);

% 3) Gaussian rate (or Dirac-Delta function-like?)

% 4) parabola rate
%Rate = (15*(ones(size(Time)))-Time).^2 /22;

% Calculate the number of RNAP on the gene
% In here, we consider the falling off rate to be same as the loading rate
% in the past (T-T_elongation)

Rate0 = Amp*(ones(size(Time)));
for i=2:length(Time)
    if Time(i)<3
        N0(i) = N0(i-1) + Rate0(i)*dt;
    elseif Time(i)>=3 && Time(i) < Tstop
        N0(i) = N0(i-1) + Rate0(i)*dt - Rate0(i-deltaT)*dt;
    else
        N0(i) = N0(i-1) - Rate0(i-deltaT)*dt;
    end
end

for i=2:length(Time)
    if Time(i)<3
        N(i) = N(i-1) + Rate(i)*dt;
    elseif Time(i)>=3 && Time(i) < Tstop
        N(i) = N(i-1) + Rate(i)*dt - Rate(i-deltaT)*dt;
    else
        N(i) = N(i-1) - Rate(i-deltaT)*dt;
    end
end

for i=2:length(Time)
    if Time(i)<3
        N2(i) = N2(i-1) + Rate2(i)*dt;
    elseif Time(i)>=3 && Time(i) < Tstop
        N2(i) = N2(i-1) + Rate2(i)*dt - Rate2(i-deltaT)*dt;
    else
        N2(i) = N2(i-1) - Rate2(i-deltaT)*dt;
    end
end

figure(1) 
hold on
plot(Time,Rate)
plot(Time,Rate2)
plot(Time,Rate0)

title('RNAP loading rate')
xlabel('Time (min)')
ylabel('RNAP loading Rate (AU)')
legend('exp1','exp2','0')
%% Plot the Number of RNAP loaded on the gene
figure(2)
hold on
plot(Time,N)
plot(Time,N2)
plot(Time,N0)

title('Loaded RNAP on the gene (MS2 fluo)')
xlabel('Time (min)')
ylabel('Number of RNAP (AU)')
legend('exp1','exp2','constant')
end