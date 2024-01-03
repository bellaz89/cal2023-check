clear
clc
disp('Start with FORW / REFL simulation')
disp(' * * * * * * * * * * * * * * * * * * * * * ')
s = tf('s');         % laplace operator for linear simulation

% Define simulation parameters 
VF_scalar(1) = 12.14;
VF_scalar(2) = 6;
VF_scalar(3) = 0;

% Klystron noise amplitude
kly_noise_std = 0.01

% Measurement noise amplitude
meas_noise_std = 0.01

t_fill = 750e-6;
t_flat = 650e-6;
t_decay = 400e-6;
fs = 10e6;

w12 = 2*pi*141.3;
dw = 0;

% build forward signal
VF(1:fs*t_fill,1) = VF_scalar(1);
VF(fs*t_fill+1:fs*(t_fill+t_flat),1) = VF_scalar(2);
VF(fs*(t_fill+t_flat)+1:fs*(t_fill+t_flat+t_decay),1) = VF_scalar(3);

VF = VF + smooth(randn(length(VF),1)*kly_noise_std,1) + smooth(randn(length(VF),1)*kly_noise_std,1) * 1j;

% VP simulation without detuning
GS = 2*w12/(s+w12-1j*dw);
t = linspace(0,(t_fill+t_flat+t_decay),length(VF)).';
VP = lsim(GS,VF,t);
% plot linear simulation
figure(1)
clf
subplot(311)
plot(abs(VF),'LineWidth',2,'DisplayName','VF')
hold on, grid on
plot(abs(VP),'LineWidth',2,'DisplayName','VP')
legend('show')
subplot(312)
plot(180/pi*angle(VF),'LineWidth',2,'DisplayName','VF')
hold on, grid on
plot(180/pi*angle(VP),'--','LineWidth',2,'DisplayName','VP')
legend('show')
subplot(313)
hold on, grid on
legend('show')

% VP simulation with detuning --> time varying simulation
% Simulation shifted in time by starting at k = 2 --> signals aligned in
% time
dw = 2*pi*(100-1*abs(VP).^2);   % variation in detuning - approach questionable since mechanical modes neglected and only steady state considered...
VP(1) = 0;
for k = 2:length(VF)
    dVP(k) = -(w12+1j*dw(k))*VP(k-1) + 2*w12*VF(k);
    VP(k) = VP(k-1)+dVP(k)/fs;
end
% compute reflected signal as difference
VR = VP-VF;
% plot NL simulation
figure(1)
subplot(311)
plot(abs(VP),'LineWidth',2,'DisplayName','VP(dw)')
plot(abs(VR),'LineWidth',2,'DisplayName','VR')
title('Ideal Probe, Forward and Reflected')
subplot(312)
plot(180/pi*angle(VP),'LineWidth',2,'DisplayName','VP(dw)')
plot(180/pi*angle(VR),'LineWidth',2,'DisplayName','VR')
subplot(313)
plot(dw/2/pi,'LineWidth',2,'DisplayName','df')

% Now we have ideal VF and ideal VP --> move to realistic values
% Define klyston / FB noise
% Comment from Andrea:
% 1) if we are talking about klystron noise (equal to an input/plant disturbance), then such a noise effectively 'drives' the cavity.
%    Thereforeit has to be added only to the forward signal and the probe/reflected has to be calculated accordingly (in simulation)
% 2) the lines below only add in-phase noise. One should add the noise also in quadrature.
dist_FORW = smooth(randn(length(VF),1)*0.01,1).*0;
dist_REFL = smooth(randn(length(VF),1)*0.01,1).*0;
% define scaling factors by directivity mismatch
% paper: |a| = 0.976, |b| = 0.145, |c| = 0.207, |d| = 0.879
a = 0.976*exp(1i*10/180*pi);
b = 0.1545*exp(1i*-80/180*pi);
c = 0.207*exp(1i*-60/180*pi);
d = 0.879*exp(1i*-5/180*pi);
disp('Coupling parameters for simulation')
disp(['a - abs / phase: ',num2str(abs(a)),' / ',num2str(180/pi*angle(a)),'deg'])
disp(['b - abs / phase: ',num2str(abs(b)),' / ',num2str(180/pi*angle(b)),'deg'])
disp(['c - abs / phase: ',num2str(abs(c)),' / ',num2str(180/pi*angle(c)),'deg'])
disp(['d - abs / phase: ',num2str(abs(d)),' / ',num2str(180/pi*angle(d)),'deg'])
% define x and y coefficient (VP = x*VF + y*VR)
x = a+c;
y = b+d;
% compute measured 
VF_meas = (VF)./(x);    % VF = VF_meas * x  --> VF_meas = VF./x;
VR_meas = (VR)./(y);    % VR = VR_meas * y  --> VR_meas = VR./y;
% Compute measured VF and VR with realistic disturbances
% might be that forward / reflected port has difference characteristic...
VF_meas_a = (VF_meas.*a+dist_FORW);
VR_meas_b = (VR_meas.*b+dist_FORW);
VF_meas_c = (VF_meas.*c+dist_REFL);
VR_meas_d = (VR_meas.*d+dist_REFL);
% Compute disturbed common forw and refl. with correlated disturbances
VF_dist = VF_meas_a+VF_meas_c; 
VR_dist = VR_meas_b+VR_meas_d;
% use inverse of coupling matrix to map disturbed common forw and refl. to
% measured forw and refl.
coeff_inv = inv([a b ; c d]);
VP_meas_dist = VP + smooth(randn(length(VF),1)*meas_noise_std,1) + smooth(randn(length(VF),1)*meas_noise_std,1) * 1j
VF_meas_dist = VF_dist*coeff_inv(1,1) + VR_dist*coeff_inv(1,2) + smooth(randn(length(VF),1)*meas_noise_std,1) + smooth(randn(length(VF),1)*meas_noise_std,1) * 1j;
VR_meas_dist = VF_dist*coeff_inv(2,1) + VR_dist*coeff_inv(2,2) + smooth(randn(length(VF),1)*meas_noise_std,1) + smooth(randn(length(VF),1)*meas_noise_std,1) * 1j;
% plot measured signals
figure(2)
clf
subplot(321)
plot(abs(VP),'LineWidth',2,'DisplayName','VP(dw)')
hold on, grid on
plot(abs(VF_meas),'LineWidth',2,'DisplayName','VF')
plot(abs(VR_meas),'LineWidth',2,'DisplayName','VR')
title('Ideal Probe, Forward and Reflected')
subplot(323)
plot(180/pi*angle(VP),'LineWidth',2,'DisplayName','VP(dw)')
hold on, grid on
plot(180/pi*angle(VF_meas),'LineWidth',2,'DisplayName','VP(dw)')
plot(180/pi*angle(VR_meas),'LineWidth',2,'DisplayName','VR')
subplot(325)
plot(dw/2/pi,'LineWidth',2,'DisplayName','df')
hold on, grid on
subplot(322)
plot(abs(VP),'LineWidth',2,'DisplayName','VP(dw)')
hold on, grid on
plot(abs(VF_meas_dist),'LineWidth',2,'DisplayName','VF dist')
plot(abs(VR_meas_dist),'LineWidth',2,'DisplayName','VR dist')
title('Ideal Probe, Forward and Reflected')
subplot(324)
plot(180/pi*angle(VP),'LineWidth',2,'DisplayName','VP(dw)')
hold on, grid on
plot(180/pi*angle(VF_meas_dist),'LineWidth',2,'DisplayName','VF dist')
plot(180/pi*angle(VR_meas_dist),'LineWidth',2,'DisplayName','VR dist')
subplot(326)
plot(dw/2/pi,'LineWidth',2,'DisplayName','df')
hold on, grid on
% do first check and compute unique x and y with measurements
coeff = [VF_meas VR_meas]\VP;
disp(' * * * * * * * * * * * * * * * * * * * * * ')
disp(['x for simulation and its estimation: ', num2str(x),' & ',num2str(coeff(1))])
disp(['y for simulation and its estimation: ', num2str(y),' & ',num2str(coeff(2))])

% Calibration as in XFEL/FLASH
% this should follow the concept used for XFEL and FLASH ...
% 1 tuning parameter for removal of df and w12 discontinuity points
k_add = 2.69595*exp(1i*118.7687*pi/180);      % tuning parameter to be adjusted w.r.t. beam phase; leave it as 1 first
% define timing parameters for pulse and decay
t1 = (t_fill)*fs;
t2 = (t_fill+t_flat)*fs;
% Define offsets for decay phase
t2offsetstart = 100;
t2offsetend   = 250*10;

CAV_complex      = VP_meas_dist;
FOR_complex      = VF_meas_dist;
REFL_complex     = VR_meas_dist;

% Select signals for decay phase and compute coupling ratio
FOR_decay_compl  = FOR_complex(t2+t2offsetstart:t2+t2offsetend,1);
REFL_decay_compl = REFL_complex(t2+t2offsetstart:t2+t2offsetend,1);
CAV_decay_compl  = CAV_complex(t2+t2offsetstart:t2+t2offsetend,1);
% define A and B to have always same computation formula, see below
A = -REFL_decay_compl;
B = FOR_decay_compl;
% get dimension of decay coupling for weighting of coefficients
S = inv(A.'*A)*A.'*B;
% compute x (=a+c) and y (=b+d) for entire pulse
A = [FOR_complex,REFL_complex];
B = CAV_complex;
coeff_xy = inv(A.'*A)*A.'*B;
x_n = coeff_xy(1);
y_n = coeff_xy(2);
% select signals for coupling computation
a1 = FOR_complex;
a2 = REFL_complex;
a3 = FOR_decay_compl;
a4 = REFL_decay_compl;
b1 = CAV_complex;
b2 = CAV_decay_compl;
% define zeros for "big" matrix
zb = zeros(size(b2));
za = zeros(size(a3));
% put all vectors together to calibration in-/output matrix
Wb = abs(S);
Wc = k_add*abs(S);
% here A and B are given by Eqn. 10 in Paper
% with B = A * x with x =[a b c d]'
B = [b1;b2;zb;abs(x_n);abs(y_n)];
A = [a1,a2,a1,a2;...
    za,za,a3,a4;...
    a3,a4,za,za;...
    abs(x_n)-(Wc) 0 inv(Wc) 0;...
    0 inv(Wb) 0 abs(y_n)-(Wb);...
    ];
% compute coefficients
coeff = inv(A'*A)*A'*B;
% get a, b, c, d
disp(' * * * * * * * * * * * * * * * * * * * * * ')
disp('Cross-coupling results:')
a_n = coeff(1,:);
b_n = coeff(2,:);
c_n = coeff(3,:);
d_n = coeff(4,:);
a;
b;
c;
d;
% clc
disp(['a simulation / estimation: ',num2str(a),' / ',num2str(a_n)])
disp(['b simulation / estimation: ',num2str(b),' / ',num2str(b_n)])
disp(['c simulation / estimation: ',num2str(c),' / ',num2str(c_n)])
disp(['d simulation / estimation: ',num2str(d),' / ',num2str(d_n)])
disp(' * * * * * * * * * * * * * * * * * * * * * ')
disp(['Abs error in a: ',num2str(abs(a-a_n)/abs(x)*100),'% (normalized to a+c)'])
disp(['Abs error in b: ',num2str(abs(b-b_n)/abs(y)*100),'% (normalized to b+d)'])
disp(['Abs error in c: ',num2str(abs(c-c_n)/abs(x)*100),'% (normalized to a+c)'])
disp(['Abs error in d: ',num2str(abs(d-d_n)/abs(y)*100),'% (normalized to b+d)'])
disp(' * * * * * * * * * * * * * * * * * * * * * ')

FOR_complex_cal =  VF_meas_dist.*a_n+VR_meas_dist.*b_n;
REF_complex_cal =  VF_meas_dist.*c_n+VR_meas_dist.*d_n;

figure(1)
% clf
subplot(311)
plot(abs(VF_meas_dist.*a_n+VR_meas_dist.*b_n),'--','LineWidth',2,'DisplayName','VF cal')
hold on, grid on
plot(abs(VF_meas_dist.*c_n+VR_meas_dist.*d_n),'--','LineWidth',2,'DisplayName','VR cal')
legend('show')
ylabel('Amplitude [MV/m]')
xlabel(['Samples at ',num2str(fs/1e6),'MHz'])
subplot(312)
plot(180/pi*angle(VF_meas_dist.*a_n+VR_meas_dist.*b_n),'--','LineWidth',2,'DisplayName','VF cal')
hold on, grid on
plot(180/pi*angle(VF_meas_dist.*c_n+VR_meas_dist.*d_n),'--','LineWidth',2,'DisplayName','VR cal')
legend('show')
ylabel('Phase [deg]')
xlabel(['Samples at ',num2str(fs/1e6),'MHz'])
subplot(313)
hold on, grid on
title(['With k_{add} = ',num2str(abs(k_add)),'* exp(1i*',num2str(angle(k_add)*180/pi),'deg)'])
legend('show')

%
dt = 1/fs ; 
smooth_span = round(0.1e-6*fs);
K = w12(1);
disp(['Smooth dV over ', num2str(smooth_span/fs*1e6),'us for dw and w12 computation'])
%disp(['Smooth dV via ', num2str(smooth_span),' samples'])
if(1)
    dVI = [0 ; smooth(diff(real(CAV_complex)),smooth_span)/dt];
    dVQ = [0 ; smooth(diff(imag(CAV_complex)),smooth_span)/dt];
else
    dVI = [0 ; diff(smooth(real(CAV_complex),smooth_span))/dt];
    dVQ = [0 ; diff(smooth(imag(CAV_complex),smooth_span))/dt];
end

dw1 = real(CAV_complex).*(dVQ - 2*K*imag(FOR_complex_cal))./abs(CAV_complex).^2;
dw2 = imag(CAV_complex).*(2*K*real(FOR_complex_cal) - dVI)./abs(CAV_complex).^2;

% dw_1EQN = imag(-(K*2*FOR_complex_cal.*(CAV_complex').' - complex(dVI,dVQ).*(CAV_complex').')./abs(CAV_complex).^2);

dw_est = (dw1 + dw2);
df_est = -dw_est./2/pi;
w12_1 = real(CAV_complex).*(2*K*real(FOR_complex_cal)-dVI)./abs(CAV_complex).^2;
w12_2 = imag(CAV_complex).*(2*K*imag(FOR_complex_cal)-dVQ)./abs(CAV_complex).^2;

w12_est = w12_1 + w12_2;

figure(1)
subplot(313)
plot(movmean(df_est,5),'LineWidth',2,'DisplayName','df est')
plot(movmean(w12_est/2/pi-141,50),'LineWidth',2,'DisplayName','f12 est - 141Hz')
ylabel('Freq [Hz]')
xlabel(['Samples at ',num2str(fs/1e6),'MHz'])
xlim([1000 length(VF)])

disp('Simulation finished')

% Save simulation dataset
save('dataset.mat')
