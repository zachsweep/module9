%% Problem 1
clear;
% close all;

% M bands
M = 20;

% M band DFT filter bank
l = (0:M-1);
k = (0:M-1).';
b_dft = exp(1j*2*pi*l.*k/M);

% % Compute frequency response for first 5 bands
% H0 = (1/M)*fft(b_dft(1,:),2048);
% H1 = (1/M)*fft(b_dft(2,:),2048);
% H2 = (1/M)*fft(b_dft(3,:),2048);
% H3 = (1/M)*fft(b_dft(4,:),2048);
% H4 = (1/M)*fft(b_dft(5,:),2048);
% 
% % Plot magnitude response for first 5 bands
% figure(1)
% hold on
% w = (0:2047)*2*pi/2048;
% plot(w/pi,20*log10(abs(H0)));
% plot(w/pi,20*log10(abs(H1)));
% plot(w/pi,20*log10(abs(H2)));
% plot(w/pi,20*log10(abs(H3)));
% plot(w/pi,20*log10(abs(H4)));
% axis([0 2 -60 10]);
% title({'20-Band DFT Filter Banks','Banks(1-5)'});
% xlabel('Normalized Frequency (x pi rad/sample)');ylabel('Magnitude (dB)');


%% Problem 2

% M-bands
M = 20;

% Prototype filter passband and stopband ripple specificaitons
% Rp is assumed to be the peak maximum passband ripple in dB
Rp = 0.1;   %dB
Rs =  40;   %dB

% Find linear value of Rp
Rp_linear = min([(1-10^(-Rp/20)) (10^(Rp/20)-1)]);

% Choose cutoff such that each sub-band overlaps the neighbor at 3dB.
% This required iterating through different transition band widths as well
% as the the transition band center
% wp = 0.048*pi;
% ws = 0.0565*pi;
offset = pi/50;
cent = (pi/M) +0.0196;
wp = cent - offset;
ws = cent + offset;

% Compute order of prototype filter needed to meet specifications
[n,fo,mo,w] = firpmord([wp/pi ws/pi], [1 0], [Rp_linear 10^(-Rs/20)]);

% Order was surprisingly over estimated for this filter
n = n-1;
pm_order = n;
fprintf('\nParks McClellan Prototype Filter Order: %d\n',pm_order);

% Compute prototype filter impulse response
b = firpm(n,fo,mo,w);

% Compute the subbands by modulating the prototype filter impulse response
b_pm = zeros(M,n+1);
for i=0:M-1
    b_pm(i+1,:) = b.*exp(1j*2*pi*(0:length(b)-1)*i/M);
end

% Compute FFT for first 5 sub-bands
% H1 = fft(b_pm(1,:),2048);
% H2 = fft(b_pm(2,:),2048);
% H3 = fft(b_pm(3,:),2048);
% H4 = fft(b_pm(4,:),2048);
% H5 = fft(b_pm(5,:),2048);
% 
% w = (0:2047)*2*pi/2048;
% 
% figure(2)
% hold on;
% plot(w/pi,20*log10(abs(H1)));
% plot(w/pi,20*log10(abs(H2)));
% plot(w/pi,20*log10(abs(H3)));
% plot(w/pi,20*log10(abs(H4)));
% plot(w/pi,20*log10(abs(H5)));
% line([0 2],[0.1 0.1],'color','red','LineStyle','--');
% line([0 2],[-0.1 -0.1],'color','red','LineStyle','--');
% line([0 2],[-3 -3],'color','red','LineStyle','--');
% axis([0 2 -60 10]);
% title({'20-Band Parks McClellan Filter Banks','Banks(1-5)'});
% xlabel('Normalized Frequency (x pi rad/sample)');ylabel('Magnitude (dB)');
% 
% % Zoomed view of passband
% figure(3)
% hold on;
% plot(w/pi,20*log10(abs(H1)));
% plot(w/pi,20*log10(abs(H2)));
% plot(w/pi,20*log10(abs(H3)));
% plot(w/pi,20*log10(abs(H4)));
% plot(w/pi,20*log10(abs(H5)));
% line([0 2],[0.1 0.1],'color','red','LineStyle','--');
% line([0 2],[-0.1 -0.1],'color','red','LineStyle','--');
% line([0 2],[-3 -3],'color','red','LineStyle','--');
% axis([0 0.045 -0.2 0.15]);
% title({'20-Band Parks McClellan Filter Banks','Zoomed Passband bank 1'});
% xlabel('Normalized Frequency (x pi rad/sample)');ylabel('Magnitude (dB)');
% 
% % Zoomed view of overlap
% figure(4)
% hold on;
% plot(w/pi,20*log10(abs(H1)));
% plot(w/pi,20*log10(abs(H2)));
% plot(w/pi,20*log10(abs(H3)));
% plot(w/pi,20*log10(abs(H4)));
% plot(w/pi,20*log10(abs(H5)));
% line([0 2],[0.1 0.1],'color','red','LineStyle','--');
% line([0 2],[-0.1 -0.1],'color','red','LineStyle','--');
% line([0 2],[-3 -3],'color','red','LineStyle','--');
% axis([0.049 0.0505 -3.08 -2.92]);
% title({'20-Band Parks McClellan Filter Banks','Zoomed view of overlap'});
% xlabel('Normalized Frequency (x pi rad/sample)');ylabel('Magnitude (dB)');
% 
% % Zoomed view of stopband
% figure(5)
% hold on;
% plot(w/pi,20*log10(abs(H1)));
% plot(w/pi,20*log10(abs(H2)));
% plot(w/pi,20*log10(abs(H3)));
% plot(w/pi,20*log10(abs(H4)));
% plot(w/pi,20*log10(abs(H5)));
% line([0 2],[0.1 0.1],'color','red','LineStyle','--');
% line([0 2],[-0.1 -0.1],'color','red','LineStyle','--');
% line([0 2],[-3 -3],'color','red','LineStyle','--');
% line([0 2],[-40 -40],'color','red','LineStyle','--');
% axis([0 0.55 -46 -35]);
% title({'20-Band Parks McClellan Filter Banks','Zoomed Stopband'});
% xlabel('Normalized Frequency (x pi rad/sample)');ylabel('Magnitude (dB)');


%% Problem 3

% Number of banks in filter bank
MM = 20;

% Choose cutoff such that each sub-band overlaps the neighbor at 3dB.
% The design of the filter that met the given specifications required
% iterating throught the design process with different values for the
% passband/stopband frequencies such that (wp + ws)/2 = pi/M was satisfied.
% offset = pi/31.4;
% offset = pi/19.86;
% wp = pi/MM - offset;
% ws = pi/MM + offset;

% This is apparently the new solution  b
offset = pi/24.5;
cent = (pi/MM) + 0.0196;
wp = cent - offset;
ws = cent + offset;

N = pm_order;

% Since the order of the prototype filter in problem 2 has an even order 
% (impulse response has odd length) and even symetry this is a type II
% linear phase FIR filter.
M = N/2;

% Stopband weight alpha = 0.2
alpha = 0.2;

% Pre-allocating impulse response vector
h = zeros(1,N+1);

% Compute c(w) and c(w)*c(w).'
syms w;
cw = cos((0:M)*w).';
cw = cw*cw.';

% Compute Ps Matrix
Ps = zeros(M+1,M+1);
w = linspace(ws,pi,1000);

for m=0:M
    for n=0:M
        % Evaluate cw(m,n) at each value of w
        c = eval(cw(m+1,n+1));

        % Handling cases where c = 1 or 0
        if c==1
            c = ones(length(w),1);
        elseif c==0
            c = zeros(length(w),1);
        end

        % Compute integral
        Ps(m+1,n+1) = (1/pi)*trapz(w,c);
    end
end

% Compute c(w) and (1-c(w))*(1-c(w)).'
syms w;
cw = cos((0:M)*w).';
cw = (1 - cw)*(1 - cw).';

% Compute Pp Matrix
Pp = zeros(M+1,M+1);
w = linspace(0,wp,1000);

for m=0:M
    for n=0:M
        % Evaluate cw(m,n) at each value of w
        c = eval(cw(m+1,n+1));

        % Handling cases where c = 1 or 0
        if c==1
            c = ones(length(w),1);    
        elseif c==0
            c = zeros(length(w),1);
        end

        % Compute Integral
        Pp(m+1,n+1) = (1/pi)*trapz(w,c);
    end
end

% Compute P matrix with weights of alpha = 0.2 for the stop band and
% (1-alpha) for the passband
P = alpha*Ps + (1-alpha)*Pp;

% Compute Eigen Vectors/Values of P
[V,D] = eig(P,'vector');

% Find index of smallest Eigen value in the Eigen value column vector
ind = find(D==min(D));

% Find Eigen Vector containing smallest Eigen value using the index
b = V(:,ind);

% Re-organize bn to get h(n)
% h(M) = b(0)
h(M+1) = b(1);

% h(n) = b(n)/2 for n = 1 to M
h(M+2:end) = b(2:M+1).'/2;
h(1:M) = flip(b(2:M+1).'/2);

% Normalizing the impulse response such that the gain of the filter is
% unity
h = h/sum(h);

% Create the 20-band filter bank by modulating the eigen filter prototype
b_eig = zeros(MM,N+1);
for i=0:MM-1
    b_eig(i+1,:) = h.*exp(1j*2*pi*(0:length(h)-1)*i/MM);
end

% Compute FFT for first 5 sub-bands
H1 = fft(b_eig(1,:),8192);
H2 = fft(b_eig(2,:),8192);
H3 = fft(b_eig(3,:),8192);
H4 = fft(b_eig(4,:),8192);
H5 = fft(b_eig(5,:),8192);

w = (0:8191)*2*pi/8192;

% Plot magnitude response for first 5 bands
figure(6)
hold on;
plot(w/pi,20*log10(abs(H1)));
plot(w/pi,20*log10(abs(H2)));
plot(w/pi,20*log10(abs(H3)));
plot(w/pi,20*log10(abs(H4)));
plot(w/pi,20*log10(abs(H5)));
line([0 2],[0.1 0.1],'color','red','LineStyle','--');
line([0 2],[-0.1 -0.1],'color','red','LineStyle','--');
line([0 2],[-3 -3],'color','red','LineStyle','--');
axis([0 2 -60 10]);
title({'20-Band Eigen Filter Banks','Banks(1-5)'});
xlabel('Normalized Frequency (x pi rad/sample)');ylabel('Magnitude (dB)');

% Plot zoomed view first 5 bands
figure(7)
hold on;
plot(w/pi,20*log10(abs(H1)));
plot(w/pi,20*log10(abs(H2)));
plot(w/pi,20*log10(abs(H3)));
plot(w/pi,20*log10(abs(H4)));
plot(w/pi,20*log10(abs(H5)));
line([0 2],[0.1 0.1],'color','red','LineStyle','--');
line([0 2],[-0.1 -0.1],'color','red','LineStyle','--');
line([0 2],[-3 -3],'color','red','LineStyle','--');
axis([0 0.25 -60 10]);
title({'20-Band Eigen Filter Banks','zoomed view of bands 1-5'});
xlabel('Normalized Frequency (x pi rad/sample)');ylabel('Magnitude (dB)');

% % Plot passband at overlap
% figure(8)
% hold on;
% plot(w/pi,20*log10(abs(H1)));
% plot(w/pi,20*log10(abs(H2)));
% plot(w/pi,20*log10(abs(H3)));
% plot(w/pi,20*log10(abs(H4)));
% plot(w/pi,20*log10(abs(H5)));
% line([0 2],[0.1 0.1],'color','red','LineStyle','--');
% line([0 2],[-0.1 -0.1],'color','red','LineStyle','--');
% line([0 2],[-3 -3],'color','red','LineStyle','--');
% axis([0 0.025 -3.5 -2.5]);
% title({'20-Band Eigen Filter Banks','zoomed view band overlap'});
% xlabel('Normalized Frequency (x pi rad/sample)');ylabel('Magnitude (dB)');

%% Problem 4
% [y,fs] = audioread('Misty Mountain Hop Snippet.wav');
% 
% len = length(y);
% 
% f1 = 'recordings/dft_bank_';
% f2 = 'recordings/pm_bank_';
% f3 = 'recordings/eig_bank_';
% 
% % Normalize the gain of the dft filter bank impulse responses
% b_dft = b_dft/sum(b_dft(:));
% 
% for i=1:5
%     figure(i+10)
%     hold on
%     
%     % Apply dft filter bank(i) to input signal y. Normalize the filtered
%     % signal between -1 and 1. Write filtered signal to file
%     fname = strcat(f1,int2str(i),'.wav');
%     y_dft = real(conv(y,b_dft(i,:),'same'));
%     y_dft = y_dft./(max(abs(y_dft)));
%     audiowrite(fname,y_dft,fs);
%     
%     % Apply Parks McClellan filter bank(i) to input signal y. Normalize the
%     % filtered signal between -1 and 1. Write filtered signal to file
%     fname = strcat(f2,int2str(i),'.wav');
%     y_pm = real(conv(y,b_pm(i,:),'same'));
%     y_pm = y_pm./(max(abs(y_pm)));
%     audiowrite(fname,y_pm,fs);
%     
%     % Apply Eigen filter bank(i) to input signal y. Normalize the filtered 
%     % signal between -1 and 1. Write filtered signal to file
%     fname = strcat(f3,int2str(i),'.wav');
%     y_eig = real(conv(y,b_eig(i,:),'same'));
%     y_eig = y_eig./(max(abs(y_eig)));
%     audiowrite(fname,y_eig,fs);
%     
%     % Compute FFT for each filter bank output
%     H_dft = fft(y_dft);
%     H_pm = fft(y_pm);
%     H_eig = fft(y_eig);
%     
%     % Compute power spectral density of signal after
%     
%     PSD_dft = H_dft.*conj(H_dft);
%     plot(PSD_dft(1:len/2));
%     
%     % Compute power spectral density of signal after applying pm filter
%     
%     PSD_pm = H_pm.*conj(H_pm);
%     plot(PSD_pm(1:len/2));
%     
%     % Compute power spectral density of signal after applying eig filter
%     
%     PSD_eig = H_eig.*conj(H_eig);
%     plot(PSD_eig(1:len/2));
%     
%     hold off
%     text = strcat('Power Spectral Density Bank ',int2str(i));
%     title(text);
% end

