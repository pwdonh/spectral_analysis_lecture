
% cd /home/peterd/Documents/MEG_workshop/

sFile = './matrix_150203_1839.mat';
sMat = load(sFile);
trial = 70;
epoch = [2, 6];
fs = 1/(sMat.Time(2)-sMat.Time(1));
iEpoch = (sMat.Time>epoch(1))&(sMat.Time<epoch(2));
signal = zscore(sMat.Value(trial,iEpoch));
t = sMat.Time(iEpoch);
iBase = t>5;
iStim = t<3.5;
tp = t-epoch(1);

%%

set(groot,'defaultAxesColorOrder',[0 0.4470 0.7410; 0.8500 0.3250 0.0980]);
set(groot,'defaultAxesColorOrder',[0 0.4470 0.7410]);

%% Linear vs. log-scale

frequencies = linspace(1,20,5);
frequencies2 = logspace(log10(1),log10(20),5);
figure(1); clf(); hold on
set(gcf, 'PaperPosition', [0 0 4 3]);
figure(2); clf(); hold on
set(gcf, 'PaperPosition', [0 0 4 3]);
for iF = 1:length(frequencies)
    figure(1)
    y{iF} = sin(2*pi*frequencies(iF)*tp(1:3:end));
    plot(tp(1:3:end), y{iF}-iF*3)
    figure(2)
    y{iF} = sin(2*pi*frequencies2(iF)*tp(1:3:end));
    plot(tp(1:3:end), y{iF}-iF*3)
end
% figure(1)
% print ./sinusoids_linear.eps -depsc2
% figure(2)
% print ./sinusoids_log.eps -depsc2

%% Fourier series

figure(1); clf(); hold on
set(gcf, 'PaperPosition', [0 0 4 3]);

plot(tp(1:3:end), y{1}+10)
plot(tp(1:3:end), y{2}+4)
plot(tp(1:3:end), y{3}-2)
plot(tp(1:3:end), (y{1}+y{2}+y{3})./3-8)

% print ./fourier_series.eps -depsc2

%% Time -> Frequency

yseries = y{1}+y{2}+y{3};
yseries2 = yseries(1:2:end);
yseries3 = yseries(1:4:end);
yseries4 = yseries(1:8:end);
fs2 = fs/3/2; fs3 = fs/3/4; fs4 = fs/3/8;

figure(1); clf(); hold on
set(gcf, 'PaperPosition', [0 0 4 3]);

plot(tp(1:3:end), yseries+10, 'o-', 'MarkerSize', 2)
plot(tp(1:3*2:end), yseries2+4, 'o-', 'MarkerSize', 2)
plot(tp(1:3*4:end), yseries3-2, 'o-', 'MarkerSize', 2)
plot(tp(1:3*8:end), yseries4-8, 'o-', 'MarkerSize', 2)

% print ./fourier_series_sampling.eps -depsc2

y_power2 = (abs(fft(yseries2)).^2)./length(yseries2);
y_power3 = (abs(fft(yseries3)).^2)./length(yseries2);
y_power4 = (abs(fft(yseries4)).^2)./length(yseries2);

figure(2); clf(); hold on
set(gcf, 'PaperPosition', [0 0 4 3]);

cla; plot(0:fs4/length(y_power4):fs4-.1, y_power4, 'o-', 'MarkerSize', 2);
ylim([-1 20]);
% print ./fourier_series_fft4full.eps -depsc2
figure(3); clf(); hold on
cla; plot(0:fs4/length(y_power4):fs4/2, y_power4(1:length(y_power4)/2+1)*4);
% print ./fourier_series_fft4.eps -depsc2
figure(4); clf(); hold on
cla; plot(0:fs3/length(y_power3):fs3-.1, y_power3/4, 'o-', 'MarkerSize', 2);
% print ./fourier_series_fft3full.eps -depsc2
figure(5); clf(); hold on
cla; plot(0:fs3/length(y_power3):fs3/2, y_power3(1:length(y_power3)/2+1));
% print ./fourier_series_fft3.eps -depsc2
figure(6); clf(); hold on
cla; plot(0:fs2/length(y_power2):fs2/2, y_power2(1:length(y_power2)/2+1)/4);
% print ./fourier_series_fft2.eps -depsc2

%% Raw signal

figure(1); clf; hold on
set(gcf, 'PaperPosition', [0 0 4 3]);
plot(t, signal, 'Color', [0 0.4470 0.7410]);
% xlim(epoch);
% plot([0 0], [-10 -5]);
plot([3.795 3.795], [-10 -5]);
% print ./signal_raw.eps -depsc2

%% FFT

N = length(signal);
signal_fft = fft(signal);
signal_fft = signal_fft(1:N/2+1);
f_fft = 0:fs/N:fs/2;
signal_psd_fft = (1/(fs*N))*abs(signal_fft).^2;
signal_psd_fft(2:end-1) = 2*signal_psd_fft(2:end-1);

figure(1); clf; hold on;
set(gcf,'PaperPosition',[0 0 4 3]);
plot(f_fft, signal_psd_fft, 'Color', [0 0.4470 0.7410]);
xlabel('Frequency (Hz)');
ylabel('Power');

% print ./signal_fft.eps -depsc2

set(gca,'XScale','log',...
    'XTick', [1.5, 3, 6, 13, 30, 70, 150, 300],...
    'XLim', [1, 300], 'LineWidth',1.5);
% print ./signal_fft_loglinear.eps -depsc2

set(gca,'YScale','log');
% print ./signal_fft_loglog.eps -depsc2

%% PSD

lwindow = [.1 .25 .5 1 2];

for iWin = 1:length(lwindow)
    [signal_psd{iWin}, f_psd{iWin}] =  bst_psd(signal, fs, lwindow(iWin), .75);
end

for iWin = 1:length(lwindow)
    figure(iWin); clf; hold on;
    set(gcf,'PaperPosition',[0 0 4 3]);
    
    l=plot(f_psd{iWin}, squeeze(signal_psd{iWin}));
    set(l, 'Color', [0 0.4470 0.7410], 'LineWidth', 1.5);
    set(gca,'XLim',[1 300]);
    xlabel('Frequency (Hz)');
    ylabel('Power');
    
%     print(['./signal_psd_linear' num2str(lwindow(iWin)) '.eps'], '-depsc2');
    
    l=plot(f_psd{3}, squeeze(signal_psd{3}));
    set(l, 'Color', [0 0.4470 0.7410], 'LineWidth', 1.5);
    set(gca,'XLim',[1 300], 'YLim', [10e-6 10e-1]);
    set(gca,'XScale','log',...
        'XTick', [1.5, 3, 6, 13, 30, 70, 150, 300]);
%     print(['./signal_psd_loglinear' num2str(lwindow(iWin)) '.eps'], '-depsc2');
    set(gca,'YScale','log');
%     print(['./signal_psd_loglog' num2str(lwindow(iWin)) '.eps'], '-depsc2');
end

%% Wavelets

f = logspace(log10(2), log10(150), 150);
w = morlet_transform(signal, t, f, 1, 4);
w_zbase = bsxfun(@rdivide, bsxfun(@minus, w, mean(w(1,iBase,:))), std(w(1,iBase,:)));
w_zstim = bsxfun(@rdivide, bsxfun(@minus, w, mean(w(1,iStim,:))), std(w(1,iStim,:)));
w_dbbase = 10.^bsxfun(@minus, log10(w), log10(mean(w(1,iBase,:))));

set(gcf,'PaperPosition',[0 0 4 3]);
figure(1); clf;
pcolor(tp, f, squeeze(w)'); shading flat
set(gca, 'YScale', 'log',...
    'YTick', [3, 6, 13, 30, 70, 150]);
colormap('jet'); cb = colorbar;
ylabel(cb, 'power');
xlabel('Time (seconds)');
ylabel('Frequency (Hz)');
% print ./signal_tfmap.png -dpng

figure(2); clf;
pcolor(tp, f, squeeze(w_zbase)'); shading flat
set(gca, 'YScale', 'log', 'CLim', [-20, 20],...
    'YTick', [3, 6, 13, 30, 70, 150]);
colormap('jet'); cb = colorbar;
ylabel(cb, 'z-score');
xlabel('Time (seconds)');
ylabel('Frequency (Hz)');
% print ./signal_tfmap_zstim.png -dpng

figure(3); clf;
pcolor(tp, f, squeeze(w_zstim)'); shading flat
set(gca, 'YScale', 'log', 'CLim', [-20, 20],...
    'YTick', [3, 6, 13, 30, 70, 150]);
colormap('jet'); cb = colorbar;
ylabel(cb, 'z-score');
xlabel('Time (seconds)');
ylabel('Frequency (Hz)');
% print ./signal_tfmap_zbase.png -dpng

figure(4); clf;
plot(f, squeeze(mean(w)), 'Color', [0 0.4470 0.7410], 'LineWidth', 1.5);
set(gca, 'XScale', 'log', 'YScale', 'log',...
    'XTick', [3, 6, 13, 30, 70, 150]);
xlim([2, 150]);
ylim([2e-5 .5e0])
% print ./signal_psd_wavelets.eps -depsc2

figure(5); clf;
plot(tp, squeeze(w(1,:,53)))
ylim([-.1 .6])
% print ./signal_wavelet_ts.eps -depsc2

figure(6); clf;
plot(tp, squeeze(w(1,:,125)))
ylim([-.001 .003])
% print ./signal_wavelet_ts_gamma.eps -depsc2

%% Wavelet different parameters

f = logspace(log10(2), log10(150), 150);
w = morlet_transform(signal, t, f, 1, 2);
w_zbase = bsxfun(@rdivide, bsxfun(@minus, w, mean(w(1,iBase,:))), std(w(1,iBase,:)));

figure(1); clf;
pcolor(tp, f, squeeze(w_zbase)'); shading flat
set(gca, 'YScale', 'log', 'CLim', [-20, 20],...
    'YTick', [3, 6, 13, 30, 70, 150]);
colormap('jet'); cb = colorbar;
ylabel(cb, 'z-score');
xlabel('Time (seconds)');
ylabel('Frequency (Hz)');
% print ./signal_tfmap2_zbase.png -dpng

%% Plot morlet wavelet

tw = -3:1/fs:3;
wavel = morlet_wavelet(tw, 1, 1);

figure(1); clf; hold on
set(gcf, 'PaperPosition', [0 0 4 3]);
plot(tw, real(wavel), 'Color', [0 0.4470 0.7410]);
% print ./morlet_wavelet.eps -depsc2

wavel = morlet_wavelet(tw, 1, .6);

figure(2); clf; hold on
set(gcf, 'PaperPosition', [0 0 4 3]);
plot(tw, real(wavel), 'Color', [0 0.4470 0.7410]);
% print ./morlet_wavelet2.eps -depsc2

%% PSD wavelet over bins

fbands = [2 4; 4 8; 8 12; 13 30; 30 60; 60 90];
for iBand = 1:size(fbands,1)
    findex = (f>=fbands(iBand,1))&(f<=fbands(iBand,2));
    w_band(iBand) = mean(mean(squeeze(w(:,:,findex))));
end

figure(1); clf;
set(gcf,'PaperPosition',[0 0 4 3]);
plot(mean(fbands,2), w_band, 'o-', 'Color', [0 0.4470 0.7410], 'LineWidth', 1.5);
plot(1:size(fbands,1), w_band, 'o-', 'Color', [0 0.4470 0.7410], 'LineWidth', 1.5);
% bar(1:size(fbands,1), w_band);
set(gca, 'YScale', 'log',...
    'XTick', 1:size(fbands,1), 'XTickLabel',{'delta','theta',...
    'alpha','beta','gamma1','gamma2'});
% set(gca, 'YScale', 'log');
xlim([1-.15, size(fbands,1)+.15]);
ylim([2e-5 .5e0])

% print ./signal_psd_wavelets_bar.eps -depsc2

%% Hilbert

signal_delta = process_bandpass('Compute', signal, fs, 2.6, 3.8);
signal_alpha = process_bandpass('Compute', signal, fs, 8.5, 10.5);
signal_gamma = process_bandpass('Compute', signal, fs, 60, 90);
env_delta = abs(hilbert(signal_delta')');
env_alpha = abs(hilbert(signal_alpha')');
env_gamma = abs(hilbert(signal_gamma')');
phase_delta = angle(hilbert(signal_delta')');
phase_alpha = angle(hilbert(signal_alpha')');
phase_gamma = angle(hilbert(signal_gamma')');

figure(1); clf; hold on; set(gcf, 'PaperPosition', [0 0 4 3]);
plot(tp, signal+8.5, 'Color', [0 0.4470 0.7410]);
plot(tp, signal_delta*1.5+2.5, 'Color', [0 0.4470 0.7410]);
plot(tp, signal_alpha*1.5-2.5, 'Color', [0 0.4470 0.7410]);
plot(tp, signal_gamma*1.5-7.5, 'Color', [0 0.4470 0.7410]);

ylim([-10, 14]);
% print ./signal_hilbert_filt.eps -depsc2

plot(tp, env_delta*1.5+2.5, 'LineWidth', 1.25, 'Color', [0.8500 0.3250 0.0980]);
plot(tp, env_alpha*1.5-2.5, 'LineWidth', 1.25, 'Color', [0.8500 0.3250 0.0980]);
plot(tp, env_gamma*1.5-7.5, 'LineWidth', 1.25, 'Color', [0.8500 0.3250 0.0980]);

ylim([-10, 14]);
% print ./signal_hilbert_amp.eps -depsc2

figure(2); clf; hold on; set(gcf, 'PaperPosition', [0 0 4 3]);
plot(tp, signal+8.5, 'Color', [0 0.4470 0.7410]);
plot(tp, phase_delta./1.75+2.5, 'Color', [0 0.4470 0.7410]);
plot(tp, phase_alpha./1.75-2.5, 'Color', [0 0.4470 0.7410]);
plot(tp, phase_gamma./1.75-7.5, 'Color', [0 0.4470 0.7410]);
ylim([-10, 14]);
% print ./signal_hilbert_phase.eps -depsc2

cla
plot(tp, real(hilbert(signal_delta')')*1.5+2.5, 'Color', [0 0.4470 0.7410]);
plot(tp, imag(hilbert(signal_delta')')*1.5+2.5, 'Color', [0.8500    0.3250    0.0980]);

plot(tp, env_delta*1.5+2.5, 'LineWidth', 1.25, 'Color', [0 0.4470 0.7410]);
plot(tp, phase_delta./2.25-1.5, 'Color', [0 0.4470 0.7410]);






