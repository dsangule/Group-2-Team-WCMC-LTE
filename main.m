%% Parameters
% LTE Parameters
bandwidth = 1.4e6;          % Bandwidth in Hz
Nrb = 6;                    % Number of resource blocks (1.4 MHz)
Nsc = 12;                   % Subcarriers per RB
Nfft = 128;                 % FFT size (for 1.4 MHz)
Ncp = 9;                    % Cyclic prefix length (~4.7 us)
numSymbols = 14;            % OFDM symbols per subframe (2 slots)
Fs = 1.92e6;                % Sampling rate
subframeDuration = 1e-3;    % 1 ms
Nsubcarriers = Nrb * Nsc;   % Total subcarriers (72)
SNR_dB = 10;                % SNR for AWGN channel
useFading = true;           % Enable multipath fading

% Cell ID Parameters
N_ID_2 = 0;                 % PSS index (0, 1, or 2)
N_ID_1 = 0;                 % SSS index (0 to 167)
PCI = 3 * N_ID_1 + N_ID_2;  % Physical Cell Identity

%% Generate PSS (Zadoff-Chu Sequence)
u = [25, 29, 34];          % Root indices for PSS (N_ID_2 = 0, 1, 2)
pssLen = 62;                % PSS sequence length
pssSeq = zeros(pssLen, 1);
for n = 0:pssLen-1
    pssSeq(n+1) = exp(-1i * pi * u(N_ID_2+1) * n * (n+1) / 63);
end
% Map PSS to 72 subcarriers (centered, DC unused)
pssGrid = zeros(Nsubcarriers, 1);
pssGrid(6:67) = pssSeq;     % Center 62 subcarriers

%% Generate Time-Domain PSS
pssFreq = zeros(Nfft, 1);
startIdx = (Nfft - Nsubcarriers) / 2 + 1; % = 29
pssFreq(startIdx:startIdx+Nsubcarriers-1) = pssGrid;
pssTime = ifft(pssFreq, Nfft);
pssTimeCP = [pssTime(end-Ncp+1:end); pssTime]; % CP + symbol

%% Generate SSS (Placeholder)
sssLen = 62;
sssSeq = exp(1i * pi * (0:sssLen-1)' / sssLen); % Simplified SSS
sssGrid = zeros(Nsubcarriers, 1);
sssGrid(6:67) = sssSeq;
% Time-domain SSS
sssFreq = zeros(Nfft, 1);
sssFreq(startIdx:startIdx+Nsubcarriers-1) = sssGrid;
sssTime = ifft(sssFreq, Nfft);
sssTimeCP = [sssTime(end-Ncp+1:end); sssTime];

%% Generate CRS (Simplified QPSK)
crsSeq = (randi([0 1], Nsubcarriers, 1) * 2 - 1) / sqrt(2);
crsPos = [1, 5, 8, 12];    % Symbols 0, 4, 7, 11 (0-based)

%% Generate Resource Grid (1 subframe)
resourceGrid = zeros(Nsubcarriers, numSymbols);
% PSS (symbol 6, slot 1)
resourceGrid(:, 7) = pssGrid;
% SSS (symbol 5, slot 1)
resourceGrid(:, 6) = sssGrid;
% CRS
for sym = crsPos
    resourceGrid(:, sym+1) = crsSeq;
end
% PBCH (symbols 7-10)
pbchData = (randi([0 1], Nsubcarriers*4, 1) * 2 - 1) / sqrt(2);
pbchGrid = reshape(pbchData, Nsubcarriers, 4);
resourceGrid(:, 8:11) = pbchGrid;
% PDSCH (remaining REs)
for sym = 1:numSymbols
    if ~ismember(sym, [6, 7, 8, 9, 10, 11, 12])
        resourceGrid(:, sym) = (randi([0 1], Nsubcarriers, 1) * 2 - 1) / sqrt(2);
    end
end

%% OFDM Modulation
txSubframe = zeros(Nfft + Ncp, numSymbols);
for sym = 1:numSymbols
    freqSym = zeros(Nfft, 1);
    freqSym(startIdx:startIdx+Nsubcarriers-1) = resourceGrid(:, sym);
    timeSym = ifft(freqSym, Nfft);
    txSubframe(:, sym) = [timeSym(end-Ncp+1:end); timeSym];
end
txFrame = txSubframe(:);

%% Channel Model
% AWGN
noisePower = 10^(-SNR_dB/10);
noise = sqrt(noisePower/2) * (randn(size(txFrame)) + 1i * randn(size(txFrame)));
rxFrame = txFrame + noise;

% Multipath Fading
if useFading
    tau = [0 1e-6]; % Delay spread
    pdb = [0 -3];   % Power delay profile
    fading = zeros(size(txFrame));
    for tap = 1:length(tau)
        gain = 10^(pdb(tap)/20) * (randn(size(txFrame)) + 1i * randn(size(txFrame))) / sqrt(2);
        delaySamples = round(tau(tap) * Fs);
        fading = fading + [zeros(delaySamples, 1); gain(1:end-delaySamples)];
    end
    rxFrame = conv(txFrame, fading, 'same') + noise;
end

%% CFO Estimation
cfoCorr = 0;
numSyms = floor(length(rxFrame) / (Nfft + Ncp));
for sym = 1:numSyms
    startIdx = (sym-1) * (Nfft + Ncp) + 1;
    cp = rxFrame(startIdx:startIdx+Ncp-1);
    tail = rxFrame(startIdx+Nfft:startIdx+Nfft+Ncp-1);
    cfoCorr = cfoCorr + sum(conj(cp) .* tail);
end
cfoEst = angle(cfoCorr) / (2 * pi * Nfft / Fs);
t = (0:length(rxFrame)-1)' / Fs;
rxFrame = rxFrame .* exp(-1i * 2 * pi * cfoEst * t);

%% PSS Detection (Time-Domain)
correlationLength = length(rxFrame) - length(pssTimeCP) + 1;
correlation = zeros(correlationLength, 1);
for i = 1:correlationLength
    segment = rxFrame(i:i+length(pssTimeCP)-1);
    if length(segment) == length(pssTimeCP)
        correlation(i) = abs(sum(conj(pssTimeCP) .* segment));
    end
end
[peakVal, peakPos] = max(correlation);
syncStart = peakPos; % Start of PSS symbol (including CP)

%% SSS Detection (Placeholder)
% Assume PSS detected, check SSS one symbol earlier
sssCorr = zeros(3, 1);
for i = 1:3
    sssTest = exp(1i * pi * (0:sssLen-1)' / sssLen * i);
    sssTestGrid = zeros(Nsubcarriers, 1);
    sssTestGrid(6:67) = sssTest;
    sssFreq = zeros(Nfft, 1);
    sssFreq(startIdx:startIdx+Nsubcarriers-1) = sssTestGrid;
    sssTestTime = ifft(sssFreq, Nfft);
    sssTestTimeCP = [sssTestTime(end-Ncp+1:end); sssTestTime];
    sssSegment = rxFrame(syncStart-(Nfft+Ncp):syncStart-(Nfft+Ncp)+length(pssTimeCP)-1);
    if length(sssSegment) == length(sssTestTimeCP)
        sssCorr(i) = abs(sum(conj(sssTestTimeCP) .* sssSegment));
    end
end
[~, detected_N_ID_1] = max(sssCorr);

%% OFDM Demodulation
rxGrid = zeros(Nsubcarriers, numSymbols);
for sym = 1:numSymbols
    startIdx = syncStart + (sym-6) * (Nfft + Ncp); % Adjust for PSS in symbol 6
    if startIdx+Ncp+Nfft-1 <= length(rxFrame)
        rxSymbol = rxFrame(startIdx+Ncp:startIdx+Ncp+Nfft-1);
        rxFFT = fft(rxSymbol, Nfft);
        rxGrid(:, sym) = rxFFT((Nfft-Nsubcarriers)/2+1:(Nfft+Nsubcarriers)/2);
    end
end

%% Channel Estimation (LS)
H_est = zeros(Nsubcarriers, numSymbols);
for sym = crsPos
    H_est(:, sym+1) = rxGrid(:, sym+1) ./ crsSeq; % LS
end
for sc = 1:Nsubcarriers
    H_est(sc, :) = interp1(crsPos+1, H_est(sc, crsPos+1), 1:numSymbols, 'linear', 'extrap');
end

%% PBCH Decoding (Placeholder)
pbchRx = rxGrid(:, 8:11) ./ H_est(:, 8:11);
pbchRx = pbchRx(:);
pbchBits = real(pbchRx) < 0;
disp('PBCH Decoded (Placeholder)');

%% PDSCH Processing (Placeholder)
pdschRx = rxGrid(:, 1:5) ./ H_est(:, 1:5);
pdschRx = pdschRx(:);
pdschBits = real(pdschRx) < 0;
disp('PDSCH Decoded (Placeholder)');

%% Performance Metrics
numTrials = 100;
pssSuccess = 0;
for trial = 1:numTrials
    rxTrial = txFrame + sqrt(noisePower/2) * (randn(size(txFrame)) + 1i * randn(size(txFrame)));
    corr = zeros(correlationLength, 1);
    for i = 1:correlationLength
        segment = rxTrial(i:i+length(pssTimeCP)-1);
        if length(segment) == length(pssTimeCP)
            corr(i) = abs(sum(conj(pssTimeCP) .* segment));
        end
    end
    [~, peak] = max(corr);
    if abs(peak - syncStart) < 5
        pssSuccess = pssSuccess + 1;
    end
end
pssProb = pssSuccess / numTrials;
fprintf('PSS Detection Probability: %.2f%%\n', pssProb * 100);

%% Plotting
figure;
plot(real(txFrame));
title('Transmitted LTE Subframe');
xlabel('Sample');
ylabel('Amplitude');

figure;
plot(correlation);
hold on;
plot(peakPos, peakVal, 'ro');
title('PSS Correlation Detection');
xlabel('Sample Index');
ylabel('Correlation Magnitude');
legend('Correlation', 'Peak');

figure;
scatter(real(rxGrid(:)), imag(rxGrid(:)), 'filled');
title('Received Constellation');
xlabel('In-phase');
ylabel('Quadrature');
grid on;

figure;
surf(1:numSymbols, 1:Nsubcarriers, abs(H_est));
title('Channel Magnitude Estimate');
xlabel('OFDM Symbol');
ylabel('Subcarrier');
zlabel('Magnitude');
