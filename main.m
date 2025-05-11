%% Parameters
% LTE Parameters for 3 MHz Bandwidth
bandwidth = 3e6;            % Bandwidth in Hz
Nrb = 15;                   % Number of Resource Blocks for 3 MHz
Nsc = 12;                   % Subcarriers per Resource Block
Nsubcarriers = Nrb * Nsc;   % Total subcarriers = 180
Nfft = 256;                 % FFT size (next power of 2 ≥ Nsubcarriers)
Ncp = 144;                  % Normal cyclic prefix length (samples)
numSymbols = 14;            % OFDM symbols per subframe (Normal CP, 2 slots)
Fs = 4.8e6;                 % Sampling rate (typically 1.6 × BW for 3 MHz)
subframeDuration = 1e-3;    % Subframe duration = 1 ms
SNR_dB = 10;                % Signal-to-Noise Ratio in dB
useFading = true;           % Enable multipath fading

% Cell ID Parameters
N_ID_2 = 0;                 % PSS index (0, 1, 2)
N_ID_1 = 0;                 % SSS index (0–167)
PCI = 3 * N_ID_1 + N_ID_2;  % Physical Cell Identity (0–503)


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
rxFrame = rxFrame .* exp(-1i * 2 * pi * cfoEst * t); % CFO correction

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
    startIdx = syncStart + (sym-6) * (Nfft + Ncp);  % Align symbol 6 with syncStart
    if startIdx > 0 && (startIdx + Ncp + Nfft - 1) <= length(rxFrame)
        rxSymbol = rxFrame(startIdx + Ncp : startIdx + Ncp + Nfft - 1);
        rxFFT = fft(rxSymbol, Nfft);
        rxGrid(:, sym) = rxFFT((Nfft - Nsubcarriers)/2 + 1 : (Nfft + Nsubcarriers)/2);
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
%% Improved Channel Estimation (Linear Interpolation)
H_est_refined = zeros(Nsubcarriers, numSymbols);

% Initial estimation using LS
for sym = crsPos
    H_est(:, sym+1) = rxGrid(:, sym+1) ./ crsSeq; % Least Squares estimation
end

% Refine the channel estimate by interpolating between CRS symbols
for sc = 1:Nsubcarriers
    % Linear interpolation for each subcarrier across the symbols
    H_est_refined(sc, :) = interp1(crsPos+1, H_est(sc, crsPos+1), 1:numSymbols, 'linear', 'extrap');
end

% Optional: Use spline interpolation for smoother results (more accurate than linear)
% Known pilot positions (symbols with CRS)
x = crsPos + 1;       % Pilot symbol indices (row vector)
xt = 1:numSymbols;    % All OFDM symbol indices (row vector)
rho = 0.9;            % Temporal correlation factor
noiseVar = 10^(-SNR_dB / 10);  % Noise variance

% Preallocate
H_est_linear = zeros(Nsubcarriers, numSymbols);
H_est_mmse = zeros(Nsubcarriers, numSymbols);

% Compute correlation matrices for MMSE
R_xx = rho .^ abs(x' - x);           % [4 x 4]
R_xt = rho .^ abs(xt(:) - x);        % [14 x 4]

for sc = 1:Nsubcarriers
    y = H_est(sc, x).';  % pilot-based estimate at subcarrier sc [4 x 1]
    
    % Linear Interpolation
    H_est_linear(sc, :) = interp1(x, H_est(sc, x), xt, 'linear', 'extrap');
    
    % MMSE Interpolation
    H_est_mmse(sc, :) = (R_xt / (R_xx + noiseVar * eye(length(x)))) * y;
end

%%

% Plot the refined channel estimates
figure;
subplot(2,1,1);
surf(1:numSymbols, 1:Nsubcarriers, abs(H_est_linear));
title('Refined Channel Estimate (Linear Interpolation)');
xlabel('OFDM Symbol');
ylabel('Subcarrier');
zlabel('Magnitude');
colorbar;

subplot(2,1,2);
surf(1:numSymbols, 1:Nsubcarriers, abs(H_est_mmse));
title('Refined Channel Estimate (MMSE Interpolation)');
xlabel('OFDM Symbol');
ylabel('Subcarrier');
zlabel('Magnitude');
colorbar;

%% PCFICH Detection and Decoding (QPSK on fixed locations)

% PCFICH uses 4 REGs (16 REs) in symbol 0
% Simplified: Use first 16 subcarriers in symbol 0
pcfichREs = 1:16;

% Extract and equalize
pcfichSymbols = rxGrid(pcfichREs, 1) ./ H_est(pcfichREs, 1);

% QPSK Demodulation
pcfichBits = zeros(32, 1);
for k = 1:16
    pcfichBits(2*k-1) = real(pcfichSymbols(k)) < 0;
    pcfichBits(2*k)   = imag(pcfichSymbols(k)) < 0;
end

% Interpret first 2 bits as control format (0 to 3)
controlFormat = bi2de(pcfichBits(1:2)', 'left-msb') + 1;
fprintf('PCFICH Decoded Control Format Indicator (CFI): %d OFDM symbols\n', controlFormat);

%% PDCCH Detection and Decoding (QPSK on symbols 1, 2, 3)
pdcchSymbols = rxGrid(1:12, 2:4) ./ H_est(1:12, 2:4);

pdcchBits = zeros(144, 1); 
for sym = 1:3
    for subcarrier = 1:12
        idx = (sym-1)*12 + subcarrier;
        pdcchBits(2*idx-1) = real(pdcchSymbols(subcarrier, sym)) < 0;
        pdcchBits(2*idx)   = imag(pdcchSymbols(subcarrier, sym)) < 0;
    end
end

disp('PDCCH Decoded (Placeholder)');
disp(pdcchBits(1:16)); 

%% PBCH Decoding and MIB Extraction 

pbchRx = rxGrid(:, 8:11) ./ H_est(:, 8:11);
pbchRx = pbchRx(:);
pbchBits = real(pbchRx) < 0;

% Assume first 24 bits are MIB (bypass BCH decoding)
if length(pbchBits) >= 24
    mibBits = pbchBits(1:24);
    systemFrameNumber = bi2de(mibBits(1:8), 'left-msb');
    phichDuration = mibBits(9);
    phichResource = bi2de(mibBits(10:11), 'left-msb');
    bandwidthConfig = bi2de(mibBits(12:14), 'left-msb');

    fprintf('MIB Decoded (No BCH Decoding):\n');
    fprintf(' - SFN: %d\n', systemFrameNumber);
    fprintf(' - PHICH Duration: %d\n', phichDuration);
    fprintf(' - PHICH Resource: %d\n', phichResource);
    fprintf(' - Bandwidth Config: %d\n', bandwidthConfig);
else
    warning('PBCH decoding failed: insufficient bits.');
end


%% PDSCH Decoding (QPSK on symbols 1-5)
pdschSymbols = rxGrid(1:12, 1:5) ./ H_est(1:12, 1:5);

pdschBits = zeros(120, 1); 
for sym = 1:5
    for subcarrier = 1:12
        idx = (sym-1)*12 + subcarrier;
        pdschBits(2*idx-1) = real(pdschSymbols(subcarrier, sym)) < 0;
        pdschBits(2*idx)   = imag(pdschSymbols(subcarrier, sym)) < 0;
    end
end

disp('PDSCH Decoded (Placeholder)');
disp(pdschBits(1:16)); 


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

%% PBCH CRC Checking
% CRC Polynomial (CRC-24C)
crcPoly = [1 0 0 1 1 1 0 1 1 1];  % CRC-24C (standard LTE polynomial)

% Check CRC of MIB bits (first 24 bits)
mibCRCValid = checkCRC(mibBits, crcPoly);

if mibCRCValid
    fprintf('PBCH CRC Check Passed: MIB is valid.\n');
else
    fprintf('PBCH CRC Check Failed: MIB is invalid.\n');
end

%% plot
% Parameters for simulation
snrRange = -10:2:30; % Range of SNR values (in dB)
numFrames = 100; % Number of frames to simulate
numSymbolsPerFrame = 100; % Number of symbols per frame
subframeDuration = 1e-3; % Duration of a subframe in seconds
totalBlocks = zeros(1, length(snrRange));
pssDetectionCount = zeros(1, length(snrRange));
sssDetectionCount = zeros(1, length(snrRange));
cellIDAccuracy = zeros(1, length(snrRange));
timingError = zeros(1, length(snrRange));
mibSuccessCount = zeros(1, length(snrRange));
pdcchSuccessCount = zeros(1, length(snrRange));
pdschSuccessCount = zeros(1, length(snrRange));
throughputBits = zeros(1, length(snrRange));

% Loop over SNR range
for snrIdx = 1:length(snrRange)
    snr = snrRange(snrIdx);
    
    % Simulate frames
    for frameIdx = 1:numFrames
        % Simulate transmission (e.g., PSS/SSS detection, cell ID detection, etc.)
        % Example placeholders for actual simulation steps
        % PSS/SSS Detection
        pssDetected = false;
        sssDetected = false;
        if rand() > 0.5 % Simulating a random detection condition
            pssDetected = true;
        end
        if rand() > 0.5 % Simulating a random detection condition
            sssDetected = true;
        end

        pssDetectionCount(snrIdx) = pssDetectionCount(snrIdx) + pssDetected;
        sssDetectionCount(snrIdx) = sssDetectionCount(snrIdx) + sssDetected;

        % Cell ID detection accuracy
        txCellID = randi([0, 1000]); % Example cell ID
        rxCellID = txCellID; % Assume detection is accurate for this example
        correctCellIDDetection = (rxCellID == txCellID);
        cellIDAccuracy(snrIdx) = cellIDAccuracy(snrIdx) + correctCellIDDetection;

        % Timing synchronization error
        estimatedOffset = randn() * 10; % Random offset
        trueOffset = 0; % Assume perfect timing for simplicity
        syncError = abs(estimatedOffset - trueOffset);
        timingError(snrIdx) = timingError(snrIdx) + syncError;

        % MIB and PBCH decoding success rate
        mibCRC = randi([0, 1]); % Simulating MIB CRC result (0 = success, 1 = failure)
        mibSuccess = mibCRC == 0;
        mibSuccessCount(snrIdx) = mibSuccessCount(snrIdx) + mibSuccess;

        % PDCCH blind decoding success probability
        dci = rand() > 0.5; % Simulating a random DCI detection
        pdcchSuccess = ~isempty(dci);
        pdcchSuccessCount(snrIdx) = pdcchSuccessCount(snrIdx) + pdcchSuccess;

        % Block Error Rate (BLER) for PDSCH
        transportBlockCRC = randi([0, 1]); % Simulating CRC check result (0 = success, 1 = failure)
        if transportBlockCRC == 0
            pdschSuccessCount(snrIdx) = pdschSuccessCount(snrIdx) + 1;
        end
        totalBlocks(snrIdx) = totalBlocks(snrIdx) + 1;

        % Throughput measurements (in bits per second)
        decodedBits = randi([0, 1], [1, 1000]); % Random decoded bits
        if transportBlockCRC == 0
            throughputBits(snrIdx) = throughputBits(snrIdx) + numel(decodedBits);
        end
    end
end

% Normalize all metrics
pssDetectionProb = pssDetectionCount / numFrames;
sssDetectionProb = sssDetectionCount / numFrames;
cellIDAccuracy = cellIDAccuracy / numFrames;
timingError = timingError / numFrames;
mibSuccessRate = mibSuccessCount / numFrames;
pdcchSuccessRate = pdcchSuccessCount / numFrames;
bler = 1 - (pdschSuccessCount ./ totalBlocks);
throughput = throughputBits / (numFrames * subframeDuration);  % bits/sec

% Plotting Results

% PSS/SSS Detection vs SNR
figure;
plot(snrRange, pssDetectionProb, '-o'); title('PSS Detection vs SNR'); xlabel('SNR (dB)'); ylabel('Probability');

% SSS Detection vs SNR
figure;
plot(snrRange, sssDetectionProb, '-o'); title('SSS Detection vs SNR'); xlabel('SNR (dB)'); ylabel('Probability');

% Cell ID Accuracy vs SNR
figure;
plot(snrRange, cellIDAccuracy, '-o'); title('Cell ID Accuracy vs SNR'); xlabel('SNR (dB)'); ylabel('Accuracy');

% Timing Synchronization Error vs SNR
figure;
plot(snrRange, timingError, '-o'); title('Timing Error vs SNR'); xlabel('SNR (dB)'); ylabel('Synchronization Error (ms)');

% MIB Success Rate vs SNR
figure;
plot(snrRange, mibSuccessRate, '-o'); title('MIB Success Rate vs SNR'); xlabel('SNR (dB)'); ylabel('Success Rate');

% PDCCH Success Rate vs SNR
figure;
plot(snrRange, pdcchSuccessRate, '-o'); title('PDCCH Success Rate vs SNR'); xlabel('SNR (dB)'); ylabel('Success Rate');

% Block Error Rate (BLER) for PDSCH vs SNR
figure;
plot(snrRange, bler, '-o'); title('BLER for PDSCH vs SNR'); xlabel('SNR (dB)'); ylabel('Block Error Rate');

% Throughput vs SNR
figure;
plot(snrRange, throughput, '-o'); title('Throughput vs SNR'); xlabel('SNR (dB)'); ylabel('Throughput (bits/sec)');

function isValid = checkCRC(bits, crcPoly)
    % Ensure bits are row vectors
    bits = bits(:).';  % Ensure bits is a row vector
    
    % Perform CRC checking using the specified polynomial
    numBits = length(bits);
    crc = [bits, zeros(1, length(crcPoly) - 1)]; % Append zeros to bits for CRC calculation
    crcLength = length(crc);
    
    for i = 1:numBits  % Perform the division process
        if crc(i) == 1
            crc(i:i+length(crcPoly)-1) = mod(crc(i:i+length(crcPoly)-1) + crcPoly, 2);
        end
    end
    
    % The CRC is valid if the remainder is zero
    isValid = all(crc(end-length(crcPoly)+2:end) == 0);  % Check remainder for validity
end