% Define eNodeB configuration (basic)
enb = struct;
enb.NDLRB = 50;                % Number of Downlink Resource Blocks (10 MHz => 50 RBs)
enb.CellRefP = 1;              % Number of cell-specific antenna ports
enb.NCellID = 0;               % Physical cell ID
enb.CyclicPrefix = 'Normal';   % Normal CP
enb.DuplexMode = 'FDD';        % Frame type 1 (FDD)
enb.NSubframe = 0;             % Starting subframe number
enb.TotSubframes = 10;         % Total 10 subframes = 1 frame

grid = lteDLResourceGrid(enb, enb.CellRefP);    % Empty resource grid

% Fill the grid with dummy (random) data
% For example, use random QPSK-modulated data
% data = lteDLDMRS(enb); % Generate Demodulation Reference Signal (optional)

% Alternatively, generate random data for all resource elements
grid(:) = (randi([0 1], size(grid)) * 2 - 1) + 1i * (randi([0 1], size(grid)) * 2 - 1);

% Perform OFDM modulation to get the time-domain waveform
[txWaveform, txInfo] = lteOFDMModulate(enb, grid);

% Visualization
figure;
imagesc(abs(grid(:,:,1))); % Plot first antenna port
xlabel('OFDM Symbols');
ylabel('Subcarriers');
title('LTE Resource Grid (FDD - 1 Antenna)');
colorbar;