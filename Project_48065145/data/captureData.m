% Full example: capturing, processing, decoding, and displaying ADS-B data
centerFreq = 1090e6;
sampleRate = 2.4e6;
samplesPerFrame = 9.6e4;
tunerGain = 25;
% Define the RTL-SDR receiver
sdrReceiver = comm.SDRRTLReceiver('CenterFrequency', centerFreq, ...
    'SampleRate', sampleRate, ...
    'EnableTunerAGC', false, ...
    'TunerGain', tunerGain, ... 
    'OutputDataType', 'single', ... 
    'SamplesPerFrame', samplesPerFrame, ...
    'FrequencyCorrection',0);

% Capture samples
[samples, len, overrun] = step(sdrReceiver);

% Release the receiver
release(sdrReceiver);

% Plot the captured samples
figure;
plot(real(samples));
title('Captured ADS-B Signal');
xlabel('Sample Index');
ylabel('Amplitude');

% Process the captured samples to detect ADS-B messages
threshold = 0.01; % Set a threshold for PPM demodulation
binarySignal = demodulateADSB(samples, threshold);

% Extract ADS-B messages
adsbMessages = extractADSBMessages(binarySignal);

% Decode ADS-B messages
decodedData = decodeADSBMessages(adsbMessages);

% Convert to table if data exists
if ~isempty(decodedData)
    adsbTable = cell2table(decodedData, 'VariableNames', {'Callsign'});
    % Display table
    disp(adsbTable);
else
    disp('No valid ADS-B messages decoded.');
end

% Functions for demodulating, extracting, and decoding ADS-B signals
function binarySignal = demodulateADSB(samples, threshold)
    % Demodulate ADS-B signal to extract binary signal
    binarySignal = samples > threshold;
end

function adsbMessages = extractADSBMessages(binarySignal)
    adsbMessages = {};
    preamble = [1 0 0 0 0 0 1 0 1 0 0 0 0 0 0]; % Known preamble pattern
    preambleLen = length(preamble);
    messageLen = 112; % ADS-B message length in bits
    i = 1;
    
    while i <= length(binarySignal) - (preambleLen + messageLen)
        if isequal(binarySignal(i:i+preambleLen-1), preamble)
            % Extract the message bits following the preamble
            messageBits = binarySignal(i + preambleLen:i + preambleLen + messageLen - 1);
            adsbMessages{end + 1} = messageBits; %#ok<AGROW>
            i = i + preambleLen + messageLen; % Move to the end of the message
        else
            i = i + 1; % Move to the next bit
        end
    end
end

function decodedData = decodeADSBMessages(adsbMessages)
    decodedData = {};
    for k = 1:length(adsbMessages)
        msg = adsbMessages{k};
        % Decode message fields (example for callsign messages)
        df = bi2de(msg(1:5)); % Downlink format
        if df == 17 % DF17 is used for ADS-B messages
            tc = bi2de(msg(33:37)); % Type code
            if tc >= 1 && tc <= 4 % Aircraft identification
                callsign = char(bi2de(reshape(msg(41:96), 6, 8)') + '0'); % Convert bits to characters
                decodedData{end + 1, 1} = callsign; %#ok<AGROW>
                % Add additional fields like position, altitude, etc.
            end
        end
    end
end
