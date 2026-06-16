clear; clc; close all;

%% Echolocation Tunnel Mapping Simulation
% Bandpass filtering + amplifier + convolution/matched filtering
% Adds material-based distance correction algorithm
% Uses chirp signal instead of single-frequency sine pulse
% Sensor moves to the next detected wall point before each pulse

%% Setup
fs = 48000;
T = 0.12;
t = 0:1/fs:T;
c = 343;

numSignals = 10;

%% Sensor movement setup
% Each pulse scans at a different angle
angles = linspace(-90, 90, numSignals);

% Sensor starts at origin
currentSensorX = 0;
currentSensorY = 0;

% Store sensor position before each pulse
sensorX = zeros(numSignals, 1);
sensorY = zeros(numSignals, 1);

% Store estimated wall point after each pulse
estimatedWallX = zeros(numSignals, 1);
estimatedWallY = zeros(numSignals, 1);

%% Transmitted chirp pulse
chirpStartFreq = 3000;
chirpEndFreq = 5000;
pulseDur = 0.004;

Npulse = sum(t <= pulseDur);
tPulse = t(1:Npulse);

% Create chirp from 3 kHz to 5 kHz
txTemplate = chirp(tPulse, chirpStartFreq, pulseDur, chirpEndFreq, 'linear');

% Apply Hann window to smooth pulse edges
txTemplate = txTemplate .* hann(Npulse)';

% Full-length transmitted pulse
txPulse = zeros(size(t));
txPulse(1:Npulse) = txTemplate;

%% Material database
materials = ["Concrete", "Brick", "Metal", "Wood", "Rock", "Wet Rock"];

reflectionCoeff = [0.75, 0.60, 0.95, 0.35, 0.70, 0.50];

materialPhase = [0, pi/6, pi/2, pi/4, pi/3, 2*pi/3];

% Correction factors for estimated distance
% Wet Rock is given a stronger correction because moisture/roughness can distort echo timing
distanceCorrectionFactor = [1.00, 1.01, 0.99, 1.03, 1.02, 1.08];

%% Bandpass filter
lowCut = 3000;
highCut = 5000;

[b, a] = butter(4, [lowCut highCut]/(fs/2), 'bandpass');

%% Amplifier setup
% Amplifier is placed after the bandpass filter.
% This boosts the cleaned echo before matched filtering.
amplifierGain = 8;

% Saturation limit to model real amplifier clipping
ampLimit = 1;

%% Storage variables
rawSignals = zeros(numSignals, length(t));
filteredSignals = zeros(numSignals, length(t));
amplifiedSignals = zeros(numSignals, length(t));
matchedOutputs = cell(numSignals, 1);

actualDistances = zeros(numSignals, 1);
estimatedDistances = zeros(numSignals, 1);
correctedDistances = zeros(numSignals, 1);

actualDelays = zeros(numSignals, 1);
estimatedDelays = zeros(numSignals, 1);

actualPhases = zeros(numSignals, 1);
estimatedPhases = zeros(numSignals, 1);

chosenMaterials = strings(numSignals, 1);
estimatedMaterials = strings(numSignals, 1);

%% Generate and process 10 received signals
for k = 1:numSignals

    %% Store current sensor position before sending pulse
    sensorX(k) = currentSensorX;
    sensorY(k) = currentSensorY;

    % Current scanning angle
    scanAngle = angles(k);

    %% Random distance and material
    distance = 1 + 9*rand;
    delay = 2*distance/c;
    delaySamples = round(delay*fs);
    
    matIndex = randi(length(materials));
    material = materials(matIndex);

    alpha = reflectionCoeff(matIndex);
    phaseShift = materialPhase(matIndex) + 0.25*randn;

    %% Create phase-shifted chirp echo pulse
    analyticChirp = hilbert(txTemplate);
    echoPulse = real(analyticChirp .* exp(1j*phaseShift));

    distanceLoss = 1/(distance^2);
    echoPulse = alpha * distanceLoss * echoPulse;

    %% Place delayed echo
    rx = zeros(size(t));

    if delaySamples + Npulse <= length(t)
        rx(delaySamples + 1 : delaySamples + Npulse) = echoPulse;
    end

    %% Add unwanted noise
    whiteNoise = 0.015*randn(size(t));
    lowNoise = 0.02*sin(2*pi*300*t + 2*pi*rand);
    highNoise = 0.015*sin(2*pi*11000*t + 2*pi*rand);

    rxNoisy = rx + whiteNoise + lowNoise + highNoise;

    %% Bandpass filter
    rxFiltered = filtfilt(b, a, rxNoisy);

    %% Amplifier after bandpass filter
    rxAmplified = amplifierGain * rxFiltered;

    % Clip the amplified signal if it exceeds the amplifier limit
    rxAmplified(rxAmplified > ampLimit) = ampLimit;
    rxAmplified(rxAmplified < -ampLimit) = -ampLimit;

    %% Convolution / matched filter
    matched = conv(rxAmplified, fliplr(txTemplate), 'same');

    [~, peakIndex] = max(abs(matched));

    estimatedDelay = t(peakIndex) - pulseDur/2;

    if estimatedDelay < 0
        estimatedDelay = 0;
    end

    estimatedDistance = c*estimatedDelay/2;

    %% Phase estimation for chirp
    analyticRx = hilbert(rxAmplified);

    rxPhaseAtPeak = angle(analyticRx(peakIndex));

    % Estimate the local transmitted chirp phase at the echo peak
    localTime = mod(t(peakIndex), pulseDur);

    chirpRate = (chirpEndFreq - chirpStartFreq)/pulseDur;

    txPhaseAtPeak = 2*pi*(chirpStartFreq*localTime + ...
        0.5*chirpRate*localTime^2);

    phaseEstimate = angle(exp(1j*(rxPhaseAtPeak - txPhaseAtPeak)));

    %% Material classification
    % Estimate material based on phase and amplified echo strength
    echoStrength = max(abs(rxAmplified));

    estimatedMaterial = classifyMaterial(phaseEstimate, echoStrength);

    %% Correct distance based on estimated material
    correctedDistance = correctDistanceByMaterial(estimatedDistance, estimatedMaterial);

    %% Calculate estimated wall position from current sensor position
    estimatedWallX(k) = currentSensorX + correctedDistance*cosd(scanAngle);
    estimatedWallY(k) = currentSensorY + correctedDistance*sind(scanAngle);

    %% Move sensor to the detected wall point before the next pulse
    currentSensorX = estimatedWallX(k);
    currentSensorY = estimatedWallY(k);

    %% Store data
    rawSignals(k, :) = rxNoisy;
    filteredSignals(k, :) = rxFiltered;
    amplifiedSignals(k, :) = rxAmplified;
    matchedOutputs{k} = matched;

    actualDistances(k) = distance;
    estimatedDistances(k) = estimatedDistance;
    correctedDistances(k) = correctedDistance;

    actualDelays(k) = delay;
    estimatedDelays(k) = estimatedDelay;

    actualPhases(k) = angle(exp(1j*phaseShift));
    estimatedPhases(k) = phaseEstimate;

    chosenMaterials(k) = material;
    estimatedMaterials(k) = estimatedMaterial;
end

%% Results table
results = table((1:numSignals)', chosenMaterials, estimatedMaterials, ...
    actualDistances, estimatedDistances, correctedDistances, ...
    sensorX, sensorY, estimatedWallX, estimatedWallY, ...
    actualPhases, estimatedPhases, ...
    'VariableNames', {'Signal', 'ActualMaterial', 'EstimatedMaterial', ...
    'ActualDistance_m', 'RawEstimatedDistance_m', 'CorrectedDistance_m', ...
    'SensorX_m', 'SensorY_m', 'WallX_m', 'WallY_m', ...
    'ActualPhase_rad', 'EstimatedPhase_rad'});

disp(results);

%% Plot transmitted chirp pulse
figure;
plot(tPulse, txTemplate, 'LineWidth', 1.5);
grid on;

xlabel('Time (s)');
ylabel('Amplitude');
title('Transmitted Chirp Pulse');

%% Plot raw noisy signals
figure;
hold on;
grid on;

for k = 1:numSignals
    plot(t, rawSignals(k, :) + 0.12*k, 'LineWidth', 1);
end

xlabel('Time (s)');
ylabel('Amplitude + Offset');
title('Raw Noisy Received Echo Signals');

%% Plot bandpass filtered signals
figure;
hold on;
grid on;

for k = 1:numSignals
    plot(t, filteredSignals(k, :) + 0.12*k, 'LineWidth', 1);
end

xlabel('Time (s)');
ylabel('Amplitude + Offset');
title('Bandpass Filtered Echo Signals');

%% Plot amplified filtered signal
figure;
hold on;
grid on;

for k = 1:numSignals
    plot(t, amplifiedSignals(k, :) + 0.12*k, 'LineWidth', 1);
end

xlabel('Time (s)');
ylabel('Amplitude + Offset');
title('Amplified Filtered Signal');

%% Plot matched filter outputs
figure;
hold on;
grid on;

for k = 1:numSignals
    matched = matchedOutputs{k};

    if max(abs(matched)) ~= 0
        matchedNorm = matched / max(abs(matched));
    else
        matchedNorm = matched;
    end

    plot(t, matchedNorm + 1.2*k, 'LineWidth', 1);
end

xlabel('Time (s)');
ylabel('Matched Filter Output + Offset');
title('Convolution / Matched Filter Outputs After Amplification');

%% Compare actual, raw estimated, and corrected distance
figure;
plot(1:numSignals, actualDistances, 'o-', 'LineWidth', 1.5);
hold on;
plot(1:numSignals, estimatedDistances, 'x--', 'LineWidth', 1.5);
plot(1:numSignals, correctedDistances, 's-.', 'LineWidth', 1.5);
grid on;

xlabel('Signal Number');
ylabel('Distance (m)');
title('Actual vs Raw Estimated vs Corrected Distance');
legend('Actual Distance', 'Raw Estimated Distance', 'Corrected Distance');

%% 2D Tunnel Mapping with Moving Sensor
figure;
hold on;
grid on;
axis equal;

% Plot estimated wall points
plot(estimatedWallX, estimatedWallY, 'o-', 'LineWidth', 2, 'MarkerSize', 8);

% Plot sensor positions before each pulse
plot(sensorX, sensorY, 'kx--', 'LineWidth', 1.5, 'MarkerSize', 8);

% Draw rays from sensor position to detected wall point
for k = 1:numSignals
    plot([sensorX(k), estimatedWallX(k)], ...
         [sensorY(k), estimatedWallY(k)], ':', 'LineWidth', 1);

    text(sensorX(k), sensorY(k), "  S" + num2str(k));
    text(estimatedWallX(k), estimatedWallY(k), "  Wall " + num2str(k));
end

xlabel('X Position (m)');
ylabel('Y Position (m)');
title('2D Tunnel Map with Moving Sensor');

legend('Estimated Wall Path', 'Sensor Position Before Each Pulse', 'Echo Ray');

%% 2D Tunnel Map with Moving Sensor and Material Labels
figure;
hold on;
grid on;
axis equal;

% Plot wall points and material labels
for k = 1:numSignals
    plot(estimatedWallX(k), estimatedWallY(k), 'o', 'MarkerSize', 10, 'LineWidth', 2);

    text(estimatedWallX(k), estimatedWallY(k), ...
        "  " + estimatedMaterials(k) + newline + ...
        "  Wall " + num2str(k));

    % Draw pulse ray from moving sensor to wall
    plot([sensorX(k), estimatedWallX(k)], ...
         [sensorY(k), estimatedWallY(k)], ':', 'LineWidth', 1);
end

% Connect wall points
plot(estimatedWallX, estimatedWallY, 'k--', 'LineWidth', 1.5);

% Plot moving sensor positions
plot(sensorX, sensorY, 'kx--', 'MarkerSize', 10, 'LineWidth', 2);

for k = 1:numSignals
    text(sensorX(k), sensorY(k), "  Sensor " + num2str(k));
end

xlabel('X Position (m)');
ylabel('Y Position (m)');
title('2D Tunnel Map with Moving Sensor and Estimated Surface Materials');

legend('Estimated Wall Point', 'Echo Ray', 'Estimated Wall Path', 'Sensor Positions');

%% Local functions

function estimatedMaterial = classifyMaterial(phaseEstimate, echoStrength)

    phaseEstimate = abs(phaseEstimate);

    if echoStrength > 0.05 && phaseEstimate > 1.1
        estimatedMaterial = "Metal";

    elseif echoStrength < 0.025 && phaseEstimate > 1.5
        estimatedMaterial = "Wet Rock";

    elseif echoStrength < 0.025 && phaseEstimate < 1.0
        estimatedMaterial = "Wood";

    elseif phaseEstimate > 0.8 && phaseEstimate <= 1.5
        estimatedMaterial = "Rock";

    elseif phaseEstimate > 0.3 && phaseEstimate <= 0.8
        estimatedMaterial = "Brick";

    else
        estimatedMaterial = "Concrete";
    end
end

function correctedDistance = correctDistanceByMaterial(estimatedDistance, estimatedMaterial)

    switch estimatedMaterial

        case "Concrete"
            correctionFactor = 1.00;

        case "Brick"
            correctionFactor = 1.01;

        case "Metal"
            correctionFactor = 0.99;

        case "Wood"
            correctionFactor = 1.03;

        case "Rock"
            correctionFactor = 1.02;

        case "Wet Rock"
            correctionFactor = 1.08;

        otherwise
            correctionFactor = 1.00;
    end

    correctedDistance = estimatedDistance * correctionFactor;
end