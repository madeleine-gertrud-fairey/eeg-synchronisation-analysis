clc; clearvars; close all;

%% ================== USER SETTINGS ==================
edfPath   = "C:\Users\Usuario\Desktop\TFG\TrainingP1Data\P1\P1_56.edf";
chanNames = ["F3","T3"];
chanIdx   = [1 2];

% Display window for raw EEG (seconds)
winSecs   = 20;            % bottom EEG window length
stepSecs  = 5;             % how far each arrow press moves

% Chunk length for recomputing MPC (seconds)
chunkSecs = 600;           % 10 minutes

% MPC parameters (as requested)
mpcWin_s  = 20;            % 20 s sliding window
mpcStep_s = 2;             % 2 s step size

% MPC bands shown in the top two panels
band1 = [8 12];            % alpha
band2 = [4 8];             % theta

filt_order = 4;

% Seizure markers (optional)
markTimes = [ ...
    datetime('12/11/2020 00:42:14','InputFormat','dd/MM/yyyy HH:mm:ss'), ...
    datetime('12/11/2020 00:49:38','InputFormat','dd/MM/yyyy HH:mm:ss') ...
];

%% ================== READ EDF + BUILD CONTINUOUS VECTORS ==================
[edfStartDT, fs] = inferEdfStartAndFsRaw(edfPath, chanIdx(1));
dt = 1/fs;

tbl = edfread(edfPath);
timeOffsets = seconds(tbl.Properties.RowTimes);

nRec = height(tbl);
sig1_all = cell(nRec,1);
sig2_all = cell(nRec,1);
t_all    = cell(nRec,1);

for i = 1:nRec
    x1 = tbl{i,chanIdx(1)}{1};
    x2 = tbl{i,chanIdx(2)}{1};

    base_time = edfStartDT + seconds(timeOffsets(i));
    local_t = base_time + seconds((0:numel(x1)-1)*dt);

    sig1_all{i} = x1(:);
    sig2_all{i} = x2(:);
    t_all{i}    = local_t(:);
end

X1 = vertcat(sig1_all{:});
X2 = vertcat(sig2_all{:});
t  = vertcat(t_all{:});

totalSamples = numel(X1);
winSamples   = round(winSecs*fs);
stepSamps    = round(stepSecs*fs);
chunkSamples = round(chunkSecs*fs);

%% ================== INITIAL CHUNK MPC COMPUTE ==================
currentChunkStart = 1;
initLen = min(chunkSamples, totalSamples);

[m1_t, m1_val] = computeMPC_forChunk(X1(1:initLen), X2(1:initLen), ...
    t(1), fs, band1, filt_order, mpcWin_s, mpcStep_s);

[m2_t, m2_val] = computeMPC_forChunk(X1(1:initLen), X2(1:initLen), ...
    t(1), fs, band2, filt_order, mpcWin_s, mpcStep_s);

%% ================== FIGURE LAYOUT ==================
f = figure('Name','EEG MPC Viewer','Units','normalized', ...
    'Position',[0.05 0.05 0.9 0.9], ...
    'KeyPressFcn',@keyControl);

uimenu(f, 'Label', 'Jump to Time', 'Callback', @jumpToTime);

tiledlayout(3,1,'TileSpacing','compact','Padding','tight');

% ---- Row 1: MPC band1 ----
ax1 = nexttile(1);
hMPC1 = plot(ax1, m1_t, m1_val, 'LineWidth', 1.2);
grid(ax1,'on'); ylim(ax1,[0 1.02]);
title(ax1, sprintf('MPC %g-%g Hz (%s–%s)', band1(1), band1(2), chanNames(1), chanNames(2)));
ylabel(ax1,'MPC');
xlim(ax1, [m1_t(1) m1_t(end)]);

% Dark band = EEG view window
hPatchEEG1 = patch(ax1, [t(1) t(1) t(1) t(1)], [0 0 1.02 1.02], ...
    'k','FaceAlpha',0.10,'EdgeColor','none');

% Light band = MPC 20s window (PROSPECTIVE: ends at "now")
hPatchMPC1 = patch(ax1, [t(1) t(1) t(1) t(1)], [0 0 1.02 1.02], ...
    'k','FaceAlpha',0.04,'EdgeColor','none');

% ---- Row 2: MPC band2 ----
ax2 = nexttile(2);
hMPC2 = plot(ax2, m2_t, m2_val, 'LineWidth', 1.2);
grid(ax2,'on'); ylim(ax2,[0 1.02]);
title(ax2, sprintf('MPC %g-%g Hz (%s–%s)', band2(1), band2(2), chanNames(1), chanNames(2)));
ylabel(ax2,'MPC');
xlim(ax2, [m2_t(1) m2_t(end)]);

hPatchEEG2 = patch(ax2, [t(1) t(1) t(1) t(1)], [0 0 1.02 1.02], ...
    'k','FaceAlpha',0.10,'EdgeColor','none');

% Light band = MPC 20s window (PROSPECTIVE: ends at "now")
hPatchMPC2 = patch(ax2, [t(1) t(1) t(1) t(1)], [0 0 1.02 1.02], ...
    'k','FaceAlpha',0.04,'EdgeColor','none');

% ---- Row 3: Raw signals ----
ax3 = nexttile(3);
idx0 = 1:min(winSamples,totalSamples);
h1 = plot(ax3, t(idx0), X1(idx0), 'LineWidth', 1);
hold(ax3,'on');
offset = -300;
h2 = plot(ax3, t(idx0), X2(idx0)+offset, 'LineWidth', 1);
hold(ax3,'off');
ylim(ax3,[-400 100]);
ax3.YTickLabel = [];
xlabel(ax3,'Time'); title(ax3,'EEG Signals');
xlim(ax3, [t(idx0(1)) t(idx0(end))]);

% ---- Initialize the shaded bars to match the first EEG window ----
tWin0 = t(idx0);
xEEG0 = [tWin0(1) tWin0(end) tWin0(end) tWin0(1)];

% PROSPECTIVE MPC window: past 20s ending at EEG window right edge
tMPC0_R = tWin0(end);
tMPC0_L = tMPC0_R - seconds(mpcWin_s);
xMPC0   = [tMPC0_L tMPC0_R tMPC0_R tMPC0_L];

set(hPatchEEG1, 'XData', xEEG0);
set(hPatchEEG2, 'XData', xEEG0);
set(hPatchMPC1, 'XData', xMPC0);
set(hPatchMPC2, 'XData', xMPC0);

% Keep patches behind the MPC line
uistack(hPatchMPC1,'bottom'); uistack(hPatchEEG1,'bottom');
uistack(hPatchMPC2,'bottom'); uistack(hPatchEEG2,'bottom');

% ---- Seizure markers on all axes ----
allAxes = [ax1, ax2, ax3];
for ax = allAxes
    for mt = markTimes
        xline(ax, mt, '--', 'LineWidth', 1.3, 'Color', [0 0 0], ...
            'Label', datestr(mt,'dd-mmm-yyyy HH:MM:SS'), ...
            'LabelOrientation','horizontal', ...
            'LabelVerticalAlignment','bottom', ...
            'LabelHorizontalAlignment','center');
    end
end

%% ---- Store data for callbacks ----
data = struct();
data.fs = fs;
data.winSamples = winSamples;
data.stepSamps  = stepSamps;
data.chunkSamples = chunkSamples;
data.totalSamples = totalSamples;

data.t  = t;
data.X1 = X1;
data.X2 = X2;
data.offset = offset;

data.currentStart = 1;
data.currentChunkStart = currentChunkStart;

data.band1 = band1;
data.band2 = band2;
data.filt_order = filt_order;
data.mpcWin_s = mpcWin_s;
data.mpcStep_s = mpcStep_s;

data.ax1 = ax1; data.ax2 = ax2; data.ax3 = ax3;
data.hMPC1 = hMPC1; data.hMPC2 = hMPC2;
data.hPatchEEG1 = hPatchEEG1; data.hPatchEEG2 = hPatchEEG2;
data.hPatchMPC1 = hPatchMPC1; data.hPatchMPC2 = hPatchMPC2;
data.h1 = h1; data.h2 = h2;

guidata(f,data);

%% ================== CALLBACKS ==================
function keyControl(src,event)
    data = guidata(src);

    switch event.Key
        case 'rightarrow'
            newStart = data.currentStart + data.stepSamps;

        case 'leftarrow'
            newStart = data.currentStart - data.stepSamps;

        case 'j'
            jumpToTime(src);        % full datetime
            return;

        case 'g'
            jumpToClockTime(src);   % HH:MM:SS only
            return;

        otherwise
            return;
    end

    newStart = max(1, min(newStart, data.totalSamples - data.winSamples + 1));

    % Determine chunk for MPC
    chunkStart = floor((newStart-1)/data.chunkSamples)*data.chunkSamples + 1;

    if chunkStart ~= data.currentChunkStart
        chunkEnd = min(chunkStart + data.chunkSamples - 1, data.totalSamples);
        x1c = data.X1(chunkStart:chunkEnd);
        x2c = data.X2(chunkStart:chunkEnd);
        t0  = data.t(chunkStart);

        [tM1, vM1] = computeMPC_forChunk(x1c, x2c, t0, data.fs, ...
            data.band1, data.filt_order, data.mpcWin_s, data.mpcStep_s);

        [tM2, vM2] = computeMPC_forChunk(x1c, x2c, t0, data.fs, ...
            data.band2, data.filt_order, data.mpcWin_s, data.mpcStep_s);

        set(data.hMPC1, 'XData', tM1, 'YData', vM1);
        set(data.hMPC2, 'XData', tM2, 'YData', vM2);

        xlim(data.ax1, [tM1(1) tM1(end)]);
        xlim(data.ax2, [tM2(1) tM2(end)]);

        data.currentChunkStart = chunkStart;
    end

    % Update raw EEG window
    data.currentStart = newStart;
    idx = newStart : newStart + data.winSamples - 1;
    tWin = data.t(idx);

    set(data.h1, 'XData', tWin, 'YData', data.X1(idx));
    set(data.h2, 'XData', tWin, 'YData', data.X2(idx) + data.offset);
    xlim(data.ax3, [tWin(1) tWin(end)]);

    % Update EEG-view dark bar
    xEEG = [tWin(1) tWin(end) tWin(end) tWin(1)];
    set(data.hPatchEEG1, 'XData', xEEG);
    set(data.hPatchEEG2, 'XData', xEEG);

    % PROSPECTIVE MPC window: past 20s ending at EEG window right edge
    tR = tWin(end);
    tL = tR - seconds(data.mpcWin_s);
    xMPC = [tL tR tR tL];
    set(data.hPatchMPC1, 'XData', xMPC);
    set(data.hPatchMPC2, 'XData', xMPC);

    guidata(src,data);
    drawnow;
end

function jumpToTime(src, ~)
    data = guidata(src);

    prompt = {'Enter date and time (dd/MM/yyyy HH:mm:ss):'};
    dlgtitle = 'Jump to Timestamp';
    definput = {datestr(data.t(data.currentStart), 'dd/mm/yyyy HH:MM:SS')};

    answer = inputdlg(prompt, dlgtitle, [1 50], definput);
    if isempty(answer), return; end

    try
        targetTime = datetime(answer{1}, 'InputFormat', 'dd/MM/yyyy HH:mm:ss');
    catch
        errordlg('Invalid datetime format. Use dd/MM/yyyy HH:mm:ss.','Error');
        return;
    end

    [~, idxNearest] = min(abs(data.t - targetTime));
    data.currentStart = max(1, min(idxNearest, data.totalSamples - data.winSamples + 1));
    guidata(src,data);

    % Force redraw using the same update logic
    fakeEvent.Key = 'rightarrow';
    keyControl(src, fakeEvent);
end

function jumpToClockTime(src, ~)
    data = guidata(src);

    prompt = {'Enter clock time (HH:MM:SS):'};
    dlgtitle = 'Go to Clock Time';
    definput = {datestr(data.t(data.currentStart), 'HH:MM:SS')};

    answer = inputdlg(prompt, dlgtitle, [1 40], definput);
    if isempty(answer), return; end

    % Parse time-of-day only
    try
        tod = datetime(answer{1}, 'InputFormat', 'HH:mm:ss');
    catch
        errordlg('Invalid format. Use HH:MM:SS (e.g., 00:42:14).','Error');
        return;
    end
    targetTOD = timeofday(tod);

    % Find next occurrence at/after current time
    curT = data.t(data.currentStart);
    dtCandidate = dateshift(curT, 'start', 'day') + targetTOD;
    if dtCandidate < curT
        dtCandidate = dtCandidate + days(1);
    end

    % Clamp into recording range
    if dtCandidate < data.t(1),   dtCandidate = data.t(1);   end
    if dtCandidate > data.t(end), dtCandidate = data.t(end); end

    [~, idxNearest] = min(abs(data.t - dtCandidate));
    data.currentStart = max(1, min(idxNearest, data.totalSamples - data.winSamples + 1));
    guidata(src,data);

    fakeEvent.Key = 'rightarrow';
    keyControl(src, fakeEvent);
end

%% ================== FUNCTIONS ==================
function [tM, mpcVals] = computeMPC_forChunk(x1, x2, t0, fs, band, order, win_s, step_s)
    f_lo = band(1); f_hi = band(2);

    y1 = bandpass_filter(x1, fs, f_lo, f_hi, order);
    y2 = bandpass_filter(x2, fs, f_lo, f_hi, order);

    p1 = angle(hilbert(y1));
    p2 = angle(hilbert(y2));
    dphi = angle(exp(1j*(p1 - p2)));

    winSamp  = max(2, round(win_s*fs));
    stepSamp = max(1, round(step_s*fs));
    n = numel(dphi);

    starts = 1:stepSamp:(n-winSamp+1);
    mpcVals = zeros(numel(starts),1);
    tM = NaT(numel(starts),1);

    for ii = 1:numel(starts)
        s = starts(ii);
        e = s + winSamp - 1;

        mpcVals(ii) = abs(mean(exp(1j*dphi(s:e))));

        % ---- LEVEL A PROSPECTIVE TIMESTAMP (end of window) ----
        tM(ii) = t0 + seconds((e-1)/fs);
    end
end

function y = bandpass_filter(x, fs, f_lo, f_hi, order)
    nyq = 0.5*fs;
    [b,a] = butter(order, [f_lo/nyq, f_hi/nyq], "bandpass");
    y = filtfilt(b,a,double(x));  % NOTE: still acausal (Level B), but OK for Level A alignment
end

function [startDT, fs] = inferEdfStartAndFsRaw(edfPath, chanIndexForFs)
    fid = fopen(edfPath, 'r', 'ieee-le');
    if fid < 0, error("Cannot open EDF: %s", edfPath); end
    cleanup = onCleanup(@() fclose(fid));

    fixed = fread(fid, 256, '*char')';

    startDateStr = strtrim(string(fixed(169:176)));      % dd.mm.yy
    startTimeStr = replace(strtrim(string(fixed(177:184))), ".", ":"); % hh.mm.ss -> hh:mm:ss

    d = datetime(startDateStr, "InputFormat","dd.MM.yy", "Locale","en_US");
    tt = datetime(startTimeStr, "InputFormat","HH:mm:ss", "Locale","en_US");
    startDT = dateshift(d,"start","day") + timeofday(tt);

    durRec = str2double(strtrim(string(fixed(245:252))));
    ns     = str2double(strtrim(string(fixed(253:256))));

    skipBytes = ns*(16 + 80 + 8 + 8 + 8 + 8 + 8 + 80);
    fseek(fid, 256 + skipBytes, 'bof');
    nsampChars = fread(fid, ns*8, '*char')';

    nsamp = zeros(ns,1);
    for k = 1:ns
        s = strtrim(string(nsampChars((k-1)*8 + (1:8))));
        nsamp(k) = str2double(s);
    end

    fs = double(nsamp(chanIndexForFs)) / double(durRec);
end
