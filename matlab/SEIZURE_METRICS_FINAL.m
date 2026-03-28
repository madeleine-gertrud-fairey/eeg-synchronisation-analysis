%% ===================== FULL PIPELINE =====================

clc; close all;

%% ===================== 0) PATHS + SETTINGS (EDIT THESE) =====================
% C:\Users\Usuario\Desktop\TFG\TrainingP1Data\P1
% D:\PATIENT2\TrainingP2
% D:\PATIENT3\TrainingP3
dataDir   = "D:\PATIENT3\TrainingP3";
excelFile = "C:\Users\Usuario\Desktop\TFG\TrainingP1Data\SeizureDatesTraining.xlsx";
sheetName = "P3";
chanIdx   = [1 2];

% Segment around seizure onset (seconds)
pre_s  = 5*60;
post_s = 5*60;

% MPC band + filter
f_lo = 8; f_hi = 12; filt_order = 4;

% MPC time-series resolution
win_s  = 20;
step_s = 2;

% Metric windows relative to onset
baseline_win = [-300 -60];
seizure_win  = [-10   60];

% Trim edges to reduce filter/Hilbert edge artifacts
edge_trim_s = max(30, win_s);

% Recovery definition (relative to baselineMean)
recovery_frac = 0.95;
recovery_consecutive = 3;

fprintf("Settings: pre=%gs post=%gs | band=%g-%g Hz | win=%gs step=%gs | edge_trim=%gs\n", ...
    pre_s, post_s, f_lo, f_hi, win_s, step_s, edge_trim_s);

%% ===================== 1) READ EXCEL =====================
T = readtable(excelFile, "Sheet", sheetName, "VariableNamingRule","preserve");
onsetDT = toDatetimeRobust(T.("onset"));

good = ~isnat(onsetDT) & ~isnan(T.("EDF"));
T = T(good,:);
onsetDT = onsetDT(good);
fprintf("Parsed seizure rows: %d\n", height(T));

%% ===================== 2) BUILD EDF INDEX FROM HEADERS (FAST) =====================
d = dir(fullfile(dataDir, "*.edf"));
if isempty(d), error("No EDF files found in %s", dataDir); end

EDF = table('Size',[numel(d) 5], ...
    'VariableTypes',["string","datetime","double","double","datetime"], ...
    'VariableNames',["path","startDT","dur_s","fs","endDT"]);

for k = 1:numel(d)
    p = string(fullfile(d(k).folder, d(k).name));
    [startDT, fs, dur_s] = readEdfHeaderSummaryRaw(p, chanIdx);

    EDF.path(k)    = p;
    EDF.startDT(k) = startDT;
    EDF.dur_s(k)   = dur_s;
    EDF.fs(k)      = fs;
    EDF.endDT(k)   = startDT + seconds(dur_s);
end

EDF = sortrows(EDF, "startDT");

fprintf("Indexed %d EDF files.\n", height(EDF));
fprintf("First EDF: %s | %s\n", EDF.path(1), string(EDF.startDT(1)));
fprintf("Last  EDF: %s | %s\n", EDF.path(end), string(EDF.startDT(end)));

%% ===================== 3) RESOLVE EACH SEIZURE ONSET TO AN EDF =====================
shifts = [0, +1, -1, +2, -2]; % timezone/DST fallback in hours

resolved = table();
resolved.onset = onsetDT;
resolved.edfPath = strings(height(T),1);
resolved.edfStartDT = repmat(datetime(NaT), height(T), 1);
resolved.usedShiftHours = NaN(height(T),1);
resolved.offset_s = NaN(height(T),1);
resolved.fs = NaN(height(T),1);

miss = 0;

for i = 1:height(T)
    tOrig = onsetDT(i);

    hitIdx = [];
    usedShift = NaN;
    tUsed = tOrig;

    for s = shifts
        t = tOrig + hours(s);
        idx = find(t >= EDF.startDT & t <= EDF.endDT, 1, 'first');
        if ~isempty(idx)
            hitIdx = idx;
            usedShift = s;
            tUsed = t;
            break;
        end
    end

    if isempty(hitIdx)
        miss = miss + 1;
        continue;
    end

    resolved.edfPath(i) = EDF.path(hitIdx);
    resolved.edfStartDT(i) = EDF.startDT(hitIdx);
    resolved.usedShiftHours(i) = usedShift;
    resolved.offset_s(i) = seconds(tUsed - EDF.startDT(hitIdx));
    resolved.fs(i) = EDF.fs(hitIdx);
end

resolved = resolved(resolved.edfPath ~= "", :);

fprintf("\nResolved %d/%d seizures to EDF files (missed %d)\n", height(resolved), height(T), miss);
disp(groupsummary(resolved, "usedShiftHours"));

resolvedSeizures = resolved; %#ok<NASGU>
edfIndex = EDF; %#ok<NASGU>

%% ===================== 4) MPC METRICS (CACHE EDF + MPC ONCE PER SEIZURE) =====================
RS = resolved;

% Cache: read each EDF only once
sigCache = containers.Map('KeyType','char','ValueType','any');

out = table();
out.onset      = RS.onset;
out.edfPath    = RS.edfPath;
out.edfStartDT = RS.edfStartDT;
out.offset_s   = RS.offset_s;
out.fs         = RS.fs;

out.baselineMean = NaN(height(RS),1);
out.seizureMean  = NaN(height(RS),1);
out.deltaMean    = NaN(height(RS),1);

out.baselineMin  = NaN(height(RS),1);
out.seizureMin   = NaN(height(RS),1);
out.deltaMin     = NaN(height(RS),1);

out.baselineMax  = NaN(height(RS),1);
out.seizureMax   = NaN(height(RS),1);
out.deltaMax     = NaN(height(RS),1);

out.tSeizMin_s = NaN(height(RS),1);
out.recovery_s = NaN(height(RS),1);

skipped = false(height(RS),1);

for i = 1:height(RS)
    edfPath = char(RS.edfPath(i));
    fs = RS.fs(i);
    onsetOffset_s = RS.offset_s(i);

    % ---- load EDF signals once ----
    if isKey(sigCache, edfPath)
        X = sigCache(edfPath);
    else
        X = readEdfChannelsQuiet(edfPath, chanIdx); % samples x channels
        sigCache(edfPath) = X;
    end

    % ---- extract window around onset ----
    [seg, tRel] = extractSegmentWithTRel(X, fs, onsetOffset_s, pre_s, post_s);

    % ---- trim ends to reduce edge artifacts ----
    keep = (tRel >= (tRel(1)+edge_trim_s)) & (tRel <= (tRel(end)-edge_trim_s));
    if nnz(keep) >= round((win_s+5)*fs)
        seg  = seg(keep,:);
        tRel = tRel(keep);
    end

    % ---- compute MPC time series ONCE ----
    ts = computeMPCtimeseries(seg(:,1), seg(:,2), fs, f_lo, f_hi, filt_order, win_s, step_s, tRel);

    baseMask = (ts.t >= baseline_win(1)) & (ts.t <= baseline_win(2));
    seizMask = (ts.t >= seizure_win(1))  & (ts.t <= seizure_win(2));

    if nnz(baseMask) < 3 || nnz(seizMask) < 3
        skipped(i)=true;
        continue;
    end

    bVals = ts.mpc(baseMask);
    sVals = ts.mpc(seizMask);

    bMean = mean(bVals,'omitnan');
    sMean = mean(sVals,'omitnan');
    bMin  = min(bVals,[],'omitnan');
    [sMin, idxMin] = min(sVals,[],'omitnan');
    bMax  = max(bVals,[],'omitnan');
    sMax  = max(sVals,[],'omitnan');

    out.baselineMean(i) = bMean;
    out.seizureMean(i)  = sMean;
    out.deltaMean(i)    = bMean - sMean;

    out.baselineMin(i)  = bMin;
    out.seizureMin(i)   = sMin;
    out.deltaMin(i)     = bMin - sMin;

    out.baselineMax(i)  = bMax;
    out.seizureMax(i)   = sMax;
    out.deltaMax(i)     = bMax - sMax;

    tCandidates = ts.t(seizMask);
    out.tSeizMin_s(i) = tCandidates(idxMin);

    % recovery: first time >= frac*baselineMean for K consecutive samples after seizure minimum
    thr = recovery_frac * bMean;
    afterMin = (ts.t >= out.tSeizMin_s(i));
    tAfter = ts.t(afterMin);
    mAfter = ts.mpc(afterMin);

    runStart = findConsecutiveRun(mAfter >= thr, recovery_consecutive);
    if ~isnan(runStart)
        out.recovery_s(i) = tAfter(runStart) - out.tSeizMin_s(i);
    end
end

outA = out(~skipped,:);
fprintf("\nAnalysed %d/%d seizures (skipped %d). Unique EDFs read: %d\n", ...
    height(outA), height(out), sum(skipped), sigCache.Count);

mpc_metrics_table = outA; %#ok<NASGU>
assignin('base','mpc_metrics_table', outA);

%% ===================== 5) OUTLIERS + EXTREME + TYPICAL (MEAN) SUBPLOT =====================
dMean = outA.deltaMean;
negIdx = find(dMean < 0);
posIdx = find(dMean > 0);

[~, ord] = sort(dMean(posIdx), 'descend');
extIdx = posIdx(ord(1:min(2,numel(ord))));

remain = setdiff(posIdx, extIdx);
rng(1);
typIdx = remain(randperm(numel(remain), min(2,numel(remain))));

negPick = negIdx(1:min(2,numel(negIdx)));

showIdx = [negPick(:); extIdx(:); typIdx(:)];
showLab = [repmat("NEG deltaMean",numel(negPick),1);
           repmat("Extreme drop",numel(extIdx),1);
           repmat("Typical drop",numel(typIdx),1)];

figure('Color','w');
tiledlayout(2,3,'TileSpacing','compact','Padding','compact');

for kk = 1:numel(showIdx)
    ii = showIdx(kk);

    edfPath = char(outA.edfPath(ii));
    fs = outA.fs(ii);
    onsetOffset_s = outA.offset_s(ii);
    X = sigCache(edfPath);

    [seg, tRel] = extractSegmentWithTRel(X, fs, onsetOffset_s, pre_s, post_s);
    keep = (tRel >= (tRel(1)+edge_trim_s)) & (tRel <= (tRel(end)-edge_trim_s));
    if nnz(keep) >= round((win_s+5)*fs)
        seg = seg(keep,:); tRel = tRel(keep);
    end

    ts = computeMPCtimeseries(seg(:,1), seg(:,2), fs, f_lo, f_hi, filt_order, win_s, step_s, tRel);

    baseMask = (ts.t >= baseline_win(1)) & (ts.t <= baseline_win(2));
    seizMask = (ts.t >= seizure_win(1))  & (ts.t <= seizure_win(2));

    bMean = mean(ts.mpc(baseMask),'omitnan');
    sMean = mean(ts.mpc(seizMask),'omitnan');

    nexttile;
    plot(ts.t, ts.mpc, '-', 'LineWidth', 1); hold on;
    yline(bMean,'--','LineWidth',1);
    yline(sMean,':','LineWidth',1);
    xline(0,'-');
    xline(baseline_win(1),'--'); xline(baseline_win(2),'--');
    xline(seizure_win(1),'--');  xline(seizure_win(2),'--');
    grid on; ylim([0 1]);
    title(sprintf('%s | idx=%d | \\DeltaMean=%.3f', showLab(kk), ii, outA.deltaMean(ii)));
    xlabel('tRel (s)'); ylabel('MPC');
end
%% ===================== 5) OUTLIERS + EXTREME + TYPICAL (MIN) SUBPLOT =====================
dMin = outA.deltaMin;
negIdx = find(dMin < 0);
posIdx = find(dMin > 0);

[~, ord] = sort(dMin(posIdx), 'descend');
extIdx = posIdx(ord(1:min(2,numel(ord))));

remain = setdiff(posIdx, extIdx);
rng(1);
typIdx = remain(randperm(numel(remain), min(2,numel(remain))));

negPick = negIdx(1:min(2,numel(negIdx)));

showIdx = [negPick(:); extIdx(:); typIdx(:)];
showLab = [repmat("NEG deltaMin",numel(negPick),1);
           repmat("Extreme drop",numel(extIdx),1);
           repmat("Typical drop",numel(typIdx),1)];

figure('Color','w');
tiledlayout(2,3,'TileSpacing','compact','Padding','compact');

for kk = 1:numel(showIdx)
    ii = showIdx(kk);

    edfPath = char(outA.edfPath(ii));
    fs = outA.fs(ii);
    onsetOffset_s = outA.offset_s(ii);
    X = sigCache(edfPath);

    [seg, tRel] = extractSegmentWithTRel(X, fs, onsetOffset_s, pre_s, post_s);
    keep = (tRel >= (tRel(1)+edge_trim_s)) & (tRel <= (tRel(end)-edge_trim_s));
    if nnz(keep) >= round((win_s+5)*fs)
        seg = seg(keep,:); tRel = tRel(keep);
    end

    ts = computeMPCtimeseries(seg(:,1), seg(:,2), fs, f_lo, f_hi, filt_order, win_s, step_s, tRel);

    baseMask = (ts.t >= baseline_win(1)) & (ts.t <= baseline_win(2));
    seizMask = (ts.t >= seizure_win(1))  & (ts.t <= seizure_win(2));

    bMean = mean(ts.mpc(baseMask),'omitnan');
    sMean = mean(ts.mpc(seizMask),'omitnan');

    nexttile;
    plot(ts.t, ts.mpc, '-', 'LineWidth', 1); hold on;
    yline(bMean,'--','LineWidth',1);
    yline(sMean,':','LineWidth',1);
    xline(0,'-');
    xline(baseline_win(1),'--'); xline(baseline_win(2),'--');
    xline(seizure_win(1),'--');  xline(seizure_win(2),'--');
    grid on; ylim([0 1]);
    title(sprintf('%s | idx=%d | \\DeltaMin=%.3f', showLab(kk), ii, outA.deltaMin(ii)));
    xlabel('tRel (s)'); ylabel('MPC');
end
%% ===================== SUBPLOT: cols=deltaMean/deltaMin/deltaMax, rows=(1) dots (2) paired mins =====================
assert(exist('outA','var')==1, 'Need outA.');

% Use all valid rows (skip any NaNs cleanly)
good = isfinite(outA.deltaMean) & isfinite(outA.deltaMin) & isfinite(outA.deltaMax) & ...
       isfinite(outA.baselineMin) & isfinite(outA.seizureMin);

dMean = outA.deltaMean(good);
dMin  = outA.deltaMin(good);
dMax  = outA.deltaMax(good);

bMin  = outA.baselineMin(good);
sMin  = outA.seizureMin(good);

N = numel(dMean);
idx = 1:N;

% Limits (edit as you like)
xlim_delta = [1, max(2,N)];
ylim_delta = [-0.6 0.5];   % <-- change to your preferred limits
ylim_mpc   = [0 1];

figure('Color','w');
tiledlayout(2,3,'TileSpacing','compact','Padding','compact');

% ---------------- Row 1: DOT plots (deltas vs index) ----------------
nexttile(1);
scatter(idx, dMean, 25, 'filled'); grid on; yline(0,'-');
xlim(xlim_delta); ylim(ylim_delta);
title(sprintf('\\DeltaMean vs index (N=%d)', N));
xlabel('Seizure index'); ylabel('\Delta');

nexttile(2);
scatter(idx, dMin, 25, 'filled'); grid on; yline(0,'-');
xlim(xlim_delta); ylim(ylim_delta);
title(sprintf('\\DeltaMin vs index (N=%d)', N));
xlabel('Seizure index'); ylabel('\Delta');

nexttile(3);
scatter(idx, dMax, 25, 'filled'); grid on; yline(0,'-');
xlim(xlim_delta); ylim(ylim_delta);
title(sprintf('\\DeltaMax vs index (N=%d)', N));
xlabel('Seizure index'); ylabel('\Delta');

% ---------------- Row 2: PAIRED line plots (baselineMin vs seizureMin) ----------------
% Same paired-min values appear under each delta column; column title links them to the delta context.

% ---------------- Row 2: PAIRED line plots (baseline vs seizure) per metric ----------------

% Column 1: Mean pairing
nexttile(4);
plot([1 2], [outA.baselineMean(good) outA.seizureMean(good)]', '-o', ...
    'LineWidth', 0.8, 'MarkerSize', 4);
grid on;
xlim([0.8 2.2]); xticks([1 2]); xticklabels({'Baseline mean','Seizure mean'});
ylim(ylim_mpc);
title('Paired MEANS (BaselineMean → SeizureMean) | under \DeltaMean');
ylabel('MPC');

% Column 2: Min pairing
nexttile(5);
plot([1 2], [outA.baselineMin(good) outA.seizureMin(good)]', '-o', ...
    'LineWidth', 0.8, 'MarkerSize', 4);
grid on;
xlim([0.8 2.2]); xticks([1 2]); xticklabels({'Baseline min','Seizure min'});
ylim(ylim_mpc);
title('Paired MINS (BaselineMin → SeizureMin) | under \DeltaMin');
ylabel('MPC');

% Column 3: Max pairing
nexttile(6);
plot([1 2], [outA.baselineMax(good) outA.seizureMax(good)]', '-o', ...
    'LineWidth', 0.8, 'MarkerSize', 4);
grid on;
xlim([0.8 2.2]); xticks([1 2]); xticklabels({'Baseline max','Seizure max'});
ylim(ylim_mpc);
title('Paired MAXS (BaselineMax → SeizureMax) | under \DeltaMax');
ylabel('MPC');

%% ===================== 8) DELTAS vs SEIZURE TIME (SUBPLOTS) =====================
onsetDT_plot = outA.onset;
if ~isdatetime(onsetDT_plot)
    onsetDT_plot = outA.edfStartDT + seconds(outA.offset_s);
end

% optional: keep only finite values per metric
mMean = isfinite(outA.deltaMean) & isfinite(onsetDT_plot);
mMin  = isfinite(outA.deltaMin)  & isfinite(onsetDT_plot);
mMax  = isfinite(outA.deltaMax)  & isfinite(onsetDT_plot);

yl = [-0.6 0.6];

figure('Color','w');
tiledlayout(3,1,'TileSpacing','compact','Padding','compact');

nexttile;
scatter(onsetDT_plot(mMean), outA.deltaMean(mMean), 25, 'filled'); hold on;
yline(0,'-'); grid on;
ylim(yl);
ylabel('\DeltaMean');
title('\DeltaMean vs seizure time');

nexttile;
scatter(onsetDT_plot(mMin), outA.deltaMin(mMin), 25, 'filled'); hold on;
yline(0,'-'); grid on;
ylim(yl);
ylabel('\DeltaMin');
title('\DeltaMin vs seizure time');

nexttile;
scatter(onsetDT_plot(mMax), outA.deltaMax(mMax), 25, 'filled'); hold on;
yline(0,'-'); grid on;
ylim(yl);
ylabel('\DeltaMax');
title('\DeltaMax vs seizure time');
xlabel('Seizure onset time');

%% ===================== 9) RANDOM "FAKE SEIZURES" (55): UNIFORM vs CLUSTERED vs REAL =====================
Nrand = 55;

% Build ISI distribution from real seizures (seconds)
onsetDTs = sort(onsetDT_plot);
isi = seconds(diff(onsetDTs));
isi = isi(isfinite(isi) & isi>0);

% EDF durations known from header index
edfList = unique(outA.edfPath);
dur_s = nan(numel(edfList),1);
fsList = nan(numel(edfList),1);
for e = 1:numel(edfList)
    row = find(EDF.path==edfList(e),1,'first');
    dur_s(e) = EDF.dur_s(row);
    fsList(e)= EDF.fs(row);
end
prob = dur_s / sum(dur_s);

margin_s = pre_s + post_s + 2*edge_trim_s + win_s;

% Uniform random events
randUniform = table('Size',[Nrand 3], 'VariableTypes',{'string','double','double'}, ...
    'VariableNames',{'edfPath','fs','offset_s'});

for r = 1:Nrand
    e = randsample((1:numel(edfList))', 1, true, prob);
    D = dur_s(e);
    lo = margin_s; hi = D - margin_s;
    offset = (hi>lo) * (lo + (hi-lo)*rand) + (hi<=lo) * (D/2);

    randUniform.edfPath(r) = edfList(e);
    randUniform.fs(r)      = fsList(e);
    randUniform.offset_s(r)= offset;
end

% Clustered events (bursts using ISIs)
randClustered = table('Size',[Nrand 3], 'VariableTypes',{'string','double','double'}, ...
    'VariableNames',{'edfPath','fs','offset_s'});

rng(2);
r = 1;
while r <= Nrand
    e = randsample((1:numel(edfList))', 1, true, prob);
    D = dur_s(e);
    lo = margin_s; hi = D - margin_s;
    offset = (hi>lo) * (lo + (hi-lo)*rand) + (hi<=lo) * (D/2);

    for j = 1:10
        if r > Nrand, break; end
        if offset < lo || offset > hi, break; end

        randClustered.edfPath(r) = edfList(e);
        randClustered.fs(r)      = fsList(e);
        randClustered.offset_s(r)= offset;

        r = r + 1;
        if isempty(isi), break; end
        offset = offset + isi(randi(numel(isi)));
    end
end

R_uni = computeDeltasForEvents(randUniform, sigCache, chanIdx, pre_s, post_s, edge_trim_s, ...
    f_lo, f_hi, filt_order, win_s, step_s, baseline_win, seizure_win);

R_clu = computeDeltasForEvents(randClustered, sigCache, chanIdx, pre_s, post_s, edge_trim_s, ...
    f_lo, f_hi, filt_order, win_s, step_s, baseline_win, seizure_win);

%% ===================== 10) FIXED 3x3 HISTOGRAM GRID (LIMITS ALWAYS WORK) =====================
figure('Color','w');
tiledlayout(3,3,'TileSpacing','compact','Padding','compact');

edges = linspace(-0.6, 0.5, 30);

datasets = {
    outA.deltaMean,  R_uni.deltaMean,  R_clu.deltaMean;
    outA.deltaMin,   R_uni.deltaMin,   R_clu.deltaMin;
    outA.deltaMax,   R_uni.deltaMax,   R_clu.deltaMax
};

titles = {
    'Real \DeltaMean', 'Uniform \DeltaMean', 'Clustered \DeltaMean';
    'Real \DeltaMin',  'Uniform \DeltaMin',  'Clustered \DeltaMin';
    'Real \DeltaMax',  'Uniform \DeltaMax',  'Clustered \DeltaMax'
};

for rr = 1:3
    for cc = 1:3
        nexttile;
        x = datasets{rr,cc};
        x = x(isfinite(x));
        histogram(x, edges, 'FaceAlpha', 0.75);
        grid on;
        title(titles{rr,cc});
        xlabel('\Delta'); ylabel('Count');
        xlim([-0.6 0.5]);
        ylim([0 25]);
        set(gca,'XLimMode','manual','YLimMode','manual');
    end
end

%% ===================== 11) FIGURE 7: BOX + HIST OVERLAY + REAL SCATTER (LIMITS FIXED) =====================
edges2 = linspace(-0.6, 0.5, 30);
ylim_hist = [0 25];
ylim_sc   = [-0.6 0.5];

makeBox = @(realX, uniX, cluX) deal( ...
    [realX(:); uniX(:); cluX(:)], ...
    [repmat("Real",numel(realX),1); repmat("Uniform",numel(uniX),1); repmat("Clustered",numel(cluX),1)] );

metrics = {
    "DeltaMean", outA.deltaMean, R_uni.deltaMean, R_clu.deltaMean;
    "DeltaMin",  outA.deltaMin,  R_uni.deltaMin,  R_clu.deltaMin;
    "DeltaMax",  outA.deltaMax,  R_uni.deltaMax,  R_clu.deltaMax
};

figure('Color','w');
tiledlayout(3,3,'TileSpacing','compact','Padding','compact');

for m = 1:3
    name = metrics{m,1};
    realX = metrics{m,2}; uniX = metrics{m,3}; cluX = metrics{m,4};

    realX = realX(isfinite(realX));
    uniX  = uniX(isfinite(uniX));
    cluX  = cluX(isfinite(cluX));

    % Row 1: Boxplot
    nexttile(m);
    [Y,G] = makeBox(realX, uniX, cluX);
    boxchart(categorical(G), Y);
    grid on; ylim(ylim_sc);
    title(name + " | Boxplot");
    ylabel('\Delta');

    % Row 2: Histogram overlay (3 distributions on same axes)
    nexttile(m+3);
    histogram(realX, edges2, 'FaceAlpha',0.4); hold on;
    histogram(uniX,  edges2, 'FaceAlpha',0.4);
    histogram(cluX,  edges2, 'FaceAlpha',0.4);
    grid on;
    xlim([-0.6 0.5]); ylim(ylim_hist);
    set(gca,'XLimMode','manual','YLimMode','manual');
    title(name + " | Histogram overlay");
    xlabel('\Delta'); ylabel('Count');
    legend({'Real','Uniform','Clustered'},'Location','best'); legend boxoff;

    % Row 3: Real scatter (structure)
    nexttile(m+6);
    scatter(1:numel(realX), realX, 25, 'filled'); hold on; yline(0,'-');
    grid on; ylim(ylim_sc);
    title(name + " | Real scatter");
    xlabel('Index'); ylabel('\Delta');
end
%% ===================== 8) LEAD vs NON-LEAD (5h rule) + METRICS PLOTS ONLY =====================
%
% Rule:
%   NON-LEAD if time since previous seizure < 5 hours (global, sorted by Excel onset)
%   LEAD otherwise (first is lead)
%

assert(exist('outA','var')==1, "Need outA from section 4.");

% -------- Sort by seizure time (Excel onset) so the 5h rule is correct --------
outA = sortrows(outA, "onset");

% -------- Lead / Non-Lead classification --------
nonlead_thresh_h = 5;

isLead = true(height(outA),1);
for i = 2:height(outA)
    dt_h = hours(outA.onset(i) - outA.onset(i-1));
    if dt_h < nonlead_thresh_h
        isLead(i) = false;
    end
end
outA.isLead = isLead;

leadIdx    = outA.isLead;
nonLeadIdx = ~outA.isLead;

fprintf("Lead/Non-Lead (5h): N=%d | Lead=%d | NonLead=%d\n", ...
    height(outA), sum(leadIdx), sum(nonLeadIdx));

% -------- Use pipeline deltas: baseline - seizure --------
dMean = outA.deltaMean;
dMin  = outA.deltaMin;

% Colours (as you requested earlier)
colLead    = [1 0 0];   % red  = lead
colNonLead = [0 0 1];   % blue = non-lead

%% ===================== 10b) 3x3 HISTOGRAM GRID (REAL = STACKED Lead vs NonLead) =====================
% REQUIREMENT: outA.isLead exists (computed from your 5h rule)

assert(isfield(outA,'isLead') || any(strcmp(outA.Properties.VariableNames,'isLead')), ...
    'Need outA.isLead. Run the Lead/Non-Lead (5h rule) block BEFORE this section.');

figure('Color','w');
tiledlayout(3,3,'TileSpacing','compact','Padding','compact');

edges = linspace(-0.6, 0.5, 30);
yl = [0 25];
xl = [-0.6 0.5];

datasets = {
    outA.deltaMean,  R_uni.deltaMean,  R_clu.deltaMean;
    outA.deltaMin,   R_uni.deltaMin,   R_clu.deltaMin;
    outA.deltaMax,   R_uni.deltaMax,   R_clu.deltaMax
};

titles = {
    'Real (stacked) \DeltaMean', 'Uniform \DeltaMean', 'Clustered \DeltaMean';
    'Real (stacked) \DeltaMin',  'Uniform \DeltaMin',  'Clustered \DeltaMin';
    'Real (stacked) \DeltaMax',  'Uniform \DeltaMax',  'Clustered \DeltaMax'
};

for rr = 1:3
    for cc = 1:3
        nexttile;

        if cc == 1
            % ---------- REAL: stacked counts Lead vs NonLead ----------
            xAll = datasets{rr,cc};
            xAll = xAll(:);

            leadMask    = outA.isLead(:);
            nonLeadMask = ~outA.isLead(:);

            % Keep only finite per group
            xLead = xAll(leadMask  & isfinite(xAll));
            xNon  = xAll(nonLeadMask & isfinite(xAll));

            % Bin counts
            nLead = histcounts(xLead, edges);
            nNon  = histcounts(xNon,  edges);

            % Bin centers for bar()
            ctr = (edges(1:end-1) + edges(2:end))/2;
            bw  = ctr(2) - ctr(1);

            % Stacked bar (NonLead + Lead)
            bar(ctr, [nNon(:) nLead(:)], 1.0, 'stacked'); % width=1.0 fills bins
            grid on;
            xlim(xl); ylim(yl);
            title(titles{rr,cc});
            xlabel('\Delta'); ylabel('Count');
            legend({'Non-Lead','Lead'},'Location','best'); legend boxoff;

        else
            % ---------- UNIFORM / CLUSTERED: standard histogram ----------
            x = datasets{rr,cc};
            x = x(isfinite(x));
            histogram(x, edges, 'FaceAlpha', 0.75);
            grid on;
            title(titles{rr,cc});
            xlabel('\Delta'); ylabel('Count');
            xlim(xl); ylim(yl);
            set(gca,'XLimMode','manual','YLimMode','manual');
        end
    end
end
%% -------- (1) PAIRED LINES: Mean (Baseline vs Seizure) --------
figure('Color','w'); hold on;

for ii = find(nonLeadIdx)'
    plot([1 2], [outA.baselineMean(ii) outA.seizureMean(ii)], '-', ...
        'Color', colNonLead, 'LineWidth', 1.5);
end
for ii = find(leadIdx)'
    plot([1 2], [outA.baselineMean(ii) outA.seizureMean(ii)], '-', ...
        'Color', colLead, 'LineWidth', 1.5);
end

xlim([0.8 2.2]); xticks([1 2]); xticklabels({'Baseline','Seizure'});
ylabel('MPC Mean');
title(sprintf('Mean MPC (Baseline vs Seizure) | NON-LEAD if < %gh', nonlead_thresh_h));
grid on;

h1 = plot(nan,nan,'-','Color',colNonLead,'LineWidth',2);
h2 = plot(nan,nan,'-','Color',colLead,'LineWidth',2);
legend([h1 h2], {'Non-Lead','Lead'}, 'Location','best'); legend boxoff;

%% -------- (2) PAIRED LINES: Min (Baseline vs Seizure) --------
figure('Color','w'); hold on;

for ii = find(nonLeadIdx)'
    plot([1 2], [outA.baselineMin(ii) outA.seizureMin(ii)], '-', ...
        'Color', colNonLead, 'LineWidth', 1.5);
end
for ii = find(leadIdx)'
    plot([1 2], [outA.baselineMin(ii) outA.seizureMin(ii)], '-', ...
        'Color', colLead, 'LineWidth', 1.5);
end

xlim([0.8 2.2]); xticks([1 2]); xticklabels({'Baseline','Seizure'});
ylabel('MPC Min');
title(sprintf('Min MPC (Baseline vs Seizure) | NON-LEAD if < %gh', nonlead_thresh_h));
grid on;

h1 = plot(nan,nan,'-','Color',colNonLead,'LineWidth',2);
h2 = plot(nan,nan,'-','Color',colLead,'LineWidth',2);
legend([h1 h2], {'Non-Lead','Lead'}, 'Location','best'); legend boxoff;

%% -------- (3) DOT DISTRIBUTIONS: ΔMean and ΔMin (baseline - seizure) --------
figure('Color','w'); hold on;
scatter(ones(sum(nonLeadIdx),1), dMean(nonLeadIdx), 60, colNonLead, 'filled');
scatter(2*ones(sum(leadIdx),1),    dMean(leadIdx),    60, colLead,    'filled');
yline(0,'-','LineWidth',1);
xlim([0.5 2.5]); xticks([1 2]); xticklabels({'Non-Lead','Lead'});
ylabel('\Delta Mean MPC (Baseline - Seizure)');
title(sprintf('\\DeltaMean distribution | NON-LEAD if < %gh', nonlead_thresh_h));
grid on;
h1 = scatter(nan,nan,60,colNonLead,'filled'); h2 = scatter(nan,nan,60,colLead,'filled');
legend([h1 h2], {'Non-Lead','Lead'}, 'Location','best'); legend boxoff;

figure('Color','w'); hold on;
scatter(ones(sum(nonLeadIdx),1), dMin(nonLeadIdx), 60, colNonLead, 'filled');
scatter(2*ones(sum(leadIdx),1),    dMin(leadIdx),    60, colLead,    'filled');
yline(0,'-','LineWidth',1);
xlim([0.5 2.5]); xticks([1 2]); xticklabels({'Non-Lead','Lead'});
ylabel('\Delta Min MPC (Baseline - Seizure)');
title(sprintf('\\DeltaMin distribution | NON-LEAD if < %gh', nonlead_thresh_h));
grid on;
h1 = scatter(nan,nan,60,colNonLead,'filled'); h2 = scatter(nan,nan,60,colLead,'filled');
legend([h1 h2], {'Non-Lead','Lead'}, 'Location','best'); legend boxoff;

%% -------- (4) Δ vs INDEX (filled dots) --------
idx = (1:height(outA))';

figure('Color','w'); hold on;
yline(0,'-','LineWidth',1);
scatter(idx(nonLeadIdx), dMean(nonLeadIdx), 50, colNonLead, 'filled');
scatter(idx(leadIdx),    dMean(leadIdx),    50, colLead,    'filled');
grid on;
xlabel('Seizure index (sorted by onset time)');
ylabel('\Delta Mean MPC (Baseline - Seizure)');
title(sprintf('\\DeltaMean vs index (N=%d) | NON-LEAD if < %gh', height(outA), nonlead_thresh_h));
h1 = scatter(nan,nan,50,colNonLead,'filled'); h2 = scatter(nan,nan,50,colLead,'filled');
legend([h1 h2], {'Non-Lead','Lead'}, 'Location','best'); legend boxoff;

figure('Color','w'); hold on;
yline(0,'-','LineWidth',1);
scatter(idx(nonLeadIdx), dMin(nonLeadIdx), 50, colNonLead, 'filled');
scatter(idx(leadIdx),    dMin(leadIdx),    50, colLead,    'filled');
grid on;
xlabel('Seizure index (sorted by onset time)');
ylabel('\Delta Min MPC (Baseline - Seizure)');
title(sprintf('\\DeltaMin vs index (N=%d) | NON-LEAD if < %gh', height(outA), nonlead_thresh_h));
h1 = scatter(nan,nan,50,colNonLead,'filled'); h2 = scatter(nan,nan,50,colLead,'filled');
legend([h1 h2], {'Non-Lead','Lead'}, 'Location','best'); legend boxoff;


%% ===================== LOCAL FUNCTIONS =====================

function X = readEdfChannelsQuiet(edfPath, chanIdx)
% Reads EDF channels and returns [samples x channels].
% Suppresses printing from edfreadUntilDone using evalc.
try
    tmp = struct();
    evalc("[tmp.hdr, tmp.rec] = edfreadUntilDone('" + string(edfPath) + "', 'targetSignals', [" + num2str(chanIdx) + "]);");
    X = double(tmp.rec.');
catch
    % Fallback: try MATLAB edfread if available
    [~, tt] = edfread(edfPath);
    vars = tt.Properties.VariableNames;
    x1 = double(tt.(vars{chanIdx(1)})); x1 = x1(:);
    x2 = double(tt.(vars{chanIdx(2)})); x2 = x2(:);
    X = [x1 x2];
end
end

function Rout = computeDeltasForEvents(Tevents, sigCache, chanIdx, pre_s, post_s, edge_trim_s, ...
    f_lo, f_hi, filt_order, win_s, step_s, baseline_win, seizure_win)

N = height(Tevents);
Rout = table();
Rout.deltaMean = nan(N,1);
Rout.deltaMin  = nan(N,1);
Rout.deltaMax  = nan(N,1);

for i = 1:N
    edfPath = char(Tevents.edfPath(i));
    fs      = Tevents.fs(i);
    onsetOffset_s = Tevents.offset_s(i);

    if isKey(sigCache, edfPath)
        X = sigCache(edfPath);
    else
        X = readEdfChannelsQuiet(edfPath, chanIdx);
        sigCache(edfPath) = X;
    end

    [seg, tRel] = extractSegmentWithTRel(X, fs, onsetOffset_s, pre_s, post_s);

    keep = (tRel >= (tRel(1)+edge_trim_s)) & (tRel <= (tRel(end)-edge_trim_s));
    if nnz(keep) < round((win_s+5)*fs), continue; end
    seg  = seg(keep,:);
    tRel = tRel(keep);

    ts = computeMPCtimeseries(seg(:,1), seg(:,2), fs, f_lo, f_hi, filt_order, win_s, step_s, tRel);

    baseMask = (ts.t >= baseline_win(1)) & (ts.t <= baseline_win(2));
    seizMask = (ts.t >= seizure_win(1))  & (ts.t <= seizure_win(2));
    if nnz(baseMask)<3 || nnz(seizMask)<3, continue; end

    bVals = ts.mpc(baseMask);
    sVals = ts.mpc(seizMask);

    Rout.deltaMean(i) = mean(bVals,'omitnan') - mean(sVals,'omitnan');
    Rout.deltaMin(i)  = min(bVals,[],'omitnan') - min(sVals,[],'omitnan');
    Rout.deltaMax(i)  = max(bVals,[],'omitnan') - max(sVals,[],'omitnan');
end
end

function dt = toDatetimeRobust(x)
if isdatetime(x), dt = x; return; end
if iscell(x), x = string(x); end
if isnumeric(x)
    try, dt = datetime(x, "ConvertFrom","excel"); return; catch, end
end
if isstring(x) || ischar(x)
    x = string(x);
    dt = NaT(size(x));
    fmts = [
        "dd/MM/yyyy HH:mm:ss"
        "dd/MM/yyyy HH:mm"
        "dd/MM/yyyy"
        "dd-MMM-yyyy HH:mm:ss"
        "dd-MMM-yyyy HH:mm"
        "dd-MMM-yyyy"
        "yyyy-MM-dd HH:mm:ss"
        "yyyy-MM-dd"
    ];
    for k = 1:numel(fmts)
        tmp = NaT(size(x));
        try, tmp = datetime(x, "InputFormat", fmts(k), "Locale","en_US"); catch, end
        dt(isnat(dt) & ~isnat(tmp)) = tmp(isnat(dt) & ~isnat(tmp));
    end
    return;
end
try, dt = datetime(x); catch, dt = NaT(size(x)); end
end

function [startDT, fs, dur_s] = readEdfHeaderSummaryRaw(edfPath, chanIdx)
fid = fopen(edfPath, 'r', 'ieee-le');
if fid < 0, error("Cannot open EDF: %s", edfPath); end
cleanup = onCleanup(@() fclose(fid));

fixed = fread(fid, 256, '*char')';

startDateStr = strtrim(string(fixed(169:176)));
startTimeStr = replace(strtrim(string(fixed(177:184))), ".", ":");

d = datetime(startDateStr, "InputFormat","dd.MM.yy", "Locale","en_US");
t = datetime(startTimeStr, "InputFormat","HH:mm:ss", "Locale","en_US");
startDT = dateshift(d,"start","day") + timeofday(t);

nRec   = str2double(strtrim(string(fixed(237:244))));
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

fs = double(nsamp(chanIdx(1))) / double(durRec);
dur_s = double(nRec) * double(durRec);
end

function [seg, tRel] = extractSegmentWithTRel(x, fs, onsetOffset_s, pre_s, post_s)
n = size(x,1);
sOn = round(onsetOffset_s * fs) + 1;
s0  = max(1, sOn - round(pre_s*fs));
s1  = min(n, sOn + round(post_s*fs));
seg = x(s0:s1,:);
iOn = sOn - s0 + 1;
tRel = ((1:size(seg,1)) - iOn) / fs;
tRel = tRel(:);
end

function ts = computeMPCtimeseries(x1, x2, fs, f_lo, f_hi, order, win_s, step_s, tRel)
y1 = bandpass_filter(x1, fs, f_lo, f_hi, order);
y2 = bandpass_filter(x2, fs, f_lo, f_hi, order);

p1 = angle(hilbert(y1));
p2 = angle(hilbert(y2));
dphi = angle(exp(1j*(p1 - p2)));

n = numel(x1);
win_samp  = max(2, round(win_s*fs));
step_samp = max(1, round(step_s*fs));
starts = 1:step_samp:(n-win_samp+1);

mpc = zeros(numel(starts),1);
t   = zeros(numel(starts),1);

for ii = 1:numel(starts)
    s = starts(ii); e = s + win_samp - 1;
    mpc(ii) = abs(mean(exp(1j*dphi(s:e))));
    mid = round((s+e)/2);
    t(ii) = tRel(mid);
end
ts = struct("t",t,"mpc",mpc);
end

function y = bandpass_filter(x, fs, f_lo, f_hi, order)
nyq = 0.5*fs;
[b,a] = butter(order, [f_lo/nyq, f_hi/nyq], "bandpass");
y = filtfilt(b,a,double(x));
end

function idx = findConsecutiveRun(v, L)
idx = NaN; count = 0;
for i = 1:numel(v)
    if v(i)
        count = count + 1;
        if count >= L
            idx = i - L + 1;
            return;
        end
    else
        count = 0;
    end
end
end
