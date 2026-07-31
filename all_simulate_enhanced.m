%% RIoT-2 attitude reconstruction, ADCS quality, and access analysis
% Enhanced portfolio version for the 2026-04-22 ADCS log.
% Outputs figures, CSV summaries, and a Markdown report.

clear; clc; close all;

%% 1. Paths and options
thisFile = mfilename("fullpath");
outDir = fileparts(thisFile);
tleDir = fullfile(outDir, "TLE");
outputDir = fullfile(outDir, "outputs_enhanced");
figureOutDir = fullfile(outputDir, "figures");

if ~isfolder(outputDir)
    mkdir(outputDir);
end
if ~isfolder(figureOutDir)
    mkdir(figureOutDir);
end

show3DViewer = true;
sampleTimeSeconds = 10;
sensorViewAngle = 30;
stabilityRateThresholdDegps = 0.08;
pointingThresholdDeg = 25;

% Assumed spacecraft body axes. Calibrate these with spacecraft CAD or ADCS
% documentation when available.
nadirBoresightBody = [0; 0; 1];
solarArrayNormalBody = [1; 0; 0];
groundStationAntennaBody = [0; 0; 1];
gsLatDeg = 24.96;
gsLonDeg = 121.22;
gsAltMeters = 0;

adcsInput = fullfile(tleDir, "1776831872 4-30-RIoT-2-log.csv");
tleInput = [
    "RIoT2"
    "1 10256U 26067BL  26111.96345576  .00003006  00000-0  29964-3 0  9993"
    "2 10256  97.7513  71.1114 0003936  79.4998 280.6664 14.92316398  3348"
];

timeColumn = "Unix time";
quatColumns = [
    "Estimated ORC quaternion Q0"
    "Estimated ORC quaternion Q1"
    "Estimated ORC quaternion Q2"
    "Estimated ORC quaternion Q3"
];
rateColumns = [
    "Estimated body rate (ORC) X component (degps)"
    "Estimated body rate (ORC) Y component (degps)"
    "Estimated body rate (ORC) Z component (degps)"
];

quatOrder = "xyzw";
attitudeDirection = "bodyToOrc";

%% 2. Resolve and read inputs
tleFile = writeTempTleFile(tleInput);
satelliteName = tleNameFromLines(tleInput);

if ~isfile(adcsInput)
    error("Cannot find ADCS log: %s", adcsInput);
end

qRaw = readtable(adcsInput, "VariableNamingRule", "preserve");
clockUnix = qRaw{:, timeColumn};
tUtc = datetime(clockUnix, "ConvertFrom", "posixtime", "TimeZone", "UTC");

validTime = ~isnat(tUtc);
qRaw = qRaw(validTime, :);
tUtc = tUtc(validTime);
[tUtc, uniqueIdx] = unique(tUtc, "stable");
qRaw = qRaw(uniqueIdx, :);
n = height(qRaw);

qData = double([qRaw{:, quatColumns(1)}, qRaw{:, quatColumns(2)}, ...
    qRaw{:, quatColumns(3)}, qRaw{:, quatColumns(4)}]);

if quatOrder == "xyzw"
    qBodyOrc = qData;
elseif quatOrder == "wxyz"
    qBodyOrc = [qData(:,2), qData(:,3), qData(:,4), qData(:,1)];
else
    error("quatOrder must be 'xyzw' or 'wxyz'.");
end

qNormRaw = vecnorm(qBodyOrc, 2, 2);
validQ = all(isfinite(qBodyOrc), 2) & qNormRaw > 0;
qBodyOrc(validQ, :) = qBodyOrc(validQ, :) ./ qNormRaw(validQ);

bodyRateDegps = double([qRaw{:, rateColumns(1)}, qRaw{:, rateColumns(2)}, qRaw{:, rateColumns(3)}]);
totalRateDegps = vecnorm(bodyRateDegps, 2, 2);
stableRateFlag = totalRateDegps <= stabilityRateThresholdDegps;

fprintf("Loaded %d ADCS rows from %s\n", n, adcsInput);
fprintf("Telemetry span: %s to %s UTC\n", string(tUtc(1)), string(tUtc(end)));
fprintf("Valid quaternion coverage: %.1f %%\n", 100 * mean(validQ));

%% 3. Satellite scenario
sc = satelliteScenario(tUtc(1), tUtc(end), sampleTimeSeconds);
sat = satellite(sc, tleFile, "Name", satelliteName);

try
    sat.Visual3DModel = "SmallSat.glb";
    sat.Visual3DModelScale = 1;
catch
    warning("SmallSat.glb was not available. Continuing with default satellite marker.");
end

coordinateAxes(sat, "Scale", 2);
gs = groundStation(sc, gsLatDeg, gsLonDeg, "Name", "Zhongli_GS");
sensor = conicalSensor(sat, "MaxViewAngle", sensorViewAngle, "Name", "Sat_Sensor");
tx = transmitter(sat, "Frequency", 2.4e9, "Power", 10);
rx = receiver(gs, "Name", "Zhongli_Rx");
lnk = link(tx, rx);
acSatGs = access(sat, gs);
acSensorGs = access(sensor, gs);

%% 4. Attitude reconstruction
qInertial = nan(n, 4);
RBodyToEciFlat = nan(n, 9);
rEciUnit = nan(n, 3);
sunHatEci = nan(n, 3);

for k = find(validQ).'
    [rk, vk] = states(sat, tUtc(k), "CoordinateFrame", "inertial");
    rk = rk(:);
    vk = vk(:);

    R_eci_from_orc = eciFromLvlh(rk, vk);
    R_orc_from_body = quatScalarLastToDcm(qBodyOrc(k, :));

    if attitudeDirection == "bodyToOrc"
        R_eci_from_body = R_eci_from_orc * R_orc_from_body;
    else
        R_eci_from_body = R_eci_from_orc * R_orc_from_body.';
    end

    RBodyToEciFlat(k, :) = reshape(R_eci_from_body, 1, []);
    qInertial(k, :) = dcm2quat(R_eci_from_body.');
    rEciUnit(k, :) = (rk ./ norm(rk)).';
    sunHatEci(k, :) = approximateSunEciUnit(tUtc(k)).';
end

attitudeTT = timetable(tUtc(validQ), qInertial(validQ, :));
pointAt(sat, attitudeTT);

%% 5. Orbit geometry and pointing errors
elevationDeg = nan(n, 1);
azimuthDeg = nan(n, 1);
rangeKm = nan(n, 1);
nadirErrorDeg = nan(n, 1);
sunErrorDeg = nan(n, 1);
gsPointingErrorDeg = nan(n, 1);
sunlitFlag = false(n, 1);

for k = 1:n
    [rk, ~] = states(sat, tUtc(k), "CoordinateFrame", "inertial");
    rk = rk(:);
    if any(~isfinite(rEciUnit(k, :)))
        rEciUnit(k, :) = (rk ./ norm(rk)).';
    end
    if any(~isfinite(sunHatEci(k, :)))
        sunHatEci(k, :) = approximateSunEciUnit(tUtc(k)).';
    end
    sunlitFlag(k) = ~isInCylindricalEclipse(rk, sunHatEci(k, :).');

    try
        [azimuthDeg(k), elevationDeg(k), rangeMeters] = aer(gs, sat, tUtc(k));
        rangeKm(k) = rangeMeters / 1000;
    catch
        [azimuthDeg(k), elevationDeg(k), rangeKm(k)] = deal(nan);
    end

    if validQ(k)
        R_eci_from_body = reshape(RBodyToEciFlat(k, :), 3, 3);
        nadirErrorDeg(k) = angularSeparationDeg(R_eci_from_body * nadirBoresightBody, -rEciUnit(k, :).');
        if sunlitFlag(k)
            sunErrorDeg(k) = angularSeparationDeg(R_eci_from_body * solarArrayNormalBody, sunHatEci(k, :).');
        end

        if isfinite(azimuthDeg(k)) && isfinite(elevationDeg(k)) && isfinite(rangeKm(k))
            gsDirectionEci = groundStationDirectionEci(gsLatDeg, gsLonDeg, gsAltMeters, sat, tUtc(k));
            gsPointingErrorDeg(k) = angularSeparationDeg(R_eci_from_body * groundStationAntennaBody, gsDirectionEci);
        end
    end
end

accessSatGsTable = accessIntervals(acSatGs);
accessSensorGsTable = accessIntervals(acSensorGs);
linkTable = linkIntervals(lnk);

telemetryVisibleFlag = elevationDeg > 0;
telemetrySensorCandidateFlag = telemetryVisibleFlag & gsPointingErrorDeg <= sensorViewAngle;

%% 6. Summary tables
summary = table( ...
    ["rows"; "span_minutes"; "quaternion_coverage_percent"; "quaternion_norm_median"; ...
     "quaternion_norm_max_abs_error"; "body_rate_rms_degps"; "body_rate_max_degps"; ...
     "stable_rate_fraction_percent"; "median_nadir_error_deg"; "median_sun_error_deg"; ...
     "median_gs_pointing_error_deg"; "visible_telemetry_samples"; "sensor_candidate_samples"; ...
     "geometric_access_windows"; "sensor_limited_access_windows"; "link_windows"], ...
    [n; minutes(tUtc(end) - tUtc(1)); 100 * mean(validQ); median(qNormRaw, "omitnan"); ...
     max(abs(qNormRaw(validQ) - 1), [], "omitnan"); rmsFinite(totalRateDegps); max(totalRateDegps, [], "omitnan"); ...
     100 * mean(stableRateFlag, "omitnan"); median(nadirErrorDeg, "omitnan"); median(sunErrorDeg, "omitnan"); ...
     median(gsPointingErrorDeg, "omitnan"); nnz(telemetryVisibleFlag); nnz(telemetrySensorCandidateFlag); ...
     height(accessSatGsTable); height(accessSensorGsTable); height(linkTable)], ...
    VariableNames=["Metric", "Value"]);

analysisTable = table(tUtc(:), qNormRaw(:), validQ(:), bodyRateDegps(:,1), bodyRateDegps(:,2), bodyRateDegps(:,3), ...
    totalRateDegps(:), stableRateFlag(:), elevationDeg(:), azimuthDeg(:), rangeKm(:), sunlitFlag(:), ...
    nadirErrorDeg(:), sunErrorDeg(:), gsPointingErrorDeg(:), telemetryVisibleFlag(:), telemetrySensorCandidateFlag(:), ...
    VariableNames=["UtcTime", "QuaternionNormRaw", "ValidQuaternion", "BodyRateXDegps", "BodyRateYDegps", ...
    "BodyRateZDegps", "TotalRateDegps", "StableRateFlag", "ElevationDeg", "AzimuthDeg", "RangeKm", ...
    "Sunlit", "NadirErrorDeg", "SunErrorDeg", "GroundStationPointingErrorDeg", "VisibleFromZhongli", ...
    "SensorCandidateAtTelemetryTime"]);

writetable(summary, fullfile(outputDir, "summary_metrics.csv"));
writetable(analysisTable, fullfile(outputDir, "attitude_access_timeseries.csv"));
writetable(accessSatGsTable, fullfile(outputDir, "geometric_access_windows.csv"));
writetable(accessSensorGsTable, fullfile(outputDir, "sensor_limited_access_windows.csv"));
writetable(linkTable, fullfile(outputDir, "communication_link_windows.csv"));

disp("====== Enhanced summary metrics ======");
disp(summary);

%% 7. Figures
fig = readableFigure("Quaternion quality", [80 80 1180 760]);
tl = tiledlayout(fig, 2, 1, "TileSpacing", "compact", "Padding", "compact");
title(tl, "ADCS quaternion quality check", "FontSize", 18, "FontWeight", "bold");
ax = nexttile(tl);
plot(tUtc, qNormRaw, "Color", [0.13 0.37 0.76]); hold on;
yline(1.0, "--", "Color", [0.25 0.25 0.25]);
yline(1.01, ":", "Color", [0.65 0.20 0.16]);
yline(0.99, ":", "Color", [0.65 0.20 0.16]);
ylabel("Raw norm");
title("Quaternion norm before normalization", "FontSize", 14);
styleTimeAxis(ax, tUtc);
ax = nexttile(tl);
plot(tUtc, qBodyOrc(:,1)); hold on;
plot(tUtc, qBodyOrc(:,2));
plot(tUtc, qBodyOrc(:,3));
plot(tUtc, qBodyOrc(:,4));
ylabel("Normalized quaternion");
xlabel("UTC time");
legend("qx", "qy", "qz", "qw", "Location", "northoutside", "Orientation", "horizontal");
title("Normalized quaternion components", "FontSize", 14);
styleTimeAxis(ax, tUtc);
saveReadableFigure(fig, figureOutDir, "01_quaternion_quality");

fig = readableFigure("Body rate stability", [100 100 1180 760]);
tl = tiledlayout(fig, 2, 1, "TileSpacing", "compact", "Padding", "compact");
title(tl, "ADCS body-rate stability", "FontSize", 18, "FontWeight", "bold");
ax = nexttile(tl);
plot(tUtc, bodyRateDegps(:,1)); hold on;
plot(tUtc, bodyRateDegps(:,2));
plot(tUtc, bodyRateDegps(:,3));
ylabel("Rate (deg/s)");
legend("X", "Y", "Z", "Location", "northoutside", "Orientation", "horizontal");
title("Estimated body rates in ORC frame", "FontSize", 14);
styleTimeAxis(ax, tUtc);
ax = nexttile(tl);
plot(tUtc, totalRateDegps, "Color", [0.20 0.55 0.30]); hold on;
yline(stabilityRateThresholdDegps, "--", "Color", [0.65 0.20 0.16]);
ylabel("Total rate (deg/s)");
xlabel("UTC time");
title("Total angular rate and stability threshold", "FontSize", 14);
styleTimeAxis(ax, tUtc);
saveReadableFigure(fig, figureOutDir, "02_body_rate_stability");

fig = readableFigure("Pointing errors", [120 120 1180 860]);
tl = tiledlayout(fig, 3, 1, "TileSpacing", "compact", "Padding", "compact");
title(tl, "Telemetry-derived pointing performance", "FontSize", 18, "FontWeight", "bold");
ax = nexttile(tl);
plot(tUtc, nadirErrorDeg, "Color", [0.13 0.37 0.76]); hold on;
yline(pointingThresholdDeg, "--", "Color", [0.65 0.20 0.16]);
ylabel("Error (deg)");
title("Nadir / Earth-pointing error", "FontSize", 14);
styleTimeAxis(ax, tUtc);
ax = nexttile(tl);
plot(tUtc, sunErrorDeg, "Color", [0.90 0.42 0.10]); hold on;
yline(pointingThresholdDeg, "--", "Color", [0.65 0.20 0.16]);
ylabel("Error (deg)");
title("Sun-pointing error during sunlight", "FontSize", 14);
styleTimeAxis(ax, tUtc);
ax = nexttile(tl);
plot(tUtc, gsPointingErrorDeg, "Color", [0.20 0.55 0.30]); hold on;
yline(sensorViewAngle, "--", "Color", [0.65 0.20 0.16]);
ylabel("Error (deg)");
xlabel("UTC time");
title("Ground-station antenna pointing error", "FontSize", 14);
styleTimeAxis(ax, tUtc);
saveReadableFigure(fig, figureOutDir, "03_pointing_errors");

fig = readableFigure("Access timeline", [140 140 1180 800]);
tl = tiledlayout(fig, 3, 1, "TileSpacing", "compact", "Padding", "compact");
title(tl, "Zhongli ground-station geometry and access checks", "FontSize", 18, "FontWeight", "bold");
ax = nexttile(tl);
plot(tUtc, elevationDeg, "Color", [0.13 0.37 0.76]); hold on;
yline(0, "--", "Color", [0.25 0.25 0.25]);
ylabel("Elevation (deg)");
title("Elevation angle from Zhongli", "FontSize", 14);
styleTimeAxis(ax, tUtc);
ax = nexttile(tl);
plot(tUtc, rangeKm, "Color", [0.34 0.49 0.18]);
ylabel("Range (km)");
title("Range from Zhongli", "FontSize", 14);
styleTimeAxis(ax, tUtc);
ax = nexttile(tl);
stairs(tUtc, double(telemetryVisibleFlag), "Color", [0.13 0.37 0.76], "LineWidth", 2.0); hold on;
stairs(tUtc, double(telemetrySensorCandidateFlag), "Color", [0.90 0.42 0.10], "LineWidth", 2.0);
ylim([-0.1 1.1]);
yticks([0 1]);
yticklabels(["No", "Yes"]);
ylabel("Available");
xlabel("UTC time");
legend("Elevation > 0", "Visible and sensor aligned", "Location", "northoutside", "Orientation", "horizontal");
title("Telemetry-time access candidates", "FontSize", 14);
styleTimeAxis(ax, tUtc);
saveReadableFigure(fig, figureOutDir, "04_access_timeline");

%% 8. Markdown report
reportFile = fullfile(outputDir, "enhanced_analysis_report.md");
writeMarkdownReport(reportFile, summary, figureOutDir, adcsInput, tleInput);
fprintf("Enhanced outputs written to: %s\n", outputDir);

%% 9. Optional 3D playback
if isequal(show3DViewer, true)
    viewer = satelliteScenarioViewer(sc, "ShowDetails", true);
    show(sat);
    show(gs);
    camtarget(viewer, sat);
    play(sc);
end

%% Local functions
function fig = readableFigure(name, position)
    fig = figure("Name", name, "NumberTitle", "off", "Color", "w", ...
        "WindowStyle", "normal", "Units", "pixels", "Position", position);
end

function saveReadableFigure(fig, outDir, baseName)
    if ~isfolder(outDir)
        mkdir(outDir);
    end
    polishFigure(fig);
    drawnow;
    exportgraphics(fig, fullfile(outDir, baseName + ".png"), "Resolution", 180, "BackgroundColor", "white");
    close(fig);
end

function polishFigure(fig)
    axs = findall(fig, "Type", "Axes");
    for i = 1:numel(axs)
        axs(i).Color = "w";
        axs(i).XColor = [0.12 0.12 0.12];
        axs(i).YColor = [0.12 0.12 0.12];
        axs(i).GridColor = [0.82 0.82 0.82];
        axs(i).GridAlpha = 0.45;
        axs(i).LineWidth = 1.0;
        axs(i).FontSize = 12;
        axs(i).Title.FontWeight = "bold";
    end
    legends = findall(fig, "Type", "Legend");
    for i = 1:numel(legends)
        legends(i).Color = "w";
        legends(i).TextColor = [0.08 0.08 0.08];
    end
end

function styleTimeAxis(ax, t)
    grid(ax, "on");
    box(ax, "on");
    if ~isempty(t)
        xlim(ax, [t(1) t(end)]);
    end
end

function tleFile = writeTempTleFile(tleLines)
    tleLines = string(tleLines);
    tleLines = tleLines(strlength(strtrim(tleLines)) > 0);
    tleFile = fullfile(tempdir, "riotr2_enhanced_selected.tle");
    writelines(tleLines, tleFile);
end

function satelliteName = tleNameFromLines(tleLines)
    satelliteName = matlab.lang.makeValidName(strtrim(string(tleLines(1))));
end

function R = eciFromLvlh(r, v)
    r = r(:);
    v = v(:);
    z = -r / norm(r);
    h = cross(r, v);
    y = -h / norm(h);
    x = cross(y, z);
    R = [x, y, z];
end

function dcm = quatScalarLastToDcm(q)
    qx = q(1); qy = q(2); qz = q(3); qw = q(4);
    dcm = [1 - 2*(qy^2 + qz^2), 2*(qx*qy - qw*qz), 2*(qx*qz + qw*qy); ...
           2*(qx*qy + qw*qz), 1 - 2*(qx^2 + qz^2), 2*(qy*qz - qw*qx); ...
           2*(qx*qz - qw*qy), 2*(qy*qz + qw*qx), 1 - 2*(qx^2 + qy^2)];
end

function angleDeg = angularSeparationDeg(a, b)
    if any(~isfinite(a)) || any(~isfinite(b)) || norm(a) == 0 || norm(b) == 0
        angleDeg = nan;
        return;
    end
    c = dot(a(:), b(:)) / (norm(a) * norm(b));
    c = min(1, max(-1, c));
    angleDeg = acosd(c);
end

function sunHat = approximateSunEciUnit(t)
    jd = juliandate(t);
    n = jd - 2451545.0;
    L = mod(280.460 + 0.9856474 * n, 360);
    g = deg2rad(mod(357.528 + 0.9856003 * n, 360));
    lambda = deg2rad(L + 1.915 * sin(g) + 0.020 * sin(2 * g));
    epsilon = deg2rad(23.439 - 0.0000004 * n);
    sunHat = [cos(lambda); cos(epsilon) * sin(lambda); sin(epsilon) * sin(lambda)];
    sunHat = sunHat / norm(sunHat);
end

function tf = isInCylindricalEclipse(rEciMeters, sunHat)
    earthRadiusMeters = 6378.137e3;
    behindEarth = dot(rEciMeters, sunHat) < 0;
    distanceFromSunLine = norm(cross(rEciMeters, sunHat));
    tf = behindEarth && distanceFromSunLine < earthRadiusMeters;
end

function directionEci = groundStationDirectionEci(latDeg, lonDeg, altMeters, sat, t)
    satEci = states(sat, t, "CoordinateFrame", "inertial");
    satEci = satEci(:);
    gsEci = geodeticToApproxEci(latDeg, lonDeg, altMeters, t);
    directionEci = gsEci - satEci;
    directionEci = directionEci / norm(directionEci);
end

function rEci = geodeticToApproxEci(latDeg, lonDeg, altMeters, t)
    earthSemiMajorMeters = 6378137.0;
    flattening = 1 / 298.257223563;
    e2 = flattening * (2 - flattening);

    lat = deg2rad(latDeg);
    lon = deg2rad(lonDeg);
    N = earthSemiMajorMeters / sqrt(1 - e2 * sin(lat)^2);

    rEcef = [
        (N + altMeters) * cos(lat) * cos(lon);
        (N + altMeters) * cos(lat) * sin(lon);
        (N * (1 - e2) + altMeters) * sin(lat)];

    theta = gmstRadians(t);
    R_eci_from_ecef = [
        cos(theta), -sin(theta), 0;
        sin(theta),  cos(theta), 0;
        0,           0,          1];
    rEci = R_eci_from_ecef * rEcef;
end

function theta = gmstRadians(t)
    jd = juliandate(t);
    centuries = (jd - 2451545.0) / 36525;
    thetaDeg = 280.46061837 + 360.98564736629 * (jd - 2451545.0) ...
        + 0.000387933 * centuries^2 - centuries^3 / 38710000;
    theta = deg2rad(mod(thetaDeg, 360));
end

function y = rmsFinite(x)
    x = x(isfinite(x));
    if isempty(x)
        y = nan;
    else
        y = sqrt(mean(x.^2));
    end
end

function v = metricValue(summary, metricName)
    idx = summary.Metric == string(metricName);
    if any(idx)
        v = summary.Value(find(idx, 1, "first"));
    else
        v = nan;
    end
end

function writeMarkdownReport(reportFile, summary, figureOutDir, adcsInput, tleInput)
    lines = strings(0, 1);
    lines(end+1) = "# RIoT-2 Enhanced ADCS and Access Analysis";
    lines(end+1) = "";
    lines(end+1) = "## Inputs";
    lines(end+1) = "- ADCS log: `" + string(adcsInput) + "`";
    lines(end+1) = "- TLE name: `" + string(tleInput(1)) + "`";
    lines(end+1) = "- TLE line 1: `" + string(tleInput(2)) + "`";
    lines(end+1) = "- TLE line 2: `" + string(tleInput(3)) + "`";
    lines(end+1) = "";
    lines(end+1) = "## Key Results";
    lines(end+1) = sprintf("- Rows analyzed: %.0f", metricValue(summary, "rows"));
    lines(end+1) = sprintf("- Quaternion coverage: %.2f %%", metricValue(summary, "quaternion_coverage_percent"));
    lines(end+1) = sprintf("- Max raw quaternion norm error: %.6f", metricValue(summary, "quaternion_norm_max_abs_error"));
    lines(end+1) = sprintf("- Body-rate RMS: %.4f deg/s", metricValue(summary, "body_rate_rms_degps"));
    lines(end+1) = sprintf("- Stable-rate fraction: %.2f %%", metricValue(summary, "stable_rate_fraction_percent"));
    lines(end+1) = sprintf("- Median nadir pointing error: %.2f deg", metricValue(summary, "median_nadir_error_deg"));
    lines(end+1) = sprintf("- Median Sun pointing error: %.2f deg", metricValue(summary, "median_sun_error_deg"));
    lines(end+1) = sprintf("- Geometric access windows: %.0f", metricValue(summary, "geometric_access_windows"));
    lines(end+1) = sprintf("- Sensor-limited access windows: %.0f", metricValue(summary, "sensor_limited_access_windows"));
    lines(end+1) = "";
    lines(end+1) = "## Figures";
    figureFiles = dir(fullfile(figureOutDir, "*.png"));
    for i = 1:numel(figureFiles)
        lines(end+1) = "- `" + string(figureFiles(i).name) + "`";
    end
    lines(end+1) = "";
    lines(end+1) = "## Interpretation Notes";
    lines(end+1) = "- Quaternion norm should remain close to 1 before normalization. Large deviations indicate invalid attitude telemetry.";
    lines(end+1) = "- Body-rate RMS and maximum rate indicate whether the satellite attitude is stable or still slewing.";
    lines(end+1) = "- Pointing errors depend on assumed body axes and should be calibrated with spacecraft CAD or ADCS documentation.";
    writelines(lines, reportFile);
end
