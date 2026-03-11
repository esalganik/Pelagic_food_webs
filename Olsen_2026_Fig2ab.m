%% Import of ALS and ROV data
clear; clc;
% ALS freeboard data from Hutter et al. (2023) doi:10.1594/PANGAEA.950896
project = 'awi-mosaic-l4-als-vq580-stere-0p50m-20200107T092304-20200107T105236-fv2p0.nc';
freeboard = ncread(project,'freeboard_estimate'); xc = ncread(project,'xc'); yc = ncread(project,'yc'); uncertainty = ncread(project,'freeboard_uncertainty');

[~,idx_x_min] = min(abs(xc--500)); [~,idx_x_max] = min(abs(xc--300)); % x-axis
[~,idx_y_min] = min(abs(yc--350)); [~,idx_y_max] = min(abs(yc--200)); % y-axis
y_als = xc(idx_x_min:idx_x_max); % selected area X
x_als = yc(idx_y_min:idx_y_max); % selected area Y
z_als = freeboard(idx_x_min:idx_x_max,idx_y_min:idx_y_max); % selected area Z
x_als = flipud(x_als);

% ROV ice draft data from Anhaus et al. (2024) doi:10.1594/PANGAEA.951077
project = "PS122_2_19_27_20200107_ROV_MULTIBEAM_v1_new.tab"; % Jan 07
data = tdfread(project);
X = data.Distance0x2C_realtive0x2C_X; Y = data.Distance0x2C_relative0x2C_Y; Z = data.Sea0x2Dice_draft_0x5Bm0x5D;

% making universal grid
for i = 1:length(project) % interpolation for chosen fixed mesh
    x = X; y = Y; z = Z;
    xv_rov = -350:0.5:350; % all X
    yv_rov = -350:0.5:350; % all Y
    [XX,YY] = meshgrid(xv_rov, yv_rov);
    ZZ_rov = griddata(x,y,z,XX,YY);
end
z_rov = ZZ_rov; x_rov = xv_rov; y_rov = yv_rov; z_rov = rot90(z_rov,2);
clearvars -except z_rov x_rov y_rov x_als y_als z_als

%% ALS / ROV colocation  —  rotation + translation
close all
% USER INPUTS
theta_init   = -35; % initial angle betwee ALD and ROV grid orientation
search_range =  20; % angle range
max_shift    =   5; % shift range

anchor_rov = [129,  -132]; % anchor point at ROV coordinates
anchor_als = [-266, -396]; % anchor point at ALS coordinates

x_rov_roi  = [40,  180]; % crop
y_rov_roi  = [-190, -70]; % crop

% CROP ROV TO ROI ONCE
xm = x_rov >= x_rov_roi(1) & x_rov <= x_rov_roi(2);
ym = y_rov >= y_rov_roi(1) & y_rov <= y_rov_roi(2);
[Xroi, Yroi] = meshgrid(x_rov(xm), y_rov(ym));
Zroi     = z_rov(ym, xm);
Zroi_vec = Zroi(:);

% grid resolution for correlation (use ALS resolution)
gres = abs(x_als(2) - x_als(1));

% Helper: rotate ROV, interpolate onto common grid, return Pearson r
% dxt/dyt = additional translation on top of anchor (default 0)

    function [r, xa, ya, Zrov_sub, Zals_sub] = eval_corr(Xroi, Yroi, Zroi_vec, ...
            theta, dxt, dyt, anchor_rov, anchor_als, x_als, y_als, z_als)
        ct = cosd(theta); st = sind(theta);
        Xt = ct*(Xroi-anchor_rov(1)) - st*(Yroi-anchor_rov(2)) + anchor_als(1) + dxt;
        Yt = st*(Xroi-anchor_rov(1)) + ct*(Yroi-anchor_rov(2)) + anchor_als(2) + dyt;
        xa = x_als(x_als >= min(Xt(:)) & x_als <= max(Xt(:)));
        ya = y_als(y_als >= min(Yt(:)) & y_als <= max(Yt(:)));
        if isempty(xa) || isempty(ya), r=NaN; Zrov_sub=[]; Zals_sub=[]; return; end
        [Xa,Ya] = meshgrid(xa, ya);
        Zals_sub = interp2(x_als', y_als, z_als, Xa, Ya, 'linear', NaN);
        F = scatteredInterpolant(Xt(:), Yt(:), Zroi_vec, 'linear', 'none');
        Zrov_sub = F(Xa, Ya);
        valid = ~isnan(Zals_sub) & ~isnan(Zrov_sub);
        if sum(valid(:)) < 100, r=NaN; return; end
        a = Zals_sub(valid); b = Zrov_sub(valid);
        r = (mean(a.*b) - mean(a)*mean(b)) / (std(a)*std(b));
    end

% STEP 1: COARSE ANGLE SEARCH
coarse_angles = (theta_init - search_range) : 1 : (theta_init + search_range);
r_coarse = nan(size(coarse_angles));
fprintf('Step 1/3 — coarse angle search (%d angles)...\n', numel(coarse_angles));
for k = 1:numel(coarse_angles)
    r_coarse(k) = eval_corr(Xroi, Yroi, Zroi_vec, coarse_angles(k), 0, 0, ...
        anchor_rov, anchor_als, x_als, y_als, z_als);
end
[~, ci] = max(r_coarse);
theta_coarse = coarse_angles(ci);
fprintf('  Coarse best: %.1f°  (r=%.4f)\n', theta_coarse, r_coarse(ci));

% STEP 2: FINE ANGLE SEARCH
fine_angles = (theta_coarse - 1) : 0.1 : (theta_coarse + 1);
r_fine = nan(size(fine_angles));
fprintf('Step 2/3 — fine angle search (%d angles)...\n', numel(fine_angles));
for k = 1:numel(fine_angles)
    r_fine(k) = eval_corr(Xroi, Yroi, Zroi_vec, fine_angles(k), 0, 0, ...
        anchor_rov, anchor_als, x_als, y_als, z_als);
end
[~, fi] = max(r_fine);
theta_opt = fine_angles(fi);
fprintf('  Fine best:   %.2f°  (r=%.4f)\n', theta_opt, r_fine(fi));

% STEP 3: TRANSLATION via normxcorr2  (FFT-based — single interpolation)
fprintf('Step 3/3 — translation search via normxcorr2...\n');

% Interpolate rotated ROV and ALS onto a common grid ONCE
[~, xa, ya, Zrov_sub, Zals_sub] = eval_corr(Xroi, Yroi, Zroi_vec, theta_opt, 0, 0, ...
    anchor_rov, anchor_als, x_als, y_als, z_als);

% Replace NaNs with mean (normxcorr2 needs no NaNs)
Zrov_c = Zrov_sub; Zrov_c(isnan(Zrov_c)) = mean(Zrov_sub(:),'omitnan');
Zals_c = Zals_sub; Zals_c(isnan(Zals_c)) = mean(Zals_sub(:),'omitnan');

% normxcorr2: template = ROV, image = ALS
C = normxcorr2(Zrov_c, Zals_c);

% Find peak within max_shift window (in pixels)
max_pix = round(max_shift / gres);
[nr, nc] = size(C);
cr = round(nr/2); cc = round(nc/2);  % zero-shift location in xcorr output
row_win = max(1, cr-max_pix) : min(nr, cr+max_pix);
col_win = max(1, cc-max_pix) : min(nc, cc+max_pix);
C_win = C(row_win, col_win);

[~, pidx] = max(C_win(:));
[pr, pc]  = ind2sub(size(C_win), pidx);

% Convert peak position to pixel shift, then to metres
dy_pix = pr - (max_pix+1);
dx_pix = pc - (max_pix+1);
dx_opt = dx_pix * gres;
dy_opt = dy_pix * gres;

fprintf('  Best shift:  dx=%.1f m  dy=%.1f m\n', dx_opt, dy_opt);

% Recompute Pearson r at optimal translation for reporting
r_final = eval_corr(Xroi, Yroi, Zroi_vec, theta_opt, dx_opt, dy_opt, ...
    anchor_rov, anchor_als, x_als, y_als, z_als);
fprintf('\n=== RESULT: theta=%.2f°  dx=%.1f m  dy=%.1f m  r=%.4f ===\n', ...
    theta_opt, dx_opt, dy_opt, r_final);

% APPLY OPTIMAL TRANSFORM FOR PLOTTING
ct = cosd(theta_opt); st = sind(theta_opt);
Xt = ct*(Xroi-anchor_rov(1)) - st*(Yroi-anchor_rov(2)) + anchor_als(1) + dx_opt;
Yt = st*(Xroi-anchor_rov(1)) + ct*(Yroi-anchor_rov(2)) + anchor_als(2) + dy_opt;
xa_f = x_als(x_als >= min(Xt(:)) & x_als <= max(Xt(:)));
ya_f = y_als(y_als >= min(Yt(:)) & y_als <= max(Yt(:)));
[Xa_f, Ya_f] = meshgrid(xa_f, ya_f);
F = scatteredInterpolant(Xt(:), Yt(:), Zroi_vec, 'linear', 'none');
Zrov_als = F(Xa_f, Ya_f);

% =========================================================================
% PLOTS
% =========================================================================
load('broc.mat');
range = -0.05:0.05:2.5;

% Figure 1 — angle search
figure('Units','inches','Position',[1 1 10 4],'Color','w');
subplot(1,2,1);
plot(coarse_angles, r_coarse,'b-o','LineWidth',1.2,'MarkerSize',5); hold on;
xline(theta_opt,'r--','LineWidth',1.5);
xlabel('Angle (deg)'); ylabel('Pearson r'); title('Coarse angle search');
grid on; set(gca,'FontSize',8);
subplot(1,2,2);
plot(fine_angles, r_fine,'b-o','LineWidth',1.2,'MarkerSize',5); hold on;
xline(theta_opt,'r--',sprintf(' %.2f°',theta_opt),'LineWidth',1.5,'FontSize',8);
xlabel('Angle (deg)'); ylabel('Pearson r'); title('Fine angle search');
grid on; set(gca,'FontSize',8);

% Figure 2 — side-by-side ALS and ROV on same grid and colorscale
% Crop ALS to the same sub-region as the colocated ROV
Zals_crop = interp2(x_als', y_als, z_als, Xa_f, Ya_f, 'linear', NaN);

figure('Units','inches','Position',[2 2 11 5],'Color','w');
t = tiledlayout(1, 2, 'TileSpacing','compact','Padding','compact');

nexttile;
imagesc(xa_f, ya_f, Zals_crop); hold on;
set(gca,'YDir','normal'); axis image;
clim([0 2.5]); colormap(broc);
xlabel('x  (ALS)','FontSize',9); ylabel('y  (ALS)','FontSize',9);
title('ALS','FontSize',10,'FontWeight','normal');
set(gca,'FontSize',8,'Box','on');

nexttile;
imagesc(xa_f, ya_f, Zrov_als);
set(gca,'YDir','normal'); axis image;
clim([0 6]); colormap(broc);
xlabel('x  (ALS)','FontSize',9);
title('ROV  (colocated)','FontSize',10,'FontWeight','normal');
set(gca,'FontSize',8,'Box','on');

cb = colorbar; cb.Layout.Tile = 'east';
ylabel(cb,'freeboard (m)','FontSize',9);
title(t, sprintf('\\theta=%.2f°   dx=%.1f m   dy=%.1f m   r=%.4f', ...
    theta_opt, dx_opt, dy_opt, r_final),'FontSize',10);

clearvars -except z_rov x_rov y_rov x_als y_als z_als xa_f ya_f Zrov_als Zals_crop

%% Final plot 2ab
close all
figure
tile = tiledlayout(1,2); tile.TileSpacing = 'compact'; tile.Padding = 'none';
c{1} = [0.0000 0.4470 0.7410]; c{2} = [0.8500 0.3250 0.0980]; c{3} = [0.9290 0.6940 0.1250]; c{4} = [0.4940 0.1840 0.5560];	c{5} = [0.4660 0.6740 0.1880]; c{6} = [0.3010 0.7450 0.9330]; c{7} = [0.6350 0.0780 0.1840]; % colors
fsz = 10; fsz_s = 9;
nexttile
dx = 304; dy = 431;
range = -0.05:0.05:2.5; % Graph accuracy
% contourf(x_als+dx,y_als+dy,z_als,range,'LabelSpacing',10,'edgecolor','none'); hold on;
contourf(xa_f+dx,ya_f+dy,Zals_crop,range,'LabelSpacing',10,'edgecolor','none'); hold on;
p = plot(36.5,34.5,'o','color','k'); p.MarkerSize = 4.0; set(p,'markerfacecolor',c{3}); % Ridge coring 10 Jan
p = plot(2.5,42,'o','color','k'); p.MarkerSize = 4.0; set(p,'markerfacecolor',c{5}); % FYI coring 3 Jan
p = plot(44.5,63.5,'o','color','k'); p.MarkerSize = 4.0; set(p,'markerfacecolor',c{5}); % FYI coring 3 Jan
p = plot(15,48,'o','color','k'); p.MarkerSize = 4.0; set(p,'markerfacecolor',c{2}); % DTC25
p = plot(20,50.5,'o','color','k'); p.MarkerSize = 4.0; set(p,'markerfacecolor',c{2}); % T60
p = plot(24,52.5,'o','color','k'); p.MarkerSize = 4.0; set(p,'markerfacecolor',c{2}); % DTC24
p = text(16,44,'DTC25'); set(p,'Color','k','HorizontalAlignment','left','FontSize',fsz_s);
p = text(21,48.5,'T60'); set(p,'Color','k','HorizontalAlignment','left','FontSize',fsz_s);
p = text(25,52.5,'DTC24'); set(p,'Color','k','HorizontalAlignment','left','FontSize',fsz_s);
p = text(38.5,34.5,'Ridge coring'); set(p,'Color','k','HorizontalAlignment','left','FontSize',fsz);
p = text(2.5,38,'Level ice'); set(p,'Color','w','HorizontalAlignment','left','FontSize',fsz);
p = text(2.5,33,'coring'); set(p,'Color','w','HorizontalAlignment','left','FontSize',fsz);
p = text(47,68.5,'Level ice'); set(p,'Color','w','HorizontalAlignment','left','FontSize',fsz);
p = text(47,63.5,'coring'); set(p,'Color','w','HorizontalAlignment','left','FontSize',fsz);
load('broc.mat'); colormap((broc));
leg = legend('','Ridge coring','Level ice coring','','IMBs','box','off','NumColumns',1);
set(leg,'textcolor','w','FontSize',fsz_s,'Location','southwest'); leg.ItemTokenSize = [30*0.4,18*0.4];
leg.BoxFace.ColorType='truecoloralpha'; leg.BoxFace.ColorData=uint8(255*[1 1 1 0.7]'); % transparent legend
hXLabel = xlabel('x (m)'); hYLabel = ylabel('y (m)'); set([hXLabel hYLabel gca],'FontSize',fsz,'FontWeight','normal');
colorbar; hBar1 = colorbar; ylabel(hBar1,'Snow freeboard (m)','FontSize',fsz); clim([0 2.5]);
xlim([0 70]); xticks(0:10:70); ylim([0 80]);

nexttile
range = -1:0.2:9; % Graph accuracy
contourf(xa_f+dx, ya_f+dy, Zrov_als,range,'-','ShowText','on','LabelSpacing',400,'LineColor','none'); hold on
colorbar; load('broc.mat'); colormap(broc); hBar1 = colorbar; ylabel(hBar1,'Ice draft (m)','FontSize',8);
p = plot(36.5,34.5,'o','color','k'); p.MarkerSize = 4.0; set(p,'markerfacecolor',c{3}); % Ridge coring 10 Jan
p = plot(2.5,42,'o','color','k'); p.MarkerSize = 4.0; set(p,'markerfacecolor',c{5}); % FYI coring 3 Jan
p = plot(44.5,63.5,'o','color','k'); p.MarkerSize = 4.0; set(p,'markerfacecolor',c{5}); % FYI coring 3 Jan
p = plot(15,48,'o','color','k'); p.MarkerSize = 4.0; set(p,'markerfacecolor',c{2}); % DTC25
p = plot(20,50.5,'o','color','k'); p.MarkerSize = 4.0; set(p,'markerfacecolor',c{2}); % T60
p = plot(24,52.5,'o','color','k'); p.MarkerSize = 4.0; set(p,'markerfacecolor',c{2}); % DTC24
p = text(16,44,'DTC25'); set(p,'Color','w','HorizontalAlignment','left','FontSize',fsz_s);
p = text(21,48.5,'T60'); set(p,'Color','k','HorizontalAlignment','left','FontSize',fsz_s);
p = text(25,52.5,'DTC24'); set(p,'Color','k','HorizontalAlignment','left','FontSize',fsz_s);
p = text(38.5,34.5,'Ridge coring'); set(p,'Color','k','HorizontalAlignment','left','FontSize',fsz);
p = text(2.5,38,'Level ice'); set(p,'Color','w','HorizontalAlignment','left','FontSize',fsz);
p = text(2.5,33,'coring'); set(p,'Color','w','HorizontalAlignment','left','FontSize',fsz);
p = text(47,68.5,'Level ice'); set(p,'Color','w','HorizontalAlignment','left','FontSize',fsz);
p = text(47,63.5,'coring'); set(p,'Color','w','HorizontalAlignment','left','FontSize',fsz);
leg = legend('','Ridge coring','Level ice coring','','IMBs','box','off','NumColumns',1);
set(leg,'textcolor','w','FontSize',fsz_s,'Location','southwest'); leg.ItemTokenSize = [30*0.4,18*0.4];
leg.BoxFace.ColorType='truecoloralpha'; leg.BoxFace.ColorData=uint8(255*[1 1 1 0.7]'); % transparent legend
hXLabel = xlabel('x (m)'); hYLabel = ylabel('y (m)'); set([hXLabel hYLabel gca],'FontSize',fsz,'FontWeight','normal');
colorbar; hBar1 = colorbar; ylabel(hBar1,'Ice draft (m)','FontSize',fsz); clim([1 6]);
xlim([0 70]); xticks(0:10:70); ylim([0 80]);

annotation('textbox',[0 .51 0.02 .51],'String','(a)','EdgeColor','none','HorizontalAlignment','center','FontSize',fsz);
annotation('textbox',[0.35 .51 0.35 .51],'String','(b)','EdgeColor','none','HorizontalAlignment','center','FontSize',fsz);
set(gcf,'Units','inches','Position',[3 4 6.2 2.2])
clearvars -except x_als y_als z_als xa_f ya_f Zals_crop Zrov_als

%% NetCDF export
ncFile = 'colocated_freeboard_draft_20200107.nc';

if isfile(ncFile)
    delete(ncFile)
end

Nx = length(xa_f);   % 301
Ny = length(ya_f);   % 344

% Sanity checks
assert(isequal(size(Zals_crop), [Ny, Nx]), ...
    'Zals_crop must have size [length(ya_f), length(xa_f)]')
assert(isequal(size(Zrov_als), [Ny, Nx]), ...
    'Zrov_als must have size [length(ya_f), length(xa_f)]')

% Create coordinate variables
nccreate(ncFile, 'xa_f', ...
    'Dimensions', {'x', Nx}, ...
    'Datatype', 'double');

nccreate(ncFile, 'ya_f', ...
    'Dimensions', {'y', Ny}, ...
    'Datatype', 'double');

% Create 2D data variables
nccreate(ncFile, 'Zals_crop', ...
    'Dimensions', {'y', Ny, 'x', Nx}, ...
    'Datatype', 'double');

nccreate(ncFile, 'Zrov_als', ...
    'Dimensions', {'y', Ny, 'x', Nx}, ...
    'Datatype', 'double');

% Write coordinate data
ncwrite(ncFile, 'xa_f', xa_f);
ncwrite(ncFile, 'ya_f', ya_f);

% Write 2D fields exactly as matrices
ncwrite(ncFile, 'Zals_crop', Zals_crop);
ncwrite(ncFile, 'Zrov_als', Zrov_als);

% Coordinate variable attributes
ncwriteatt(ncFile, 'xa_f', 'long_name', 'horizontal coordinate x');
ncwriteatt(ncFile, 'xa_f', 'units', 'm');

ncwriteatt(ncFile, 'ya_f', 'long_name', 'vertical coordinate y');
ncwriteatt(ncFile, 'ya_f', 'units', 'm');

% Data variable attributes
ncwriteatt(ncFile, 'Zals_crop', 'long_name', 'snow freeboard');
ncwriteatt(ncFile, 'Zals_crop', 'units', 'm');
ncwriteatt(ncFile, 'Zals_crop', 'coordinates', 'xa_f ya_f');

ncwriteatt(ncFile, 'Zrov_als', 'long_name', 'ice draft');
ncwriteatt(ncFile, 'Zrov_als', 'units', 'm');
ncwriteatt(ncFile, 'Zrov_als', 'coordinates', 'xa_f ya_f');

% Global attributes
ncwriteatt(ncFile,"/","title","Colocated snow freeboard and ice draft for 7 Jan 2020 for MOSAiC Expedition");
ncwriteatt(ncFile,"/","Conventions","CF-1.7");
ncwriteatt(ncFile,"/","contributor_name", ...
    "Evgenii Salganik; Philipp Anhaus; Nils Hutter; Stefan Hendricks; Arttu Jutila; Gerit Birnbaum; Luisa von Albedyll; Robert Ricker; Christian Haas; Christian Katlein; Ilkka Matero; Marcel Nicolaus; Stefanie Arndt; Daniela Krampe; Benjamin Allen Lange; Julia Regnery; Jan Rohde; Martin Schiller; Daniela Ransby");
ncwriteatt(ncFile,"/","contributor_email","evgenii.salganik@awi.de");
ncwriteatt(ncFile,"/","institution","Alfred Wegener Institute for Polar and Marine Research");
ncwriteatt(ncFile,"/","creator_name","Evgenii Salganik");
ncwriteatt(ncFile,"/","creator_email","evgenii.salganik@awi.de");
ncwriteatt(ncFile,"/","project","Arctic PASSION");
ncwriteatt(ncFile,"/","summary","Colocated snow freeboard and ice draft gridded data for 7 Jan 2020 collected during the MOSAiC Expedition");
ncwriteatt(ncFile,"/","license","CC-0");

ncdisp(ncFile)