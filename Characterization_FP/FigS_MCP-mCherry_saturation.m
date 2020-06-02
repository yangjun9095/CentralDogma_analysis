clear
close all
addpath('../utilities');

type_name = 'MCPmCherry';
% specify quantiles to use for analysis
quantiles_in = [15:15:90 98]/100;
tlims_in = [5 25];
nc14_flag_in = true;
% define basic ID variables

% define the file paths
DropboxFolder =  'S:\YangJoon\Dropbox';
MainFigPath = 'S:\YangJoon\Dropbox\Garcia Lab\Figures\CentralDogmaFigures';

FigPath = [MainFigPath, filesep, 'MCPmCherry_Saturation']
mkdir(FigPath);

% Use the getProjectPrefixes to read out the Prefix names from the DataStatus.xlsx tab
prefix_cell = {};
prefix_cell = getProjectPrefixes('MCPmCherry');


%% Part1. Nick and Meghan's script (using quantiles, across all AP axis)
% intialize data arrays
spot_mean_array = NaN(numel(quantiles_in),numel(prefix_cell));
spot_se_array = NaN(numel(quantiles_in),numel(prefix_cell));

offset_mean_array = NaN(numel(quantiles_in),numel(prefix_cell));
offset_se_array = NaN(numel(quantiles_in),numel(prefix_cell));

% call function to calculate offset and fluo quantiles for each prefix

disp('calculating quantiles...')
for p = 1:numel(prefix_cell)    
    Prefix = prefix_cell{p};
    disp(['analyzing ' Prefix '...'])
    [spot_mean_array(:,p), spot_se_array(:,p), offset_mean_array(:,p), offset_se_array(:,p)...
        , nc14_flag, quantiles, tlims] = ...
                    mcp_saturation_analysis(Prefix,'quantiles',quantiles_in,'nc14_flag',nc14_flag_in,'tlims',tlims_in);
end

% Make Figure
nrows = size(offset_mean_array,1);
inc = 1/nrows;
close all
pct_fig = figure;
hold on
hm_cm = flipud(brewermap(nrows,'Spectral'));
colormap(hm_cm);
plot(offset_mean_array,spot_mean_array,'-o','Color',[0 0 0 .3])
for i = 1:nrows
    scatter(offset_mean_array(i,:),spot_mean_array(i,:),'MarkerFaceColor',hm_cm(i,:),'MarkerEdgeColor','black')
    pause
end
xlabel('MCP offset (au)')
ylabel('spot fluorescence (au)')
h = colorbar('YTick',inc/2:inc:1,'YTickLabel',quantiles*100);
ylabel(h,'quantile (%)')
grid on
box on
set(gca,'Fontsize',14)
saveas(pct_fig,[FigPath,filesep, type_name '_prctile_scatters.tif'])
saveas(pct_fig,[FigPath,filesep, type_name '_prctile_scatters.pdf'])


% save data structure with analysis details
detail_struct.nc14_flag = nc14_flag;
detail_struct.tlims = tlims;
detail_struct.prefix_cell = prefix_cell;
detail_struct.quantiles = quantiles;

detail_struct.offset_mean_array = offset_mean_array;
detail_struct.offset_se_array = offset_se_array;
detail_struct.spot_mean_array = spot_mean_array;
detail_struct.spot_se_array = spot_se_array;

save([FigPath, filesep, type_name '_data.mat'],'detail_struct')

%% Part2. generate time-traces of the mean ms2 spot fluorescence 
% This is for the FigS1(A).
% plot the averaged spot fluo (over ON nuclei) over some time window
% of nc14 (0-30 min) on top of each other. We will try hard to synchronize
% the different datasets by synchronizing with the time point where it
% has the peak fluo.

% Let's use the LoadMS2Sets to extract useful info from each dataset, then save as a structure
% mean spot fluo, SD, SEM, Time, nc
%MCPdata = LoadMS2Sets('MCPmCherry');
% Note that some datasets are not listed in the MovieDatabase, perhaps in
% the old Tandem server. I need to make these options such that the
% LoadMS2Sets can be flexible in that option.

% define the parameters for plotting
% maximum time point in nc14
T_max = 20; % min
% pick one APbin, then plot the mean spot fluo on top of each other.
APbin = 11;
% I'll set the T_peak as 10min in here for an easier synchronization.
T_peak = 10; %min


% flag the genotypic dosage for coloring
% Note. This can be read out from the DataStatus.xlsx 
dosage_genotype_flag = {'2','2','2','3','3','3','4','4','4'};
dosage_gen_flag = [1,1,1,2,2,2,6,6,6];

% color map
hm_cm = flipud(brewermap(nrows,'Spectral'));
% colormap(hm_cm);

%% Loop over all dataset
traceFig = figure;
hold on
for setNum = 1:length(MCPdata)
    % initialize some fields
    clear fluo
    clear frameWindow
    
    % extract useful fields
    Time = MCPdata(setNum).ElapsedTime;
    nc14 = MCPdata(setNum).nc14;
    
    fluo_mean = MCPdata(setNum).MeanVectorAP;
    fluo_sd = MCPdata(setNum).SDVectorAP;
    num_particles = MCPdata(setNum).NParticlesAP;
    fluo_sem = fluo_sd./sqrt(num_particles);
    
    % calculate the frame window using the T_max
    tRes = median(diff(Time));
    num_frames = ceil(T_max/tRes);
    if nc14+num_frames <= length(Time)
        frameWindow = nc14:nc14+num_frames;
    else
        frameWindow = nc14:length(Time);
    end
    
    % 20sec resolution -> plot every other frame
    if tRes < 0.5
        frameWindow = frameWindow(1:2:end);
    end
    
    % find the peak time
    fluo = fluo_mean(frameWindow,APbin);
    
    % if not all the data points are NaNs
    if sum(isnan(fluo)) ~= length(frameWindow)
        index_time_peak = (fluo==max(fluo));
        time_trunc = Time(frameWindow);
        time_peak(setNum) = time_trunc(index_time_peak);
        time_sync(setNum) = T_peak - time_peak(setNum);
    
        % color
        colorflag = dosage_gen_flag(setNum);

        plothandle(setNum) = errorbar(Time(frameWindow) + time_sync(setNum), fluo_mean(frameWindow,APbin),...
                    fluo_sem(frameWindow,APbin), 'Color',hm_cm(colorflag,:))
    end
end

% edit plot format
xlim([0 25])
xlabel('time into nc14 (min)')
ylabel('MS2 spot fluorescence (AU)')


legend([plothandle(1),plothandle(5),plothandle(7)],['2 copies';'3 copies';'4 copies'])
StandardFigure(traceFig, traceFig.CurrentAxes)

% Save the plot
%saveas(traceFig, [FigPath,filesep, 'time_trace_MS2fluo_Dosages_',num2str((APbin-1)*2.5),'%.tif'])
%saveas(traceFig, [FigPath,filesep, 'time_trace_MS2fluo_Dosages_',num2str((APbin-1)*2.5),'%.pdf'])
%% Part3. Plotting the maximum fluo over MCP offset
maxfluo_offset_fig = figure;
hold on
for setNum = 1:length(MCPdata)-1 % Let's get rid of the last dataset as it does not contain the peak
    % initialize some fields
    clear fluo
    clear frameWindow
    
    % extract useful fields
    Time = MCPdata(setNum).ElapsedTime;
    nc14 = MCPdata(setNum).nc14;
    
    fluo_mean = MCPdata(setNum).MeanVectorAP;
    fluo_sd = MCPdata(setNum).SDVectorAP;
    num_particles = MCPdata(setNum).NParticlesAP;
    fluo_sem = fluo_sd./sqrt(num_particles);
    
    % offset of spots
    offset = MCPdata(setNum).MeanOffsetVector;
    
    % calculate the frame window using the T_max
    tRes = median(diff(Time));
    num_frames = ceil(T_max/tRes);
    if nc14+num_frames <= length(Time)
        frameWindow = nc14:nc14+num_frames;
    else
        frameWindow = nc14:length(Time);
    end
    
    % find the maximum fluo and sem corresponding to that frame
    Fluo = fluo_mean(frameWindow,APbin);
    Fluo_sem = fluo_sem(frameWindow,APbin);
    
    Fluo_max = max(Fluo);
    if ~isnan(Fluo_max)
        index_time_peak = (Fluo==Fluo_max);

        Fluo_max_sem = Fluo_sem(index_time_peak);

        % calculate the MCP offset during that frame window
        offset_mean = nanmean(offset(frameWindow));

        % color
        colorflag = dosage_gen_flag(setNum);
        plothandle(setNum) = errorbar(offset_mean, Fluo_max, Fluo_max_sem,'ok',...
                    'MarkerFaceColor',hm_cm(colorflag,:),'MarkerEdgeColor','black');
    end
end

% edit plot format
xlim([0 18])
ylim([0 1300])
yticks([0 200 400 600 800 1000])
xlabel('MCP offset (AU)')
ylabel('MS2 spot fluorescence (AU)')
% 
% 
legend([plothandle(2),plothandle(5),plothandle(7)],['2 copies';'3 copies';'4 copies'])
StandardFigure(maxfluo_offset_fig, maxfluo_offset_fig.CurrentAxes)
% 
% % Save the plot
saveas(maxfluo_offset_fig, [FigPath,filesep, 'max_MS2fluo_Offset_',num2str((APbin-1)*2.5),'%.tif'])
saveas(maxfluo_offset_fig, [FigPath,filesep, 'max_MS2fluo_Offset_',num2str((APbin-1)*2.5),'%.pdf'])