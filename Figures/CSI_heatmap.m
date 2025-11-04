% Define file paths
excel_file = 'D:\V1_Cortex_Paper\Cortex\CSI_excel.xlsx';
save_folder = 'D:\V1_Cortex_Paper\MultichannelDataTanks\CSI_heatmap';

% Load table
data = readtable(excel_file);

% Check required columns
required_cols = {'CSI', 'Auditory', 'Significance'};
missing_cols = setdiff(required_cols, data.Properties.VariableNames);
if ~isempty(missing_cols)
    error(['Missing columns: ', strjoin(missing_cols, ', ')]);
end

% Extract CSI values
CSI_all = data.CSI;
CSI_auditory = data.CSI(data.Auditory == 1);
CSI_significant = data.CSI(data.Significance == 1);

% Create custom colormap (white to indigo)
indigo = [153, 0, 51]/255;
custom_cmap = [linspace(1,indigo(1),256)', linspace(1,indigo(2),256)', linspace(1,indigo(3),256)'];

% Plotting helper function with fixed color scale
function plot_and_save_strip(csi_data, title_str, filename_prefix, cmap, folder)
    figure('Visible', 'off');
    imagesc(csi_data);
    colormap(cmap);
    caxis([0 1]);  % <-- Fixed color scale from 0 to 1
    colorbar;
    axis off;
    title(title_str, 'Interpreter', 'none');

    % File paths
    pdf_file = fullfile(folder, [filename_prefix, '.pdf']);
    tif_file = fullfile(folder, [filename_prefix, '.tif']);

    % Save
    exportgraphics(gcf, pdf_file, 'ContentType', 'vector');
    exportgraphics(gcf, tif_file, 'Resolution', 600);
    close;
end

% Plot and save all three versions
plot_and_save_strip(CSI_all, 'CSI - All', 'CSI_all', custom_cmap, save_folder);
plot_and_save_strip(CSI_auditory, 'CSI - Auditory = 1', 'CSI_auditory', custom_cmap, save_folder);
plot_and_save_strip(CSI_significant, 'CSI - Significance = 1', 'CSI_significance', custom_cmap, save_folder);


%% % Ensure save folder exists
if ~exist(save_folder, 'dir'); mkdir(save_folder); end

% ===== Save the vectors used for the heat maps to Excel =====
out_xlsx = fullfile(save_folder, 'CSI_heatmap_vectors.xlsx');

% If you want a fresh file each run, uncomment the next two lines:
% if exist(out_xlsx, 'file'); delete(out_xlsx); end

% Build tables (keep row indices for filtered subsets)
T_all = table((1:numel(CSI_all))', CSI_all(:), ...
    'VariableNames', {'RowIndex', 'CSI'});

aud_idx = find(data.Auditory == 1);
T_aud = table(aud_idx(:), CSI_auditory(:), ...
    'VariableNames', {'RowIndex', 'CSI'});

sig_idx = find(data.Significance == 1);
T_sig = table(sig_idx(:), CSI_significant(:), ...
    'VariableNames', {'RowIndex', 'CSI'});

% Write each set to its own sheet
writetable(T_all, out_xlsx, 'Sheet', 'CSI_all', 'FileType', 'spreadsheet');
writetable(T_aud, out_xlsx, 'Sheet', 'CSI_auditory', 'FileType', 'spreadsheet');
writetable(T_sig, out_xlsx, 'Sheet', 'CSI_significant', 'FileType', 'spreadsheet');

% Also save the colormap and color scale for reproducibility
CmapTbl = array2table(custom_cmap, 'VariableNames', {'R','G','B'});
writetable(CmapTbl, out_xlsx, 'Sheet', 'Colormap', 'FileType', 'spreadsheet');

MetaTbl = table(0, 1, "white→indigo (153,0,51)/255", ...
    'VariableNames', {'caxis_min','caxis_max','cmap_note'});
writetable(MetaTbl, out_xlsx, 'Sheet', 'Metadata', 'FileType', 'spreadsheet');

% ===== Plot & save strips (unchanged) =====
plot_and_save_strip(CSI_all, 'CSI - All', 'CSI_all', custom_cmap, save_folder);
plot_and_save_strip(CSI_auditory, 'CSI - Auditory = 1', 'CSI_auditory', custom_cmap, save_folder);
plot_and_save_strip(CSI_significant, 'CSI - Significance = 1', 'CSI_significance', custom_cmap, save_folder);
