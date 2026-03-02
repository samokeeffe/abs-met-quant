% load metabolites and Zipfs fitted concentration data 
load('zipfs_fit_10000.mat')
load('Recon3DModel_301.mat') %https://www.vmh.life/#downloadview
%% Pull original data and indices for each 

% Helper to define “has a numeric value” for a given column
hasValue = @(col) ~isnan(str2double(string(concentration_data{:,col})));

%keggCol = 15;

S = struct();

% CD4 (col 7)
S.CD4.mask = hasValue(7);
S.CD4.met  = concentration_data{S.CD4.mask, 1};
S.CD4.KEGG = concentration_data{S.CD4.mask, keggCol};

% CD8 (col 8)
S.CD8.mask = hasValue(8);
S.CD8.met  = concentration_data{S.CD8.mask, 1};
S.CD8.KEGG = concentration_data{S.CD8.mask, keggCol};

% Jurkat (col 5)
S.Jurkat.mask = hasValue(5);
S.Jurkat.met  = concentration_data{S.Jurkat.mask, 1};
S.Jurkat.KEGG = concentration_data{S.Jurkat.mask, keggCol};

%%
%[~, order] = sort(str2double(string(CD4_ordered_max(:,3))));
CD4_mets_ordered = string(S.CD4.met(CD4_ordered_max(:,3)));
CD8_mets_ordered = string(S.CD8.met(CD8_ordered_max(:,3)));
Jurkat_mets_ordered = string(S.Jurkat.met(Jurkat_ordered_max(:,3)));

%%
fig = figure;


% Pre-create axes handles for each tile (in order)
ax = gobjects(1,3);
tl = tiledlayout(fig, 1, 3, 'TileSpacing','compact', 'Padding','compact');


ax(1,1) = nexttile(tl, 1);
ax(1,2) = nexttile(tl, 2);
ax(1,3) = nexttile(tl, 3);


PlotHistograms(S.CD4, CD4_ordered_max, Recon3DModel ,ax(1,1),35,50);

PlotHistograms(S.CD8, CD8_ordered_max, Recon3DModel,ax(1,2), 32,50);
% % 
PlotHistograms(S.Jurkat, Jurkat_ordered_max, Recon3DModel, ax(1,3), 26, 50);


%end

 %% find metabolite indices
function idx = kegg2model(queryKEGG, modelKEGG)
    idx = arrayfun(@(k) find(modelKEGG == queryKEGG(k)), 1:numel(queryKEGG), 'UniformOutput', false);
end


function PlotHistograms(CellStructure, ordered_max, Recon3DModel,ax1,N,BinSize)


    CellStructure.ModelIndices = kegg2model(string(CellStructure.KEGG), string(Recon3DModel.metKEGGID));
    
    Smat = Recon3DModel.S;
    CellStructure.RxnIndices = cell(numel(CellStructure.ModelIndices),1);
    
    for k = 1:numel(CellStructure.ModelIndices)
    
        mets = CellStructure.ModelIndices{k};   % model metabolite indices
        if isempty(mets)
            CellStructure.RxnIndices{k} = [];
            continue
        end
    
        % Find reactions for all mapped metabolites
        rxnMask = any(Smat(mets,:) ~= 0, 1);
        CellStructure.RxnIndices{k} = find(rxnMask);
    end
    
    CellStructure.RxnCount = cell(numel(CellStructure.ModelIndices),1);
    
    for k = 1:numel(CellStructure.ModelIndices)
        CellStructure.RxnCount{k} = size(CellStructure.RxnIndices{k},2);
    end
    
    CellStructure.Ranks = cell(numel(CellStructure.ModelIndices),1);
    for k = 1:numel(CellStructure.KEGG)
        CellStructure.Ranks{k} = find(ordered_max(:,3)== k);
    end

    %% sort everything by metabolite concentration rank
    ranks = str2double(string(CellStructure.Ranks));
    [~, order] = sort(ranks);
    
    
    rxnCount = str2double(string(CellStructure.RxnCount(order)));
    metNames = CellStructure.met(order);
    ranks    = str2double(string(CellStructure.Ranks(order)));
    concentration = ordered_max(:,2);
    
    mdl = fitlm(concentration, rxnCount);
    r2 = mdl.Rsquared.Ordinary;
    

    
    %% version with cofactors removed
    
    %% ---- REMOVE SPECIFIC METABOLITES BY NAME ----
    removeMetNames = string({ ...
        %'atp', 'coa', 'gtp','utp','gthrd','gthox','adp','nad','udp','ctp','amp','gdp','dgmp','nadph','gmp','nadh','imp','nadp','camp','fad' ...   % <-- put the exact names you want removed here
    });
    
    % keep = NOT in removeMetNames
    keepMask = ~ismember(metNames, removeMetNames);
    %keepMask = ismember(metNames, removeMetNames);
    
    % Apply to all aligned vectors
    metNames_f = metNames(keepMask);
    rxnCount_f = rxnCount(keepMask);
    ranks_f    = ranks(keepMask);
    concentration_f = concentration(keepMask);
    
    mdl_f = fitlm(concentration_f, rxnCount_f);
    r2_f = mdl_f.Rsquared.Ordinary;
   
    
    hold(ax1, "off")         % reset hold for THIS tile
    histogram(ax1, rxnCount_f(1:N), 'facealpha',  0.5, 'BinWidth',BinSize);
    hold(ax1, "on")
    rxnCount_Zipf = rxnCount_f(1:N);
    rxnCount_NonZipf = rxnCount_f(N+1:end);
    log_rxnCount_zipf = log10(rxnCount_Zipf(find(rxnCount_Zipf>0)));
    log_rxnCount_nonzipfs= log10(rxnCount_NonZipf(find(rxnCount_NonZipf>0)));


    
    
    histogram(ax1, rxnCount_f(N+1:end), 'facealpha',  0.5, 'BinWidth',BinSize);
    
    hold(ax1,"off")
    
    ylabel(ax1, 'Count');
   
    xlabel(ax1, 'Reaction Count')
   
    legend(ax1,'Zipfian','Non-Zipfian')
  
    [h,p,ksteststat] = kstest2(log_rxnCount_zipf,log_rxnCount_nonzipfs);
    txt = sprintf('p = %.4f', p);
    xpos = 400;
    ypos = 10;
    text(ax1,xpos, ypos, txt, 'FontSize', 12, 'FontWeight','bold','Color','k');
  
    
 
end

%%
fig = figure;


% Pre-create axes handles for each tile (in order)
ax = gobjects(1,3);
tl = tiledlayout(fig, 1, 3, 'TileSpacing','compact', 'Padding','compact');


ax(1,1) = nexttile(tl, 1);
ax(1,2) = nexttile(tl, 2);
ax(1,3) = nexttile(tl, 3);


PlotCDF(S.CD4, CD4_ordered_max, Recon3DModel ,ax(1,1),35);

PlotCDF(S.CD8, CD8_ordered_max, Recon3DModel,ax(1,2), 32);
% % 
PlotCDF(S.Jurkat, Jurkat_ordered_max, Recon3DModel, ax(1,3), 26);

function PlotCDF(CellStructure, ordered_max, Recon3DModel,ax1,N)


    CellStructure.ModelIndices = kegg2model(string(CellStructure.KEGG), string(Recon3DModel.metKEGGID));
    
    Smat = Recon3DModel.S;
    CellStructure.RxnIndices = cell(numel(CellStructure.ModelIndices),1);
    
  % ---- preallocate before the loop (once) ----
        CellStructure.GeneIndices = cell(numel(CellStructure.ModelIndices),1);
        CellStructure.Genes       = cell(numel(CellStructure.ModelIndices),1);
        CellStructure.GeneCount   = cell(numel(CellStructure.ModelIndices),1);

    for k = 1:numel(CellStructure.ModelIndices)
    
        mets = CellStructure.ModelIndices{k};   % model metabolite indices
        if isempty(mets)
            CellStructure.RxnIndices{k} = [];
            continue
        end
    
        % Find reactions for all mapped metabolites
        rxnMask = any(Smat(mets,:) ~= 0, 1);
        CellStructure.RxnIndices{k} = find(rxnMask);
        rxns = CellStructure.RxnIndices{k};

    end
    CellStructure.RxnCount = cell(numel(CellStructure.ModelIndices),1);
    %%
    for k = 1:numel(CellStructure.ModelIndices)
        CellStructure.RxnCount{k} = size(CellStructure.RxnIndices{k},2);
    end
    %%
    CellStructure.Ranks = cell(numel(CellStructure.ModelIndices),1);
    for k = 1:numel(CellStructure.KEGG)
        CellStructure.Ranks{k} = find(ordered_max(:,3)== k);
    end
    
    %% sort everything by metabolite concentration rank
    ranks = str2double(string(CellStructure.Ranks));
    [~, order] = sort(ranks);
    
    
    rxnCount = str2double(string(CellStructure.RxnCount(order)));
    metNames = CellStructure.met(order);
    ranks    = str2double(string(CellStructure.Ranks(order)));
    concentration = ordered_max(:,2);
   
    %% ---- REMOVE SPECIFIC METABOLITES BY NAME ----
    removeMetNames = string({ ...
        %'atp', 'coa', 'gtp','utp','gthrd','gthox','adp','nad','udp','ctp','amp','gdp','dgmp','nadph','gmp','nadh','imp','nadp','camp','fad' ...   % <-- put the exact names you want removed here
    });
    
    % keep = NOT in removeMetNames
    keepMask = ~ismember(metNames, removeMetNames);
    %keepMask = ismember(metNames, removeMetNames);
    
    % Apply to all aligned vectors
    metNames_f = metNames(keepMask);
    rxnCount_f = rxnCount(keepMask);
    ranks_f    = ranks(keepMask);
    concentration_f = concentration(keepMask);
   
    
    hold(ax1, "off")         % reset hold for THIS tile
   
    hold(ax1, "on")
    rxnCount_Zipf = rxnCount_f(1:N);
    rxnCount_NonZipf = rxnCount_f(N+1:end);
    pos_rxnCount_Zipf = rxnCount_Zipf(rxnCount_Zipf>0);
    pos_rxnCount_NonZipf = rxnCount_NonZipf(rxnCount_NonZipf>0);
    log_rxnCount_zipf = log10(rxnCount_Zipf(find(rxnCount_Zipf>0)));
    log_rxnCount_nonzipfs= log10(rxnCount_NonZipf(find(rxnCount_NonZipf>0)));

    %[zipfs_cdf,x1] = ecdf(rxnCount_f(1:N));
    [zipfs_cdf,x1] = ecdf(pos_rxnCount_Zipf);
    %[nonzipfs_cdf,x2] = ecdf(rxnCount_f(N+1:end));
    [nonzipfs_cdf,x2] = ecdf(pos_rxnCount_NonZipf);
    zipfs_cdf = [0 ; zipfs_cdf];
    nonzipfs_cdf = [ 0; nonzipfs_cdf];
    x1 = [1; x1];
    x2 = [1; x2];
   
    hold(ax1,"off")
    stairs(ax1,x1, zipfs_cdf,'LineWidth',1.5)
 
    hold(ax1,"on")
    stairs(ax1,x2,nonzipfs_cdf,'LineWidth',1.5)
    box(ax1, "off");
    xscale(ax1,"log")
    ylabel(ax1,'Empirical CDF')
    xlabel(ax1, 'Reaction Count')
    %xscale("log")
    %xlim(ax1, [10-7 10-2])
    legend(ax1,'Zipfian','Non-Zipfian')
  
    [h,p,ksteststat] = kstest2(log_rxnCount_zipf,log_rxnCount_nonzipfs);
    txt = sprintf('p = %.4f', p);
    xpos = 2;
    ypos = 0.8;
    text(ax1,xpos, ypos, txt, 'FontSize', 12, 'FontWeight','bold','Color','k');
    axis(ax1, 'square')
    set(ax1,'TickDir','out')
 
  

    
end
%%
fig = figure;


% Pre-create axes handles for each tile (in order)
ax = gobjects(1,3);
tl = tiledlayout(fig, 1, 3, 'TileSpacing','compact', 'Padding','compact');


ax(1,1) = nexttile(tl, 1);
ax(1,2) = nexttile(tl, 2);
ax(1,3) = nexttile(tl, 3);

ScatterRxnCountFlip(S.CD4, CD4_ordered_max(:,2),CD4_ordered_max(:,3),Recon3DModel,ax(1,1))
ScatterRxnCountFlip(S.CD8, CD8_ordered_max(:,2),CD8_ordered_max(:,3),Recon3DModel,ax(1,2))
ScatterRxnCountFlip(S.Jurkat,Jurkat_ordered_max(:,2),Jurkat_ordered_max(:,3),Recon3DModel, ax(1,3))


function ScatterRxnCountFlip(CellStructure, ordered_max, ordered_vec, Recon3DModel,ax1)


    CellStructure.ModelIndices = kegg2model(string(CellStructure.KEGG), string(Recon3DModel.metKEGGID));
    
    Smat = Recon3DModel.S;
    CellStructure.RxnIndices = cell(numel(CellStructure.ModelIndices),1);
    
    for k = 1:numel(CellStructure.ModelIndices)
    
        mets = CellStructure.ModelIndices{k};   % model metabolite indices
        if isempty(mets)
            CellStructure.RxnIndices{k} = [];
            continue
        end
    
        % Find reactions for all mapped metabolites
        rxnMask = any(Smat(mets,:) ~= 0, 1);
        CellStructure.RxnIndices{k} = find(rxnMask);
    end
    
    CellStructure.RxnCount = cell(numel(CellStructure.ModelIndices),1);
    %%
    for k = 1:numel(CellStructure.ModelIndices)
        CellStructure.RxnCount{k} = size(CellStructure.RxnIndices{k},2);
    end
    %%
    CellStructure.Ranks = cell(numel(CellStructure.ModelIndices),1);
    for k = 1:numel(CellStructure.KEGG)
        CellStructure.Ranks{k} = find(ordered_vec== k);
    end
    
    %% sort everything by metabolite concentration rank
    ranks = str2double(string(CellStructure.Ranks));
    [~, order] = sort(ranks);
    
    
    rxnCount = str2double(string(CellStructure.RxnCount(order)));
    metNames = CellStructure.met(order);
    ranks    = str2double(string(CellStructure.Ranks(order)));
    concentration = ordered_max;
    
    mdl = fitlm(concentration, rxnCount);
    r2 = mdl.Rsquared.Ordinary;
    
   
    
    %% ---- REMOVE SPECIFIC METABOLITES BY NAME ----
    removeMetNames = string({ ...
        %'atp', 'coa', 'gtp','utp','gthrd','gthox','adp','nad','udp','ctp','amp','gdp','dgmp','nadph','gmp','nadh','imp','nadp','camp','fad' ...   % <-- put the exact names you want removed here
    });

    
    
    % keep = NOT in removeMetNames
    keepMask = ~ismember(metNames, removeMetNames);
    keepMask_cf = ismember(metNames, removeMetNames);
    
    % Apply to all aligned vectors
    metNames_cf = metNames(keepMask_cf);
    rxnCount_cf = rxnCount(keepMask_cf);
    ranks_cf    = ranks(keepMask_cf);
    concentration_cf = concentration(keepMask_cf);
   
    metNames_f = metNames(keepMask);
    rxnCount_f = rxnCount(keepMask);
    ranks_f    = ranks(keepMask);
    concentration_f = concentration(keepMask);

    % remove metabolites missing from Recon3DModel
    metNames_f = metNames_f(find(rxnCount_f>0));
    rxnCount_f = rxnCount_f(find(rxnCount_f>0));
    %ranks_f = ranks_f(find(rxnCount_f>0));
    concentration_f = concentration_f(find(rxnCount_f>0));
    
    mdl_f = fitlm(log(rxnCount_f),log(concentration_f));
    r2_f = mdl_f.Rsquared.Ordinary
 
    hold(ax1,"off")
    %yyaxis(ax1,"right")
    scatter(ax1, rxnCount_f, concentration_f,'filled','b')
    %hold on
    hold(ax1, "on")
    xscale(ax1,"log")
    yscale(ax1,"log")
    ylabel(ax1,'Concentration (M)');
    xlabel(ax1, 'Reaction Count');
    %title(ax1, 'Distributions of number of reaction involved with each metabolite')
    % Label points above threshold
    threshold = 150;
    idx = find(rxnCount_f > threshold);

    for i = idx'
        text(ax1,concentration_f(i), rxnCount_f(i), metNames_f(i), ...
            'FontSize', 8, ...
            'VerticalAlignment', 'bottom', ...
            'HorizontalAlignment', 'right');
    end

    % Add R^2 label in top-left corner of axes
    txt = sprintf('R^2 = %.3f', r2_f);
    xpos = 2;
    ypos = 10^-2;

    text(ax1,xpos, ypos, txt, 'FontSize', 12, 'FontWeight','bold','Color','b');

    cofactorMask = ismember(metNames, removeMetNames);

    % Apply to all aligned vectors
    metNames_cf = metNames(cofactorMask);
    rxnCount_cf = rxnCount(cofactorMask);
    ranks_cf    = ranks(cofactorMask);
    concentration_cf = concentration(cofactorMask);

    mdl_cf = fitlm(log(rxnCount_cf),log(concentration_cf));
    r2_cf = mdl_cf.Rsquared.Ordinary;

    
    scatter(ax1,  rxnCount_cf, concentration_cf,'filled','r');
    xscale(ax1,"log")
   
    ylim(ax1, [10^-7 10^-2])


    % Label points above threshold
    threshold = 500;
    idx = find(rxnCount_cf > threshold);

    for i = idx'
        text(ax1,concentration_cf(i), rxnCount_cf(i), metNames_cf(i), ...
            'FontSize', 8, ...
            'VerticalAlignment', 'bottom', ...
            'HorizontalAlignment', 'right');
    end

    % Add R^2 label in top-left corner of axes
    txt = sprintf('R^2 = %.3f', r2_cf);
    xpos = 2;
    ypos = 10^-3;

    text(ax1,xpos, ypos, txt, 'FontSize', 12, 'FontWeight','bold' , 'Color','r');
    legend(ax1, 'Non Cofactor Reactions', 'Cofactor Reactions','Location','southoutside')
    axis(ax1, 'square')
    set(ax1,'TickDir','out')
  

end



