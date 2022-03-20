function RSA_heatmap_plot (Data_Array, saveName)
    
    script_path = '/home/chendanni/Documents/Norms/analysis/MyScripts';
    addpath (genpath(script_path));
    
    imAlpha=ones(size(Data_Array));
    imAlpha(isnan(Data_Array))=0;
    
    
    CustomXLabels = string ([1:8]);
    for i = 1:7 CustomXLabels(i) = ""; end
    CustomXLabels(1) = 'In-Higher';
    CustomXLabels(2) = 'In-Lower';
    CustomXLabels(3) = 'In-Consistent';
    CustomXLabels(4) = 'Out-Higher';
    CustomXLabels(5) = 'Out-Lower';
    CustomXLabels(6) = 'Out-Consistent';
    % CustomXLabels(7) = 'Control';
    CustomXLabels(7) = 'Morph';

    figure;
    set (gcf, 'Position', [300 300 800 700]);
    
    imagesc(Data_Array,'AlphaData',imAlpha);
    cmp = customcolormap_preset('red-yellow-blue');
    colormap (cmp)
    colorbar
    title (saveName)
    set(gca,'xtick',[5:10:80],'xticklabel',CustomXLabels,...
        'ytick',[5:10:80],'yticklabel',CustomXLabels)
    xtickangle(45)
    hold on 
    for i = [0:10:80]
        for j = [0:10:80]
            rectangle ('Position', [i+.5 j+.5  10 10], 'LineWidth', 2);
        end
    end
    
    saveas(gcf,strcat(saveName, '.png'));
    pause (0.1)
    close all;
    
end