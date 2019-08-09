clear all;
close all;
clc;

cd('HASIL\MAT FILE');
direk = uigetdir('*.mat*', 'Pilih folder data');
nomor = 1;

if ~isequal(direk,0)
    
    %% Load Data
    data = dir(fullfile(direk, '*.mat*'));
    A = struct2table(data);
    sortedT = sortrows(A, 'date');
    L_data = table2struct(sortedT);
    
    tempx = [];
    tempx2 = [];
    tempx3 = [];
    
    tempy = [];
    tempy2 = [];
    tempy3 = [];
    
    for k = 1 : numel(L_data)
        % Take add data
        nama_data = L_data(k).name;
        f_data = fullfile(direk, nama_data);
        % Baca data
        load(f_data);
        % Take data R-1 R-1/2 R-3/4
        tempx{k} = x;
        tempx2{k} = x2;
        tempx3{k} = x3;
        
        tempy{k} = y;
        tempy2{k} = y2;
        tempy3{k} = y3;
        
    end
    
    six = cell2mat(tempx');
    six2 = cell2mat(tempx2');
    six3 = cell2mat(tempx3');
    
    siy = cell2mat(tempy');
    siy2 = cell2mat(tempy2');
    siy3 = cell2mat(tempy3');
    
    cd ..
    % ====================
    %       Figure
    % ====================
    
    % R = 1
    FigW=6;
    FigH=5.6;
    set(figure,'defaulttextinterpreter','latex',...
        'PaperUnits','inches','Papersize',[FigW,FigH],...
        'Paperposition',[0,0,FigW,FigH],'Units','Inches',...
        'Position',[0,0,FigW,FigH])
    for i = 1 : size(tempx,2)
        semilogy([0;tempx{i}],tempy{i})
        hold on;
    end
    hold off;
    str = {'R = 1', '200.000 Channel Realization'};
    text(max(round(unique(six)))/2, tempy{1}(1)+0.01,str,'FontSize',10,'FontName','Arial')
    xlabel('Capacity');
    ylabel('Probability');
    axis([min(six) max(six) min(siy) max(siy)]);
    set(gca,...
        'FontSize',10,'FontName','Arial');
    
    lgd = legend({'0','1','2','3','4','5','6','7','8','9','10','11' ...
        ,'12','13','14','15','16','17','18','19','20'},'Location','southeast');
    title(lgd,'Eb/N0 (dB)');
    lgd.NumColumns = 2;
    
    grid on;
    grid minor;
    out_name = sprintf('Capacittempy_R=1');
    out_name = fullfile(out_name);
    print ('-dpng','-r500', out_name);
    print ('-dpdf','-r500', out_name);
    
    
    % R = 1/2
    FigW=6;
    FigH=5.6;
    set(figure,'defaulttextinterpreter','latex',...
        'PaperUnits','inches','Papersize',[FigW,FigH],...
        'Paperposition',[0,0,FigW,FigH],'Units','Inches',...
        'Position',[0,0,FigW,FigH])
    for i = 1 : size(tempx2,2)
        semilogy([0;tempx2{i}],tempy2{i})
        hold on;
    end
    hold off;
    str = {'R = 1/2', '200.000 Channel Realization'};
    text(max(round(unique(six2)))/2, tempy2{1}(1)+0.01,str,'FontSize',10,'FontName','Arial')
    xlabel('Capacity');
    ylabel('Probability');
    axis([min(six2) max(six2) min(siy2) max(siy2)]);
    set(gca,...
        'FontSize',10,'FontName','Arial');
    
    lgd = legend({'0','1','2','3','4','5','6','7','8','9','10','11' ...
        ,'12','13','14','15','16','17','18','19','20'},'Location','southeast');
    title(lgd,'Eb/N0 (dB)');
    lgd.NumColumns = 2;
    
    grid on;
    grid minor;
    out_name = sprintf('Capacittempy_R=1per2');
    out_name = fullfile(out_name);
    print ('-dpng','-r500', out_name);
    print ('-dpdf','-r500', out_name);
    
    % R = 3/4
    FigW=6;
    FigH=5.6;
    set(figure,'defaulttextinterpreter','latex',...
        'PaperUnits','inches','Papersize',[FigW,FigH],...
        'Paperposition',[0,0,FigW,FigH],'Units','Inches',...
        'Position',[0,0,FigW,FigH])
    for i = 1 : size(tempx3,2)
        semilogy([0;tempx3{i}],tempy3{i})
        hold on;
    end
    hold off;
    str = {'R = 3/4', '200.000 Channel Realization'};
    text(max(round(unique(six3)))/2, tempy3{1}(1)+0.01,str,'FontSize',10,'FontName','Arial')
    xlabel('Capacity');
    ylabel('Probability');
    axis([min(six3) max(six3) min(siy3) max(siy3)]);
    set(gca,...
        'FontSize',10,'FontName','Arial');
    
    lgd = legend({'0','1','2','3','4','5','6','7','8','9','10','11' ...
        ,'12','13','14','15','16','17','18','19','20'},'Location','southeast');
    title(lgd,'Eb/N0 (dB)');
    lgd.NumColumns = 2;
    
    grid on;
    grid minor;
    out_name = sprintf('Capacittempy_R=3per4');
    out_name = fullfile(out_name);
    print ('-dpng','-r500', out_name);
    print ('-dpdf','-r500', out_name);
end