close all;
clear all;
clc;

[fname lname] = uigetfile('.xlsx', 'Pilih data PDP');

if ~isequal(fname, 0) || ~isequal(lname, 0);
    
    ambil_data = fullfile(lname, fname);
    
    data = xlsread(ambil_data);
    
    [r c] = size(data);
    hasil = [];
    for m = 1 : c
        if ~isempty(nonzeros(data(:,m)))
            [a{m} b{m}] = cdfcalc(nonzeros(data(:,m)));
            dataP(m) = prctile(a{m},90);
            avrgnum(m) = sum(data(:,m))/5000;
            avrgdB(m) = 10*log10(avrgnum(m));
        else
            continue
        end        
    end
    gabavrg = [avrgnum; avrgdB];

    %Th
    avrgdBth = avrgdB > -150;
    n = avrgdB.*avrgdBth;
    n(n==0)=[];
            
    % Convert to Numerik
%     n = 10.^(n/10);
    
    %Normalisasi
    avrgdBfinal = (n - max(n));
    
    % Convert to dB
%     avrgdBfinal  = 10*log10(avrgdBfinal);
    
    x = 1 : numel(avrgdBfinal);
    y = avrgdBfinal;
   
    Figure1=figure(1);
    FigW=6;
    FigH=5.6;
    set(Figure1,'defaulttextinterpreter','latex',...
        'PaperUnits','inches','Papersize',[FigW,FigH],...
        'Paperposition',[0,0,FigW,FigH],'Units','Inches',...
        'Position',[0,0,FigW,FigH])
    hStem = stem(x,y,'BaseValue',-90,'Color','blue','MarkerFaceColor','blue');
    set(gca,'Xtick', 1 : 1 : 15)
    set(gca,'XtickLabel', 10 : 10 : 150)
    
    
    X_data = get(hStem, 'XData');
    Y_data = get(hStem, 'YData');
    Labels = num2cell(Y_data);
    for labelID = 1 : 15
        text(X_data(labelID), Y_data(labelID),Labels(labelID) , ...
            'top');
    end

    axis([1-0.5 15+0.8 -90 0]);
    set(gca,...
    'FontSize',10,...
    'FontName','Arial');
    ylabel('Received Power (dB)');
    xlabel('Delay (ms)');
    grid on
    grid minor
    
    out_folder = 'Tahap_4_PDP_FINAL';
    if ~exist(out_folder, 'dir')
        mkdir(out_folder);
    end
    
    T = table(y);
    
    nama_data = sprintf('PDP_CILACAP.xlsx');
    full_data = fullfile(out_folder, nama_data);
    writetable(T, full_data);
    
    out_name_pdf = sprintf('Representative PDP.pdf');
    out_name_pdf = fullfile(out_folder, out_name_pdf);
    out_name_png = sprintf('Representative PDP.png');
    out_name_png = fullfile(out_folder, out_name_png);
    print ('-dpng','-r500', out_name_png);
    print ('-dpdf','-r500', out_name_pdf);
    
end