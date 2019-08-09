clear all;
close all;
clc;

cd('HASIL');
[fname, pname] = uigetfile('*.xlsx', 'Pilih data excel');

if ~isequal(fname, 0) || ~isequal(pname, 0)
    
    %% Input Data
    data = fullfile(pname,fname);
    % Baca data
    a = xlsread(data);
    % Ambil data
    sumbu_x = a(1,:);
    sumbu_r1 = a(3,:); %R = 1
    sumbu_r2 = a(4,:); %R = 1/2
    sumbu_r3 = a(5,:); %R = 3/4
    
    EbNo = 0 : 30; %dalam bentuk dB
    
    % ====================
    %       Figure
    % ====================
    
    figure;
    semilogy(sumbu_x,sumbu_r1, '-.*', sumbu_x,sumbu_r2, '-.diamond', ...
        sumbu_x,sumbu_r3, '-o');
    axis([0 16 10^-4 inf]);
    xlabel('Average E_{b}/N_{0} (dB)');
    ylabel('Outage Probability');
    Figure1=figure(1);
    grid on;
    grid minor;
    legend('R = 1','R = 1/2','R = 3/4')
    str = {'Modulation = BPSK', 'Block length = 128' , ...
        'Trial = 10000'};
    text(6.2,7*10^-1.5,str)
    
    FigW=6;
    FigH=5.6;
    set(Figure1,'defaulttextinterpreter','latex',...
        'PaperUnits','inches','Papersize',[FigW,FigH],...
        'Paperposition',[0,0,FigW,FigH],'Units','Inches',...
        'Position',[0,0,FigW,FigH])
    set(gca,...
        'FontSize',10,...
        'FontName','Arial');
    
    xlabel('Capacity (b/s/Hz)');
    ylabel('Cumulative Distribution Functions');
    
    out_name_png = sprintf('Grafik_R');
    out_name_png = fullfile(out_name_png);
    print ('-dpng','-r500', out_name_png);
end