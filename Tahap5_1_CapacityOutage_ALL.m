close all;
clear all;
clc;


cd('Tahap_4_PDP_FINAL');
[fname, pname] = uigetfile('*.xlsx', 'Pilih data excel');

if ~isequal(fname, 0) || ~isequal(pname, 0)
    
    banyak = 500000; % Banyaknya pengulangan (TRIAL)
    
    %% Input Data
    data = fullfile(pname,fname);
    % Baca data
    inputdb = xlsread(data);
    % Banyaknya elemen input
    n = numel(inputdb);
    % Cari index input tidak sama dengan 0
    % [~ c] = find(input~=0);
    % Convert Input to Numerik
    input =10.^(inputdb/10);
    
    %% Banyaknya nilai Eb/N0
    EBN = 0:30; % dB
    % Convert dB ke Numerik
    ebn =10.^(EBN/10);
    
    %% Modulasi BPSK
    m = 2;
    % Indeks modulasi
    M = log2(m);
    
    %% Transformasi (FFT)
    N = 128; %FFT Size
    % Bulatkan keatas berdasarkan numerology 3
    Q = ceil (0.57/8.33*N);
    % Rate broadband
    R = 1*((N*M)/((N*M)+Q));
    R2 = (1/2)*((N*M)/((N*M)+Q));
    R3 = (3/4)*((N*M)/((N*M)+Q));
    
    % hh = ?
    hh = [];
    for i = 1 : n
        hh(i) = sqrt(input(i));
    end
    lh = length(hh); % == length(input)
    tb = N-lh;
    % Buat Matrix FFT
    F = dftmtx(N)/sqrt(N);
    % Invers FFT Matrix (hermitian matriks)
    Fh = F';
    
    yy1 = [];
    yy2 = [];
    yy3 = [];
    H = [];
    Hc = [];
    c = [];
    c2 = [];
    c3 = [];
    
    %% Calculation
    for kasus = 1 : length(ebn)
        fprintf('EBN : %d\n', EBN(kasus));
        for ulang = 1 : banyak
            for jj = 1 : lh
                H(jj) = hh(jj)*(1/sqrt(2)*[randn(1,1)+1i*randn(1,1)]);
            end
            
            % Circulant matriks
            input = reshape(([H zeros(1,tb)]),[],1);
            L = length(input);
            for ii = 1 : L
                Hc(:,ii) = circshift(input,ii-1);
            end
            % Mencari nilai \psi
            psi=F*Hc*Fh;
            for i = 1 : N;
                % Cari nilai capacity
                c(i) = log2(1+(abs(psi(i,i)))^2*ebn(kasus)*M*R);
                c2(i) = log2(1+(abs(psi(i,i)))^2*ebn(kasus)*M*R2);
                c3(i) = log2(1+(abs(psi(i,i)))^2*ebn(kasus)*M*R3);
            end
            % Cari nilai capa
            capa = sum(c)/N;
            capa2 = sum(c2)/N;
            capa3 = sum(c3)/N;
            % Simpan seluruh 'capa' dalam Hasilakhir sebanyak 'ulang'
            Hasilakhir(ulang) = capa;
            Hasilakhir2(ulang) = capa2;
            Hasilakhir3(ulang) = capa3;
            
        end
        % Hitung cumulative dist function dan cari indexnya
        [y x] = cdfcalc(Hasilakhir);
        [y2 x2] = cdfcalc(Hasilakhir2);
        [y3 x3] = cdfcalc(Hasilakhir3);
        
        %         %% Plotting
        %         semilogy([0;x],y)
        %         hold on
        pos1 = find(x(:,1)>=(R*M));
        pos2 = find(x2(:,1)>=(R2*M));
        pos3 = find(x3(:,1)>=(R3*M));
        
        allpos = {pos1' ; pos2' ; pos3'};
        [b ~] = size(allpos);
        int = 1;
        hasil1 = [];
        hasiltot = [];
        for m = 1 : b
            if isempty(allpos{m})
                hasil1 = 0;
            else
                hasil1 = y((allpos{m}(1,1)),1);
            end
            hasiltot{m,1} = hasil1;
        end
        yy1(1,kasus) = hasiltot{1,:};
        yy2(1,kasus) = hasiltot{2,:};
        yy3(1,kasus) = hasiltot{3,:};
    end
    fprintf('Proses Selesai ... \n');
    
    cd ..
    % EXPORT PNG
    out_folder = 'HASIL';
    if ~exist(out_folder, 'dir');
        mkdir(out_folder);
    end
  
    figure
    semilogy(EBN,yy1,'--* r');
    hold on;
    semilogy(EBN,yy2,'-.diamond b');
    hold on;
    semilogy(EBN,yy3,'-o gr');
    xlabel('Average E_{b}/N_{0} (dB)');
    ylabel('Outage Probability');
    Figure1=figure(1);
    grid on;
    grid minor;
    legend('R = 1','R = 1/2','R = 3/4')
    str = {'Modulation = BPSK', 'Block length = 128' , ...
        'Trial = 550000'};
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
    
    out_name_png = sprintf('Hasil 1');
    out_name_png = fullfile(out_folder, out_name_png);
    print ('-dpng','-r500', out_name_png);
    
    % EXPORT MAT
    cd ('HASIL');
    out_folder = 'MAT FILE ALL';
    if ~exist(out_folder, 'dir');
        mkdir(out_folder);
    end
    
    cd('MAT FILE ALL');
    filename = sprintf('Capacity_ALL.mat');
    save (filename);
    
    
    %     Figure2=figure(2);
    %     grid on;
    %     grid minor;
    %     str = {'Rate= 1', 'Modulation= BPSK', 'Block length= 128' ,'Trial= 10'};
    %     text(6.2,7*10^-1.5,str)
    %     FigW=6;
    %     FigH=5.6;
    %     set(Figure2,'defaulttextinterpreter','latex',...
    %         'PaperUnits','inches','Papersize',[FigW,FigH],...
    %         'Paperposition',[0,0,FigW,FigH],'Units','Inches',...
    %         'Position',[0,0,FigW,FigH])
    %     set(gca,...
    %         'FontSize',10,...
    %         'FontName','Arial');
    %
    %
    %     out_name_png = sprintf('Hasil 2');
    %     out_name_png = fullfile(out_folder, out_name_png);
    %     print ('-dpng','-r500', out_name_png);    
end