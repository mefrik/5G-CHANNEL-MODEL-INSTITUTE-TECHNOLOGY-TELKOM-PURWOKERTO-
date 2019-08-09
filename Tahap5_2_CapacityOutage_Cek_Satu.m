close all;
clear all;
clc;

cd('Tahap_4_PDP_FINAL');
[fname, pname] = uigetfile('*.xlsx', 'Pilih data excel');

if ~isequal(fname, 0) || ~isequal(pname, 0)
    
    EBN = 8; % dB
    banyak = 550000; % Banyaknya pengulangan (TRIAL)
    
    tic;
    awal = tic;
    
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
    fprintf('EBN : %d\n', EBN);
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
            c(i) = log2(1+(abs(psi(i,i)))^2*ebn*M*R);
            c2(i) = log2(1+(abs(psi(i,i)))^2*ebn*M*R2);
            c3(i) = log2(1+(abs(psi(i,i)))^2*ebn*M*R3);
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
    yy1(1,:) = hasiltot{1,:};
    yy2(1,:) = hasiltot{2,:};
    yy3(1,:) = hasiltot{3,:};
    
    waktu = round(toc(awal)/60,3);
    fprintf('>> Waktu Komputasi : %.3f detik\n' , waktu)
    
    T = table(yy1, yy2, yy3);
    T.Properties.VariableNames = {'R1','R2','R3'}
    
    cd ..

    % EXPORT XLSX
    out_folder = 'HASIL';
    if ~exist(out_folder, 'dir');
        mkdir(out_folder);
    end
    
    cd('HASIL');
    
    data = {'EbN0';'Trial';'R = 1'; 'R = 1/2';'R = 3/4'; 'Waktu Komputasi'};
    filename = 'hasil.xlsx';
    sheet = 1;
    xlRange = 'A1';
    xlswrite(filename, data, sheet, xlRange);
    
    switch 1
        case EBN == 0
            data = [EBN; banyak; yy1; yy2; yy3; waktu];
            filename = 'hasil.xlsx';
            sheet = 1;
            xlRange = 'B1';
            xlswrite(filename, data, sheet, xlRange);
        case EBN == 1
            data = [EBN; banyak; yy1; yy2; yy3; waktu];
            filename = 'hasil.xlsx';
            sheet = 1;
            xlRange = 'C1';
            xlswrite(filename, data, sheet, xlRange);
        case EBN == 2
            data = [EBN; banyak; yy1; yy2; yy3; waktu];
            filename = 'hasil.xlsx';
            sheet = 1;
            xlRange = 'D1';
            xlswrite(filename, data, sheet, xlRange);
        case EBN == 3
            data = [EBN; banyak; yy1; yy2; yy3; waktu];
            filename = 'hasil.xlsx';
            sheet = 1;
            xlRange = 'E1';
            xlswrite(filename, data, sheet, xlRange);
        case EBN == 4
            data = [EBN; banyak; yy1; yy2; yy3; waktu];
            filename = 'hasil.xlsx';
            sheet = 1;
            xlRange = 'F1';
            xlswrite(filename, data, sheet, xlRange);
        case EBN == 5
            data = [EBN; banyak; yy1; yy2; yy3; waktu];
            filename = 'hasil.xlsx';
            sheet = 1;
            xlRange = 'G1';
            xlswrite(filename, data, sheet, xlRange);
        case EBN == 6
            data = [EBN; banyak; yy1; yy2; yy3; waktu];
            filename = 'hasil.xlsx';
            sheet = 1;
            xlRange = 'H1';
            xlswrite(filename, data, sheet, xlRange);
        case EBN == 7
            data = [EBN; banyak; yy1; yy2; yy3; waktu];
            filename = 'hasil.xlsx';
            sheet = 1;
            xlRange = 'I1';
            xlswrite(filename, data, sheet, xlRange);
        case EBN == 8
            data = [EBN; banyak; yy1; yy2; yy3; waktu];
            filename = 'hasil.xlsx';
            sheet = 1;
            xlRange = 'J1';
            xlswrite(filename, data, sheet, xlRange);
        case EBN == 9
            data = [EBN; banyak; yy1; yy2; yy3; waktu];
            filename = 'hasil.xlsx';
            sheet = 1;
            xlRange = 'K1';
            xlswrite(filename, data, sheet, xlRange);
        case EBN == 10
            data = [EBN; banyak; yy1; yy2; yy3; waktu];
            filename = 'hasil.xlsx';
            sheet = 1;
            xlRange = 'L1';
            xlswrite(filename, data, sheet, xlRange);
        case EBN == 11
            data = [EBN; banyak; yy1; yy2; yy3; waktu];
            filename = 'hasil.xlsx';
            sheet = 1;
            xlRange = 'M1';
            xlswrite(filename, data, sheet, xlRange);
        case EBN == 12
            data = [EBN; banyak; yy1; yy2; yy3; waktu];
            filename = 'hasil.xlsx';
            sheet = 1;
            xlRange = 'N1';
            xlswrite(filename, data, sheet, xlRange);
        case EBN == 13
            data = [EBN; banyak; yy1; yy2; yy3; waktu];
            filename = 'hasil.xlsx';
            sheet = 1;
            xlRange = 'O1';
            xlswrite(filename, data, sheet, xlRange);
        case EBN == 14
            data = [EBN; banyak; yy1; yy2; yy3; waktu];
            filename = 'hasil.xlsx';
            sheet = 1;
            xlRange = 'P1';
            xlswrite(filename, data, sheet, xlRange);
        case EBN == 15
            data = [EBN; banyak; yy1; yy2; yy3; waktu];
            filename = 'hasil.xlsx';
            sheet = 1;
            xlRange = 'Q1';
            xlswrite(filename, data, sheet, xlRange);
        case EBN == 16
            data = [EBN; banyak; yy1; yy2; yy3; waktu];
            filename = 'hasil.xlsx';
            sheet = 1;
            xlRange = 'R1';
            xlswrite(filename, data, sheet, xlRange);
        case EBN == 17
            data = [EBN; banyak; yy1; yy2; yy3; waktu];
            filename = 'hasil.xlsx';
            sheet = 1;
            xlRange = 'S1';
            xlswrite(filename, data, sheet, xlRange);
        case EBN == 18
            data = [EBN; banyak; yy1; yy2; yy3; waktu];
            filename = 'hasil.xlsx';
            sheet = 1;
            xlRange = 'T1';
            xlswrite(filename, data, sheet, xlRange);
        case EBN == 19
            data = [EBN; banyak; yy1; yy2; yy3; waktu];
            filename = 'hasil.xlsx';
            sheet = 1;
            xlRange = 'U1';
            xlswrite(filename, data, sheet, xlRange);
        case EBN == 20
            data = [EBN; banyak; yy1; yy2; yy3; waktu];
            filename = 'hasil.xlsx';
            sheet = 1;
            xlRange = 'V1';
            xlswrite(filename, data, sheet, xlRange);
        case EBN == 21
            data = [EBN; banyak; yy1; yy2; yy3; waktu];
            filename = 'hasil.xlsx';
            sheet = 1;
            xlRange = 'W1';
            xlswrite(filename, data, sheet, xlRange);
        case EBN == 22
            data = [EBN; banyak; yy1; yy2; yy3; waktu];
            filename = 'hasil.xlsx';
            sheet = 1;
            xlRange = 'X1';
            xlswrite(filename, data, sheet, xlRange);
        case EBN == 23
            data = [EBN; banyak; yy1; yy2; yy3; waktu];
            filename = 'hasil.xlsx';
            sheet = 1;
            xlRange = 'Y1';
            xlswrite(filename, data, sheet, xlRange);
        case EBN == 24
            data = [EBN; banyak; yy1; yy2; yy3; waktu];
            filename = 'hasil.xlsx';
            sheet = 1;
            xlRange = 'Z1';
            xlswrite(filename, data, sheet, xlRange);
        case EBN == 25
            data = [EBN; banyak; yy1; yy2; yy3; waktu];
            filename = 'hasil.xlsx';
            sheet = 1;
            xlRange = 'AA1';
            xlswrite(filename, data, sheet, xlRange);
        case EBN == 26
            data = [EBN; banyak; yy1; yy2; yy3; waktu];
            filename = 'hasil.xlsx';
            sheet = 1;
            xlRange = 'AB1';
            xlswrite(filename, data, sheet, xlRange);
        case EBN == 27
            data = [EBN; banyak; yy1; yy2; yy3; waktu];
            filename = 'hasil.xlsx';
            sheet = 1;
            xlRange = 'AC1';
            xlswrite(filename, data, sheet, xlRange);
        case EBN == 28
            data = [EBN; banyak; yy1; yy2; yy3; waktu];
            filename = 'hasil.xlsx';
            sheet = 1;
            xlRange = 'AD1';
            xlswrite(filename, data, sheet, xlRange);
        case EBN == 29
            data = [EBN; banyak; yy1; yy2; yy3; waktu];
            filename = 'hasil.xlsx';
            sheet = 1;
            xlRange = 'AE1';
            xlswrite(filename, data, sheet, xlRange);
        case EBN == 30
            data = [EBN; banyak; yy1; yy2; yy3; waktu];
            filename = 'hasil.xlsx';
            sheet = 1;
            xlRange = 'AF1';
            xlswrite(filename, data, sheet, xlRange);
    end
    
    % EXPORT MAT
    out_folder = 'MAT FILE';
    if ~exist(out_folder, 'dir');
        mkdir(out_folder);
    end
    
    cd('MAT FILE');
    filename = sprintf('EBN %d.mat',EBN);
    save (filename);
    fprintf('Proses Selesai ... \n');
    cd ..
    cd ..
end