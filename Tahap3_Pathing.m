close all;
clear all;
clc;

direk = uigetdir('*.*', 'Pilih folder datasets');
nomor = 1;
if ~isequal(direk,0)
    
    a = dir(fullfile(direk, '*.*'));
    A = struct2table(a);
    sortedT = sortrows(A, 'date');
    ambil_data = table2struct(sortedT);
    
    
    data2 = [];
    for n = 3 : numel(ambil_data)
        if ambil_data(n).isdir == 0
            nama_data = ambil_data(n).name;
            f_data = fullfile(direk, nama_data);
            data = xlsread(f_data);

            PA = [data(:,5) data(:,6)];
            
            for m = 1 : max(PA)
                [r c] = find(PA == m);
                temp1 = sum(PA(r,2));
                data2(n,m) = temp1;
            end

            fprintf('Proses data ke-%d \n', nomor);
            nomor = nomor + 1;
        end
    end
    
    out_folder = 'Tahap_3_CDF_PDP';
    if ~exist(out_folder, 'dir')
        mkdir(out_folder);
    end
    
    T = table(data2);
    T.Properties.VariableNames = {'All_Path'};
    
    nama_data = sprintf('CDF_PDP.xlsx');
    full_data = fullfile(out_folder, nama_data);
    writetable(T, full_data);
    fprintf('Proses selesai');
    
    
end