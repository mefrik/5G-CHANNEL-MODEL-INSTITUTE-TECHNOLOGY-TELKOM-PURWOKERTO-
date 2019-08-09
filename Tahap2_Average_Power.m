close all;
clear all;
clc;

direk = uigetdir('*.*', 'Pilih folder datasets');
nomor = 1;
alpha = 40;
if ~isequal(direk,0)
    
    a = dir(fullfile(direk, '*.*'));
    A = struct2table(a);
    sortedT = sortrows(A, 'date');
    ambil_data = table2struct(sortedT);
    
    for n = 3 : numel(ambil_data)
        if ambil_data(n).isdir == 0
            nama_data = ambil_data(n).name;
            f_data = fullfile(direk, nama_data);
            t_data = xlsread(f_data);
            t_data = abs(t_data);
            
            % Average Power
            % Dilakukan Th
            thdata = t_data < 150;
            npdata = t_data.*thdata;
            [r c] = find(npdata(:,2) ~= 0);
            
            nave = [];
            for i = 1 : size(npdata,1)
                if npdata(i,2) ~= 0
                    power(i) = npdata(i,3);
                else
                    power(i) = 0;
                end
            end
            power = power';
            %         power = padarray(power,size(npdata,1)-max(r),0,'post');
            d = npdata(:,5);
            for m = 1 : max(d)
                [r ~] = find(d == m);
                temp1 = power(r);
                temp2 = sum(temp1);
                temp3 = temp2 / alpha;
                for k = 1 : length(d);
                    if d(r) == d(k)
                        nave(k) = temp3;
                        break
                    elseif isempty(d(r))
                        break
                    end
                end
            end
            
            nave = nave';
            
            if length(nave) ~= size(npdata,1)
                nave = padarray(nave,abs(length(nave) - size(npdata,1)),0,'post');
            end
            
            out_folder = 'Tahap_2_AVERAGE_PDP';
            if ~exist(out_folder, 'dir')
                mkdir(out_folder);
            end
            
            T = table(npdata(:,1), npdata(:,2), npdata(:,3), npdata(:,4), npdata(:,5), nave);
            T.Properties.VariableNames = {'DELAY' 'POWER' 'Numerik' 'OMNI' 'Grouping_index' 'Average'};
            
            nama_data = sprintf('%d.xlsx', nomor);
            full_data = fullfile(out_folder, nama_data);
            writetable(T, full_data);
            
            fprintf('Proses data ke-%d \n', nomor);
            nomor = nomor + 1;
        end
    end
    fprintf('Proses selesai');
end

