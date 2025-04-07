% Define the filename
filename = '3d_array.xlsx';
%filename = 'time_data/3d_array_MBC3.xlsx';


% Preallocate a 3D matrix to store the data
data = zeros(6, 11, 8);


% Loop over each sheet and read the data
for sheetIdx = 1:8
    % Read the data from the current sheet
    sheetData = readmatrix(filename, 'Sheet', sheetIdx);
    
    % Store the data in the 3D matrix
    data(:, :, sheetIdx) = sheetData;
end

%MBC2_total_spk_time_mean = flip([0.24345713, 0.18276299, 0.14061479, 0.12289202, 0.10208609, 0.08648075, 0.06745413, 0.0542521]);
%MBC2_total_spk_time_std = flip([0.05264071, 0.02154329, 0.01921525, 0.02026622, 0.01523622,0.00995559, 0.01591446, 0.01701071]);

AUT1_total_spk_time_mean = flip([0.29122745, 0.18590999, 0.13827705, 0.1084794 , 0.08625309, 0.07869753, 0.0589961 , 0.05215938]);
AUT1_total_spk_time_std = flip([0.08891172, 0.03762442, 0.02698712, 0.02369915, 0.02391646, 0.02080413, 0.01662813, 0.01662845]);
%MBC3_total_spk_time_mean = flip([0.26597142, 0.18978178, 0.15148321, 0.11619238, 0.09437194, 0.07797848, 0.06229439, 0.0419264]);
%MBC3_total_spk_time_std = flip([0.06649339, 0.0325574 , 0.02644631, 0.01707068, 0.02176769, 0.0233892 , 0.02589208, 0.02326048]);

total_spk_time_mean = AUT1_total_spk_time_mean;
total_spk_time_std = AUT1_total_spk_time_std;

%total_spk_time_mean = MBC3_total_spk_time_mean;
%total_spk_time_std = MBC3_total_spk_time_std;


%each epoch use different std for normally distributed c_l

M = 8; %number of speaker
avg = 1;
total_speaking_time_l = zeros(M*(M-1),M);
ct = 1;


for k = 1:M

for j = 1:M

if j~=k

s = rand(M,1);
s(k) = 1;
s(j) = 0;

%s = [1;0;0.3];
s_ind = k; %index of the current speaker
%s = [1;0;0.26];
ds = zeros(M);
ds_new = zeros(M);

power = 4.2;


p = linspace(0.1,0.9,M);
mu = 1;
sigma = eps;
cl = avg-(power-1)/(power-2)+(1-p).^(1/(1-power))';

cs = 2*avg - cl;

alpha = 0.9;

N = 10^5;

t = zeros(1,N);
s_list = zeros(length(s),N);
s_list(:,1) = s;
s_ind_list = zeros(1,N);
s_ind_list(1) = s_ind;
ds_list = zeros(1,N);
total_speaking_time = zeros(1,M);

for i = 1:N-1
    tmp = (1-s)./(cl+cs(s_ind)); %intersecting times for speaker with each listener 
    tmp(s_ind) = max(tmp); %the intersecting time for speaker is infinite
    [deltat, ind] = min(tmp); %find the smallest interecting time
    t(i+1) = t(i) + deltat;
    s = alpha * (s + deltat*cl);
    s(s_ind) = 0; %The previous speaker has 0 desire to speak
    s(ind) = 1; %The current speaker has 1 desire to speak
    total_speaking_time(s_ind) = total_speaking_time(s_ind)+deltat; %update the total speaking time of each speaker

    s_ind = ind; %update speaker
    s_list(:,i+1) = s;
    s_ind_list(i+1) = ind;

    ds_new(s_ind_list(i),s_ind_list(i+1)) = s_list(s_ind_list(i),i) - s_list(s_ind_list(i+1),i); 
    ds_list(i+1) = norm(ds_new-ds);
    ds = ds_new;
end

total_speaking_time_l(ct,:) = total_speaking_time/sum(total_speaking_time);
ct = ct + 1;
end
end
end

% Calculate the mean, min, and max of each column
column_means = mean(total_speaking_time_l, 1); % Mean of each column
column_mins = min(total_speaking_time_l, [], 1); % Minimum value of each column
column_maxs = max(total_speaking_time_l, [], 1); % Maximum value of each column

% Calculate the error bars
lower_error = column_means - column_mins; % Difference between mean and min
upper_error = column_maxs - column_means; % Difference between max and mean

model_column_means = column_means;

figure(3)
hold on
%fill([1:M, fliplr(1:M)], ...
%     [column_maxs, fliplr(column_mins)], ...
%     [0.8, 0.8, 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.5); % Shaded region
%plot(1:M,total_speaking_time_l(min_index,:),'r','HandleVisibility','off')
for i = 1:size(data, 2)
    tmp = data(1,i,:);
    tmp = reshape(tmp,[8,1]);
    scatter(flip(1:M),tmp,'blue','filled')
end
scatter(1:M, column_means, 80, 'red','*');
plot(1:M, column_means, 'r--*')
hold off
box on
xlabel('Speaker index','FontSize',20)
ylabel('Fraction of time','FontSize',20)
%legend('model','data','Location','northwest','FontSize',12)
xlim([0.5,M+0.5])
ax=gca;
ax.FontSize = 20;
% Get the current y-ticks
current_ticks = yticks();

% Get the minimum and maximum y-axis values
y_min = min(current_ticks);
y_max = max(current_ticks);

% Create new y-tick labels
new_labels = cell(size(current_ticks)); % Start with empty labels
new_labels{1} = num2str(y_min);        % Set label for the min
new_labels{end} = num2str(y_max);      % Set label for the max

% Apply the new labels
yticklabels(new_labels);

xticklabels([]);
