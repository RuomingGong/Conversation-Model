% Define the filename
filename = '/data/speakingtime_MBC2_new.xlsx';

% Preallocate a 3D matrix to store the data
data = zeros(6, 14, 8); % mbc2 data shape: stats_dim / num_session / num_speaker 


% Loop over each sheet and read the data
for sheetIdx = 1:8
    % Read the data from the current sheet
    sheetData = readmatrix(filename, 'Sheet', sheetIdx);
    
    % Store the data in the 3D matrix
    data(:, :, sheetIdx) = sheetData;
end

total_spk_time_mean = squeeze(mean(data(1,:,:)));
total_spk_time_std = squeeze(std(data(1,:,:)));

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

s_ind = k; %index of the current speaker
ds = zeros(M);
ds_new = zeros(M);

power = 4.4;


%p = linspace(1/(2*M),1-1/(2*M),M);
p = linspace(0.1, 0.9, M);
mu = 1;
sigma = eps;
cl = avg-(power-1)/(power-2)+(1-p).^(1/(1-power))';

%cs = 1./cs;
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



% Calculation by stochastic model 

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

power = 4.4;


%p = linspace(1/(2*M),1-1/(2*M),M);
p = linspace(0.1,0.9,M);
mu = 1;
sigma = eps;
cl = avg-(power-1)/(power-2)+(1-p).^(1/(1-power))';

%cs = 1./cl;
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

    %random switch
    [max1, idx1] = max(s);
    array_excluding_max1 = s;
    array_excluding_max1(idx1) = 0; % Replace the smallest element with Inf
    [max2, idx2] = max(array_excluding_max1);
    top_two_indices = [idx1, idx2];
    random_index = top_two_indices(randi(length(top_two_indices)));
    ind = random_index;


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

figure(3)
darkGreen = [0, 0.3, 0];
hold on
for i = 1:size(data, 2)
    tmp = data(1,i,:);
    tmp = reshape(tmp,[8,1]);
    scatter(1:M,tmp,'blue','filled')
end
plot(1:M, flip(column_means), 'k-', LineWidth=3)
plot(1:M, flip(model_column_means), 'k--', LineWidth=3)
hold off
box on
xlabel('Speaker index','FontSize',20)
ylabel('Fraction of time','FontSize',20)
xlim([0.5,M+0.5])
ylim([0,0.5])
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

xticks(1:M)

% Create new y-tick labels
new_labels = cell(size(current_ticks)); % Start with empty labels
new_labels{1} = num2str(1);        % Set label for the min
new_labels{M} = num2str(M);      % Set label for the max
xticklabels(new_labels);

h = gcf;
set(h, 'PaperUnits','centimeters');
set(h, 'Units','centimeters');
pos=get(h,'Position');
set(h, 'PaperSize', [pos(3) pos(4)]);
set(h, 'PaperPositionMode', 'manual');
set(h, 'PaperPosition',[0 0 pos(3) pos(4)]);

saveas(gcf,'\results\modelfit.pdf')

