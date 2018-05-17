clear; clc; close all;

%% Entering input files

directory = {  'autocorrelation/ecoli_aut/',
    'autocorrelation/rhodo_aut/',
    'autocorrelation/ecoli_pert_aut/',
    'autocorrelation/rhodo_pert_aut/'};

description = {'Ecoli (20 C)',
    'Rhodospirillum (20 C)',
    'Ecoli (15 C)',
    'Rhodospirillum (15 C)'};


trainingFraction = .80; % Fraction of data set used for training
validationFraction = 1 - trainingFraction;

totalBacteria = numel(directory);

%% Filtering parameters

fc = 10;
fs = 30;
[b,a] = butter(6, fc/(fs/2));

%% Loading data and calculating distance and phase autocorrelation


for i=1:totalBacteria
    
    files{i} = dir(strcat(directory{i},'*.csv'));
    j=1;
    
    for file = files{i}'
        data{j} = dlmread(strcat(directory{i},file.name),',',1,0);
        j=j+1;
    end
    
    csv{i}=data;
    clear data;
    
    trainingSize{i} = fix(numel(files{i})*trainingFraction);
    
    for j=1:numel(csv{i})
        csv{1,i}{j}(:,8) = filter(b, a, csv{1,i}{j}(:,8));
        csv{1,i}{j}(:,9) = filter(b, a, csv{1,i}{j}(:,9));
        distance{1,i}{j} = calculateDistance(csv{1,i}{j}(:,8),csv{1,i}{j}(:,9)); % With DC blocker
        %distance{1,i}{j} = filter(b, a, distance{1,i}{j});
        phase{1,i}{j} = calculatePhase(csv{1,i}{j}(:,8),csv{1,i}{j}(:,9)); % With DC blocker
        distanceAutocorr{1,i}{j} = autocorr(distance{i}{1,j},299,1);
        phaseAutocorr{1,i}{j} = autocorr(phase{i}{1,j},299,1);
    end
    
    distanceAutocorrAvg{i} = zeros(300,1);
    phaseAutocorrAvg{i} = zeros(300,1);
    
    [distanceTraining{i}, trainingSet{i}]=datasample(distanceAutocorr{i},trainingSize{i},'Replace',false);
    phaseTraining{i}=phaseAutocorr{i}(:,trainingSet{i}); % Use the same set of data to train for distance and phase
    
    for j=1:trainingSize{i}
        distanceAutocorrAvg{i} = distanceAutocorrAvg{i} + distanceTraining{i}{1,j} / trainingSize{i};
        phaseAutocorrAvg{i} = phaseAutocorrAvg{i} + phaseTraining{i}{1,j} / trainingSize{i};
    end
    clear file;
end

%% Visualization of Distance and Phase Autocorrelation

figure();
set(gcf,'color','white')

subplot(1,2,1);
hold on;
title('Distance Autocorrelation');
for i=1:totalBacteria
    plot(distanceAutocorrAvg{i},'-o','MarkerSize',2);
end
xlim([0 100]);
ylim([-0.1 0.5]);
legend(description);
set(gca,'FontSize',16);
xlabel('Distance (lags)');
ylabel('Autocorrelation');

subplot(1,2,2);
hold on;
title('Phase Autocorrelation');
for i=1:totalBacteria
    plot(phaseAutocorrAvg{i},'-o','MarkerSize',2);
end
xlim([0 30]);
ylim([-0.2 0.2]);
legend(description);
set(gca,'FontSize',16);
xlabel('Phase (lags)');
ylabel('Autocorrelation');

%% Visualization of bacteria movement for selected samples

for i=1:totalBacteria
    trackVisual=datasample(csv{i},5,'Replace',false);
    figure(); hold on;
    set(gcf,'color','white')
    set(gca,'FontSize',16);
    title(description{i}); xlabel('X (\mum)'); ylabel('Y (\mum)');
    for j=1:5
        plot(trackVisual{1,j}(:,4)*0.24,trackVisual{1,j}(:,5)*0.24);
    end
end

%% Calculating cross-correlation

distanceMaxLag = 30;
phaseMaxLag = 20;

x={0,0,0,0};

for i=1:totalBacteria
    
    figure('position', [100, 100, 1200, 1000]); set(gcf,'color','white'); hold on;
    
    for j=1:numel(files{i})
        
        if ~ismember(j,trainingSet{i})
            
            x{i}=x{i}+1;
            
            for k=1:totalBacteria
                % Max lag: 50
                distanceCrossCorrelation{i}{k}{x{i}} = xcorr(distanceAutocorr{1,i}{x{i}},distanceAutocorrAvg{k},distanceMaxLag,'coeff');
                % Max lag: 20
                phaseCrossCorrelation{i}{k}{x{i}} = xcorr(phaseAutocorr{1,i}{x{i}},phaseAutocorrAvg{k},phaseMaxLag,'coeff');
                
                subplot(4,2,2*k-1); hold on; plot(distanceCrossCorrelation{i}{k}{x{i}});
                title(strcat(description{i},{' '}, '\otimes',{' '}, description{k}));
                xlim([0 distanceMaxLag*2]); xlabel('Distance (lags)'); ylabel('Cross-correlation'); set(gca,'FontSize',14);
                
                subplot(4,2,2*k); hold on; plot(phaseCrossCorrelation{i}{k}{x{i}});
                title(strcat(description{i},{' '}, '\otimes',{' '}, description{k}));
                xlim([0 phaseMaxLag*2]); xlabel('Phase (lags)'); ylabel('Cross-correlation');  set(gca,'FontSize',14);
                
            end
        end
    end
end

%% Distance statistics

for i=1:totalBacteria
    
    
    %averageDistance{i} = 0;
    distanceCumulative{i} = 0;
    
    for j=1:numel(distance{i})
        
        distanceCumulative{i} = [distanceCumulative{i} ; distance{i}{j}];
        
        
    end
    distanceCumulative{i }= distanceCumulative{i} * 0.24; % convert in mum
    distanceAvg{i} = mean(distanceCumulative{i});
    distanceStd{i} = std(distanceCumulative{i});
    distanceSkew{i}{j}=skewness(distanceCumulative{i});
    %errorbar(i,distanceAvg{i},distanceStd{i});
    
end
figure(); hold on; set(gcf,'color','white');set(gca,'FontSize',16);
x1 = distanceCumulative{1};
x2 = distanceCumulative{2};
x3 = distanceCumulative{3};
x4 = distanceCumulative{4};
x = [x1;x2;x3;x4];
g = [ones(size(x1)); 2*ones(size(x2)); 3*ones(size(x3)); 4*ones(size(x4))];
boxplot(x,g);
%hLegend = legend(findall(gca,'Tag','Box'), description, 'Orientation', 'horizontal');
%% Plot histogram

for i=1:totalBacteria
    distanceHistogram{i}=datasample(distance{i},trainingSize{i},'Replace',false);
    figure();
    set(gcf,'color','white');
        
    subplot(1,2,1); hold on;
    title(description{i}); set(gca,'FontSize',16);
    xlabel('Distance per step (\mum)'); xlim([0 inf]);
    ylabel('Frequency');
    for j=1:numel(distanceHistogram{i})
        distanceHistogram{i}{j} = distanceHistogram{i}{j} * 0.24; % convert in mum
        histogram(distanceHistogram{i}{j},50);
    end
    
    subplot(1,2,2); hold on;
    title(description{i}); set(gca,'FontSize',16);
    xlabel('Log of Distance per step (\mum)'); xlim([0 inf]);
    ylabel('Frequency');
    for j=1:numel(distanceHistogram{i})
        histogram(log(distance{i}{j}),50);
    end
    
end

clear i j file;