% Code modified From 'A parametric model for saccadic eye movement'
% IEEE x_signal Processing in Medicine and Biology Symposium (SPMB), December 2016.
% Modified by Seth Konig 3/14/17

clear
rng(170221,'twister'); %so can regenerate exactly every time

%save settings
savepath = '/Users/ben/Desktop/artificial_scanpaths';
if ~exist(savepath); mkdir(savepath); end

distribution_type = 'uniform';%uniform or gamma

num_simulations = 100;
noise_levels = [0:0.1:1]+.00001;
min_fix_dur = 10; %very minimum fixation duration allowable, may make longer later
min_sac_amp = 2;
max_sac_amp = 6; %if no limit set at >> value e.g. 1000
NSAC = 20; %number of saccades to generate, more saccades -> more eye data -> slower Cluster Fix

sname = [num2str(num_simulations) '_Generated_ScanPaths_' ...
    distribution_type '_NSAC_' num2str(NSAC) '_Saccade_Amplitude_Range_'...
    num2str(min_sac_amp) '_' num2str(max_sac_amp) '.mat'];


%---model parameters---%
P1min = (550-100)/1000; P1max = (550+100)/1000; P1range = P1max - P1min;
P2min = 6-1.5; P2max = 6+1.5;  P2range = P2max - P2min;

xy_buffer = 500; %how long to add no eye movements at begining and end

%---Data To Save---%
xy_amplitude_distribution = cell(1,num_simulations); %effective saccade amplitude size
eyedat = cell(length(noise_levels),num_simulations); %scan paths
true_reference = cell(1,num_simulations); %true saccade times

for sim = 1:num_simulations
    disp(['Simulation #' num2str(sim)])
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%---Generatte Scan Paths---%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %determine whether horizontal, vertical, or diagonal saccade
    randvals = randi(2,2,NSAC);
    generate_xy_saccades = zeros(2,NSAC);
    generate_xy_saccades(:,sum(randvals) == 4) =1; %call diagonal saccades
    generate_xy_saccades(1,randvals(1,:) > randvals(2,:)) = 1; %horizontal saccades
    generate_xy_saccades(2,randvals(2,:) > randvals(1,:)) = 1; %vetical saccades
    generate_xy_saccades(1,sum(generate_xy_saccades) == 0) = 1;%rest call horizontal since more common
    if any(sum(generate_xy_saccades,1) < 1)
        error('What saccade direction?')
    end
    
    %---Generate Horizontal Saccades---%
    x_signal = zeros(xy_buffer,1);
    x_signal_vel = zeros(xy_buffer,1);
    x_reference = x_signal;
    
    y_signal = zeros(xy_buffer,1);
    y_signal_vel = zeros(xy_buffer,1);
    y_reference = y_signal;
    
    for i = 1:NSAC
        Tmax = round(rand(1)*500) + 100;
        Tmin = -100;
        t = Tmin:Tmax;
        t = 2*t'; % 500 Hz, in ms
        eta = rand(1)*P1range + P1min;
        c = rand(1)*P2range + P2min;
        
        mn = min_sac_amp;
        mx = max_sac_amp;

        % tau, corresponds to saccade amplitude
        if strcmpi(distribution_type,'gamma')
            tau = gamrnd(2.4,3.0);
        elseif  strcmpi(distribution_type,'uniform')
            tau = rand(1)*mx + 2;
        end

        p = [eta, c, tau, 0, 0];
        model_amp = 0;
        while model_amp < mn || model_amp > mx
            [model, model_amp, model_pvel, model_vel] = sacc_model_type1('ms', t, p);
            eta = rand(1)*P1range + P1min;
            c = rand(1)*P2range + P2min;
            %get new tau in case doesn't exit while loop
            if strcmpi(distribution_type,'gamma')
                tau = gamrnd(2.4,3.0);
            elseif  strcmpi(distribution_type,'uniform')
                tau = rand(1)*mx + 2;
            end
            p = [eta, c, tau, 0, 0];
        end
        SP = find(model_vel>5, 1, 'first');
        EP = find(model_vel>5, 1, 'last');
        %figure; plot(model); plotcueline('x',[SP,EP]); pause; close(gcf)
        mov_type = zeros(size(model));
        mov_type(SP:EP) = 1;
        model = sign(randn(1)) * model; %direction of saccade
        if abs(x_signal(end) + model(end)) > 30
            model = -model;
        end
        
        model = model'; 
        mov_type = mov_type';
        model_vel = model_vel';
        
        if  generate_xy_saccades(1,i) == 1 %generate horizontal saccade
            model2 = x_signal(end) + model;
            x_signal = [x_signal; model2];
            x_signal_vel = [x_signal_vel; model_vel];
            x_reference = [x_reference; mov_type];
        else
            x_signal = [x_signal;  x_signal(end)*ones(length(model),1)];
            x_signal_vel = [x_signal_vel; zeros(length(model),1)];
            x_reference = [x_reference; zeros(length(model),1)];
        end
        
        if  generate_xy_saccades(2,i) == 1 %generate vertical saccade
            model2 = y_signal(end) + model;
            y_signal = [y_signal; model2];
            y_signal_vel = [y_signal_vel; model_vel];
            y_reference = [y_reference; mov_type];
        else
            y_signal = [y_signal; y_signal(end)*ones(length(model),1)];
            y_signal_vel = [y_signal_vel; zeros(length(model),1)];
            y_reference = [y_reference; zeros(length(model),1)];
        end
    end
    
    %---Add buffer at end so simulation ends with a fixation---%
    x_signal = [x_signal; ones(xy_buffer,1)*x_signal(end)];
    x_signal_vel = [x_signal_vel; zeros(xy_buffer,1)];
    x_reference = [x_reference; zeros(xy_buffer,1)];
    y_signal = [y_signal; ones(xy_buffer,1)*y_signal(end)];
    y_signal_vel = [y_signal_vel; zeros(xy_buffer,1)];
    y_reference = [y_reference; zeros(xy_buffer,1)];
    
    %---Combined X and Y signals---%
    xy = [x_signal'; y_signal']; %combine x and y traces
    xy_reference =  x_reference+y_reference; %combine x and y refrences
    xy_reference(xy_reference > 1) = 1;
    true_reference{sim} = xy_reference;
    
    %---Make generic white noise---%
    white_noise = wgn(2,size(xy,2),1);
    white_noise(1,:) = white_noise(1,:)/std(white_noise(1,:)); %1dva standard devation white noise
    white_noise(2,:) = white_noise(2,:)/std(white_noise(2,:)); %1dva standard devation white noise

    
    %---Make various noise level scan paths---%
    for nl = 1:length( noise_levels)
        eyedat{nl,sim} = xy+noise_levels(nl).*white_noise;
        if 0 && nl==1; %numel(noise_levels)
            x = eyedat{nl,sim}(1,:);
            y = eyedat{nl,sim}(2,:);
            pad = numel(x)-numel(model) - xy_buffer;
            x = x(pad:end);
            y = y(pad:end);
            figure; 
            plotyy(1:numel(x),x,1:numel(x),y); 
            plotcueline('x',[SP,EP]); 
            pause; close(gcf)
        end
    end
    
    %---Calculate Saccade Amplitude from non-noise Path---%
    sac_reference = xy_reference;
    sac_reference = find(sac_reference == 1);
    gaps = findgaps(sac_reference);
    saccade_amplitude = NaN(1,size(gaps,1));
    for g = 1:size(gaps,1)
        gp = gaps(g,:);
        gp(gp == 0) = [];
        if isempty(gp)
            error('what, were did that data go?!')
        elseif all(xy_reference(gp) ~= 1);
            error('what? should be a saccade!')
        end

        %get estimate of saccade amplitude in a way that we can measure later
        fix_xy_before = mean(xy(:,gp(1)-min_fix_dur:gp(1)-1),2);
        fix_xy_after = mean(xy(:,gp(end)+1:gp(end)+min_fix_dur),2);
        saccade_amplitude(g) = sqrt(sum((fix_xy_before-fix_xy_after).^2));
    end
    xy_amplitude_distribution{sim} = saccade_amplitude;
end

disp('saveing...')
save([savepath '/' sname],...
    'xy_amplitude_distribution','eyedat','true_reference')

%% Visualize Generated Saccade Amplitude Distribution
all_amps = [];
for sim = 1:num_simulations
    all_amps = [all_amps  xy_amplitude_distribution{sim}];
end
figure
hist(all_amps,25)
xlabel('Saccade Amplitude (dva)')
ylabel('Saccade Count')