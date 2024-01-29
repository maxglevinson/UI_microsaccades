function scratch_investigate_single_trial(trial_data)
% looking into single trial eye position / movement

%trial_data should be: trial_data from normal full eyelink analysis
%and have trial event indices computed
% subj 202, trial 122 (part of block 2)
% subj 201, trial 140 (part of block 3)

xPos = trial_data(2:3, :);
yPos = trial_data(4:5, :);
       
         %% engbert
%         [EyeSaccades, BinocularSaccades] = microsaccades_engbert_kliegel(data);
%         
%         % plot x and y positions
%         figure;
%         p = plot(xPos', yPos', 'Color', [1; 0; 0]);
%         p(1).Color = [1 0 0];
%         p(2).Color = [0 0 1];
%         xlim([0, 2560]);
%         ylim([0, 1440]);
%         % bold microsaccades
%         hold on;
%         for ms = 1:numel(BinocularSaccades.msStarts)
%             ms_indices = BinocularSaccades.msStarts(ms):(BinocularSaccades.msEnds(ms));
%         ms_xPos = xPos(:, ms_indices);
%         ms_yPos = yPos(:, ms_indices);
%         ms_p = plot(ms_xPos', ms_yPos', 'linewidth', 2);
%         ms_p(1).Color = [1 0 0];
%         ms_p(2).Color = [0 0 1];
%         end
%         
%         
%         figure; px = plot(xPos'); hold on; py = plot(yPos'); % both pos over time
%         px(1).Color = [1 0 0]; py(1).Color = [1 0 0];
%         px(2).Color = [0 0 1]; py(2).Color = [0 0 1];
%         % bold microsaccades
%         for ms = 1:numel(BinocularSaccades.msStarts)
%                         ms_indices = BinocularSaccades.msStarts(ms):(BinocularSaccades.msEnds(ms));
%         ms_xPos = xPos(:, ms_indices);
%         ms_yPos = yPos(:, ms_indices);
%         ms_px = plot(ms_indices, ms_xPos', 'linewidth', 2); hold on; ms_py = plot(ms_indices, ms_yPos', 'linewidth', 2);
%         ms_px(1).Color = [1 0 0]; ms_py(1).Color = [1 0 0];
%         ms_px(2).Color = [0 0 1]; ms_py(2).Color = [0 0 1];
%         end
        
        %% martinez conde
        [microsaccade_details, bi_microsaccade_details, xyVelocity] = microsaccades_martinez_conde_2006(trial_data);
        
        % plot x and y positions
        figure;
        p = plot(xPos', yPos');
        p(1).Color = [1 0 0];
        p(2).Color = [0 0 1];
        xlim(2560/2 + [-50, 50]);%([0, 2560]);
        ylim(1440/2 + [-50, 50]);%([0, 1440]);
        % bold microsaccades
        hold on;
        for ms = 1:numel(bi_microsaccade_details.start)
            ms_indices = bi_microsaccade_details.start(ms):(bi_microsaccade_details.start(ms)+bi_microsaccade_details.duration(ms));
        ms_xPos = xPos(:, ms_indices);
        ms_yPos = yPos(:, ms_indices);
        ms_p = plot(ms_xPos', ms_yPos', 'linewidth', 2);
        ms_p(1).Color = [1 0 0];
        ms_p(2).Color = [0 0 1];
        end
        
        
%         % plot both pos over time
%         t = trial_data(1, :);
%         figure; px = plot(t, xPos'); hold on; py = plot(t, yPos');
%         px(1).Color = [1 0 0]; py(1).Color = [1 0 0];
%         px(2).Color = [0 0 1]; py(2).Color = [0 0 1];
%         % bold microsaccades
%         for ms = 1:numel(bi_microsaccade_details.start)
%                         ms_indices = bi_microsaccade_details.start(ms):(bi_microsaccade_details.start(ms)+bi_microsaccade_details.duration(ms));
%         ms_xPos = xPos(:, ms_indices);
%         ms_yPos = yPos(:, ms_indices);
%                                 ms_indices = ms_indices + trial_data(1)-1; % align to first t=~-1000
%         ms_px = plot(ms_indices, ms_xPos', 'linewidth', 2); hold on; ms_py = plot(ms_indices, ms_yPos', 'linewidth', 2);
%         ms_px(1).Color = [1 0 0]; ms_py(1).Color = [1 0 0];
%         ms_px(2).Color = [0 0 1]; ms_py(2).Color = [0 0 1];
%         end
%         ys = get(gca, 'ylim');
%         plot(repmat(stim_start_idx - trial_start_idx + trial_data(1), 1, 2), ys, 'k'); % plot stim onset time
%         plot(repmat(button_press_idx - trial_start_idx + trial_data(1), 1, 2), ys, 'k'); % plot button press time