%%
function indices_out = get_event_indices(trialstruct, buffer_pre, buffer_post, event)
% buffers: how many timepoints before and after blink/saccade start/end to include
indices_out = [];
noevent = 1; % default to not extracting anything, unless events exist
switch event
    case 'blink'
        if isfield(trialstruct.Blinks, 'sttime')
            data = trialstruct.Blinks;
            noevent = 0;
        end
    case 'saccade'
        if isfield(trialstruct.Saccades, 'sttime')
            data = trialstruct.Saccades;
            noevent = 0;
        end
end
if ~noevent
    events2use = data.sttime > 0; % avoid bug of blink at t=0
    if strcmp(event, 'saccade')
        events2use = events2use & data.ampl > 2; % also only use saccades > 2 dva
    end
    for ev = 1:numel(data.sttime)
        if events2use(ev)
            indices_out = [indices_out, (double(data.sttime(ev))-buffer_pre):(double(data.entime(ev))+buffer_post)];
        end
    end
if data.sttime(1) == 0 % catch bug where eyes weren't tracked at all in the trial. Gives a blink lasting whole trial essentially
    if data.time(1) > (trialstruct.KeyEvents.STIM_OFFSET - trialstruct.KeyEvents.STIM_ONSET - 1000) % if blink lasts basically whole stim period
        indices_out = [indices_out, 0:double(data.entime(1))];
    end
end
else
    indices_out = [];
end