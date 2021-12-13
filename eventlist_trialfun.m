function [trl] = eventlist_trialfun(cfg)

% this function requires the following fields to be specified
% cfg.dataset
% cfg.trialdef.eventtype
% cfg.trialdef.eventvalue
% cfg.trialdef.prestim
% cfg.trialdef.poststim
%%%%

hdr=ft_read_header(cfg.dataset);
evt=ft_read_event(cfg.dataset);

trlidx = find(ismember({evt.type},cfg.trialdef.eventtype) & ismember([evt.value],cfg.trialdef.eventvalue));
trl = [];
for i=1:length(trlidx)
    % add this to the trl definition
    begsample      = -cfg.trialdef.prestim*hdr.Fs+evt(trlidx(i)).sample;
    endsample      = cfg.trialdef.poststim*hdr.Fs+evt(trlidx(i)).sample;
    offset         = cfg.trialdef.prestim*hdr.Fs;
    trl(i,:)       = [round([begsample endsample offset])];
    %event{i} = ScoringLabels{i};
end
