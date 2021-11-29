function trl = contseg_trialfun(cfg)

% this function requires the following fields to be specified
% cfg.dataset
% cfg.trialdef.eventtype
% cfg.trialdef.eventvalue
% cfg.trialdef.prestim
% cfg.trialdef.poststim
%%%%

hdr=ft_read_header(cfg.dataset);
numEpochs=floor(hdr.nSamples/(cfg.trialdef.lengthSegments*hdr.Fs));
trl = [];

for i=1:numEpochs
    % add this to the trl definition
    begsample     = (i-1)*cfg.trialdef.lengthSegments*hdr.Fs+1;
    endsample     = i*cfg.trialdef.lengthSegments*hdr.Fs;
    offset        = 0;
    trl(end+1, :) = [round([begsample endsample offset])];
end

begsample     = (i)*cfg.trialdef.lengthSegments*hdr.Fs+1;
endsample     = hdr.nSamples;

if endsample>begsample
    offset        = 0;
    trl(end+1, :) = [round([begsample endsample offset])];
end
