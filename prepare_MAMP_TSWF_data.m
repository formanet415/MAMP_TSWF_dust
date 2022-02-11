file = 'solo_dust.txt';
fileid = fopen(file, 'r');

tab = readtable(file);

tswf_epoch = [];
tswf_samp_rate = [];
tswf_waveform = zeros(3000, 3, 32768);
tswf_srf = zeros(3000, 2, 32768);
tswf_samps_per_ch = [];
tswf_channel_ref = zeros(3000,3);

mamp_data = zeros(3000,1000,4);
mamp_epoch = zeros(3000,1000);
mamp_recs = [];
mamp_samp_rate = zeros(3000,1000);


for r=1:length(tab.Var1)
    date = datenum(tab.Var1(r));
    mamp = tdscdf_load_l2_surv_mamp(date, 1);
    if isempty(mamp)
        continue;
    end
    tswf = tdscdf_load_l2_surv_tswf(date, 1);
    if isempty(tswf)
        continue;
    end
    idx = strsplit(convertCharsToStrings(tab.Var2{r}), ',');
    for i=1:length(idx)
        inds = str2num(strrep(idx(i),'-',':'));
        for ind = inds
            if ind>length(tswf.epoch)
                continue
            end
            t0 = tswf.epoch(ind);
            begindex = length(mamp.epoch(mamp.epoch<t0));
            t1 = addtodate(t0,ceil(1e3*tswf.samples_per_ch(ind)/tswf.samp_rate(ind)),'millisecond');
            endindex = 1+length(mamp.epoch)-length(mamp.epoch(mamp.epoch>t1));
            if endindex>length(mamp.epoch) || begindex == 0
                continue
            end

            tswf_epoch(end+1) = tswf.epoch(ind);
            len = length(tswf_epoch);
            tswf_samp_rate(end+1) = tswf.samp_rate(ind);
            spch = tswf.samples_per_ch(ind);
            tswf_samps_per_ch(end+1) = spch;
            tswf_waveform(len, :, 1:spch) = tswf.data(:,1:spch,ind);
            tswf_srf(len, :, 1:spch) = convert_to_SRF(tswf,ind);
            tswf_channel_ref(len, :) = tswf.channel_ref(:,ind);
            
            mamp_recs(end+1) = endindex + 1 - begindex;
            mamp_data(len, 1:mamp_recs(end), end) = nan;
            mamp_data(len, 1:mamp_recs(end), 1:size(mamp.data, 2)) = mamp.data(begindex:endindex, :);
            mamp_epoch(len, 1:mamp_recs(end)) = mamp.epoch(begindex:endindex);
            mamp_samp_rate(len, 1:mamp_recs(end)) = mamp.samp_rate(begindex:endindex);
        end
    end
    
end

save('autosave.mat')