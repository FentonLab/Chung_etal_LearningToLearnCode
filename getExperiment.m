function data = getExperiment(ValidStruct,channel,stimulus,experiment,channel32)
    
if nargin<5||channel32==0

    if experiment~=1
        index_stimulus=stimulus/ValidStruct.stimulus_step+1;
        index_num_stimulus=(ValidStruct.stimulus_experiment_number*(index_stimulus-1)+experiment-index_stimulus);
        data=ValidStruct.data(channel,:,index_num_stimulus);
    else
       data=[];
       error('You can not access to the first experiment');
    end
else
    index_stimulus=stimulus/ValidStruct.stimulus_step+1;
    index_num_stimulus=(ValidStruct.stimulus_experiment_number*(index_stimulus-1)+experiment);
    data=ValidStruct.data(channel,:,index_num_stimulus);
end

end