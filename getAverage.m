function data = getAverage(ValidStruct,channel,stimulus,channel32)
    

    data=0;
    for i=2:ValidStruct.stimulus_experiment_number
        if nargin<4
            data = data+getExperiment(ValidStruct,channel,stimulus,i);
        else
            data = data+getExperiment(ValidStruct,channel,stimulus,i,channel32);
        end
    end
    
    data= data/ValidStruct.stimulus_experiment_number;

end