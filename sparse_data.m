function data_out = sparse_data(data,start_ind,sparse_ind)
%This function is use to generate sparse data copy that used for plotting

if nargin<2 || isempty(start_ind)
   start_ind=1; 
end

if nargin<3 || isempty(sparse_ind)
    sparse_ind=1;
end

data_out=data;
for data_ind=start_ind:length(data)
    %%make a copy of old data and append them behind the current data
    data_out{length(data)+data_ind-start_ind+1}=data{data_ind};
    data_out{data_ind}.x=data{data_ind}.x(1:sparse_ind:length(data{data_ind}.x));
    data_out{data_ind}.y=data{data_ind}.y(1:sparse_ind:length(data{data_ind}.y));
end

end

