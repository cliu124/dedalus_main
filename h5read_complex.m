function data=h5read_complex(h5_name,data_name)
%%This is a version of h5read that make reading the complex data easy. I do
%%not need to construct the complex data anymore but directly convert them
%%into complex data
%Update 2021/12/06
data=h5read(h5_name,data_name);
try
    data=data.r+1i*data.i;
    disp(['Complex data of ',data_name]);
catch 
    disp(['Real data of ',data_name]);
end
end