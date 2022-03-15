function node_name=get_node_name(root_folder_name)

%%get the node name (branch point, Hopf point and fold point) within the
%%current root_folder_name
file=dir(root_folder_name);
node_name.bpt_list={}; node_name.hpt_list={}; node_name.fpt_list={}; node_name.pt_list={};
node_name.pt_last='';
pt_ind=1;
pt_num=0;
for file_ind=1:length(file)
    file_name_1=file(file_ind).name(1);
    switch file_name_1
        case {'b','h','f'}
            ind=regexp(file(file_ind).name,'\d*','Match');
            node_name.([file_name_1,'pt_list']){str2num(ind{1})}=file(file_ind).name;
        case 'p'
            ind=regexp(file(file_ind).name,'\d*','Match');
            node_name.(['pt_list']){pt_ind}=file(file_ind).name;
            pt_ind=pt_ind+1;
            if str2num(ind{1})>pt_num
                pt_num=str2num(ind{1});
            end
    end
    node_name.pt_last=['pt',num2str(pt_num)];
end

end