classdef my_default
    %MY Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        no_ylabel=0;
        folder_name=0;
        plot_config=struct;
        Nu=0;
        Nu_S=0;
        obj=struct;
        par_num=4;
        point_list=10000;
        bif_type={'bpt'};
        depth=1;
    end
    
    methods
        function obj = my_default(obj)
            obj.no_ylabel=0;
            obj.folder_name=0;
            obj.plot_config=struct;
            obj.plot_config.visible=0;
            obj.plot_config.print=0;
            obj.plot_config.post=0;
            obj.Nu=0;
            obj.Nu_S=0;
            obj.obj=struct;
            obj.par_num=4;
            obj.depth=1;
            obj.bif_type={'bpt'};
            obj.point_list=10000;
            %MY Construct an instance of this class
            
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

