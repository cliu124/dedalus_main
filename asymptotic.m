classdef asymptotic
    %ASYMPTOTIC Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        indepedent;%these are space or time independent variable.
        variable;%These are dependent variable
        parameter;%This is the parameter that are not related to dependent variable
        equation; 
        
        epsilon;
        
        
    end
    
    methods
        function obj = asymptotic(inputArg1,inputArg2)
            %ASYMPTOTIC Construct an instance of this class
            %   Detailed explanation goes here
            obj.Property1 = inputArg1 + inputArg2;
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

