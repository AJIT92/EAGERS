function A = ComponentProperty(varargin)
global modelParam Tags
property = varargin{1};
if length(varargin) ==1
    if ischar(property)
        r = strfind(property,'.');
        if isempty(r)
            if isfield(modelParam,property)
                A = modelParam.(property);
            else A = strcat('no field:',property);
            end
        else
            if strcmp(property(1:r(1)-1),'Tags')
                if length(r) ==1
                    if isfield(Tags,property(1:r(1)-1)) && isfield(Tags.(property(1:r(1)-1)),property(r(1)+1:end))
                        A = Tags.(property(1:r(1)-1)).(property(r(1)+1:end));
                    else A = strcat('no field:',property);
                    end
                elseif length(r) ==2
                    A = Tags.(property(1:r(1)-1)).(property(r(1)+1:r(2)-1)).(property(r(2)+1:end));
                elseif length(r) ==3
                    A = Tags.(property(1:r(1)-1)).(property(r(1)+1:r(2)-1)).(property(r(2)+1:r(3)-1)).(property(r(3)+1:end));
                end
            else
                if length(r) ==1
                    if isfield(modelParam,property(1:r(1)-1)) && isfield(modelParam.(property(1:r(1)-1)),property(r(1)+1:end))
                        A = modelParam.(property(1:r(1)-1)).(property(r(1)+1:end));
                    else A = strcat('no field:',property);
                    end
                elseif length(r) ==2
                    A = modelParam.(property(1:r(1)-1)).(property(r(1)+1:r(2)-1)).(property(r(2)+1:end));
                elseif length(r) ==3
                    A = modelParam.(property(1:r(1)-1)).(property(r(1)+1:r(2)-1)).(property(r(2)+1:r(3)-1)).(property(r(3)+1:end));
                end
            end
        end
    else A = property;
    end
else
    A = [];
    value = varargin{2};
    r = strfind(property,'.');
    if isempty(r)
        modelParam.(property) = value;
    elseif length(r) ==1
        modelParam.(property(1:r(1)-1)).(property(r(1)+1:end)) = value;
    elseif length(r) ==2
        modelParam.(property(1:r(1)-1)).(property(r(1)+1:r(2)-1)).(property(r(2)+1:end)) = value;
    elseif length(r) ==3
        modelParam.(property(1:r(1)-1)).(property(r(1)+1:r(2)-1)).(property(r(2)+1:r(3)-1)).(property(r(3)+1:end)) = value;
    end
end