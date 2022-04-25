classdef mySuperClass < hgsetget & dynamicprops
%     This is the abstract interface class for all Classes
%     it contains routines that apply to all objects that inherit from this
%     class
%     The HGSETGET class is an abstract class that provides an HG-style
%     property set and get interface.  HGSETGET is a subclass of HANDLE, so 
%     any classes derived from HGSETGET are handle classes.
%     (set,get,setdisp,getdisp)
%     The DYNAMICPROPS class is an abstract class that provides support
%     for dynamic properties for MATLAB objects.  Dynamic properties can
%     be used to attach temporary data to MATLAB objects.
%     (addprop)
    properties (Dependent = true)
        Type;% This Property enables to detect an object class when this latter has been converted to a structure;
    end
    methods %Constructor
        function obj = mySuperClass(varargin)
        end
    end
    methods% Special Setting & Getting
        function type = get.Type(obj)
            type = class(obj);
        end
        function obj = set.Type(obj,type) %#ok<INUSD>
            return;%dummy
        end
    end
    methods% Interface Functions
        function delete(obj) %#ok<INUSD>
        end
        function cobj = clone(obj)
            cobj = eval(class(obj));
            copy(obj,cobj);
        end
        function out = validH(obj,handle)
            % out = validHandle(obj,handle)
            % Property Field to test (handle) must be given as a String
            H = obj.(handle);% store in variable to eliminate risk for recursions
            out = mySuperClass.isH(H);% isH is a function defined in objectclass Utilities            
        end
        function struc = saveobj(obj)
            struc = obj2struc(obj);
        end
        function copy(obj,cobj)
            if~strcmp(obj.Type,cobj.Type), error('Copy requires objects of the same class'); end
        end
        function struc = obj2struc(obj)%needed to save objects
            struc.Type = obj.Type;
        end
    end
    methods (Abstract = true)
        obj = struc2obj(obj,struc);% needed to load objectss       
    end
    methods (Static = true)
        function obj = loadobj(struc)
            obj = struc2obj(eval(struc.Type),struc);
        end
        function out = isH(in)
            % test whether in is a valid handle
            if isempty(in), out = false; return; end% empty check
            if ~isa(in,'handle'), out = false; return; end% handle check
            if ~isvalid(in), out = false; return; end% check whether handle has been deleted 
            out = true;% succeeded all tests
        end
    end
end % classdef