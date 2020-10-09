% testing estimation of properties
% load in the data
close all;clear all;clear classes;
load('C:\MATLAB\Work\Projects\PropertyManipulator\AngleShape\AngleData.mat');
load('C:\MATLAB\Work\Projects\PropertyManipulator\AngleShape\RefFace\RefFace.mat');
RefScan = RefFace;
obj = shapePCA;
obj.RefScan = RefScan;
getAverage(obj,Data);
Dshape = alignShape(obj,Data);
getModel(obj,Dshape);
ShapeModel = obj;
stripPercVar(ShapeModel,98);
%%
PropModel = propertyPCA;
PropModel.Model = ShapeModel;
getModel(PropModel,Data);
%%
FM = FaceMorpher;
FM.Model = PropModel;
%%
out = estimatePropFromShapeCoeff(PropModel,ShapeModel.Tcoeff(1,:)');

%% Validate prop shape covariance
 prop = (40);
 trueprop = nan*zeros(1,PropModel.n);
 estprop = nan*zeros(1,PropModel.n);
 indexall = (1:PropModel.n);
 f = statusbar('Progress');
 for i=1:1:PropModel.n
     tmp = Coeff2Struc(PropModel,PropModel.Tcoeff(i,:)');
     trueprop(i) = tmp.Prop(prop);
     clear tmp;
     tmpShapeModel = clone(ShapeModel);
     tmpShapeModel.Tcoeff = tmpShapeModel.Tcoeff(setdiff(indexall,i),:);
     tmpData = reduceData(Data,setdiff(indexall,i));
     tmpPropModel = propertyPCA;
     tmpPropModel.Model = tmpShapeModel;
     getModel(tmpPropModel,tmpData);
     out = estimatePropFromShapeCoeff(PropModel,ShapeModel.Tcoeff(i,:)');
     estprop(i) = out(prop);
     clear out tmpData;
     delete(tmpPropModel);
     delete(tmpShapeModel);
     statusbar(i/PropModel.n,f);
 end
 delete(f);
 %%
 %[h,p] = kstest(trueprop-estprop);
 [h,p] = lillietest(trueprop-estprop,0.05);% P = 0.1754 P = 0.1214
 [H,P,CI] = ttest(trueprop,estprop);% P = 0.9993 P = 0.9945
 
%% ALL OK 
 
% H must be zero meaning no significant difference! for a P value of 0.05
