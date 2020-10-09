function out = assessIdentity(in)
         out = in;
         out.Assessment = compareMorphs(out.AmplFace,out.BaseFace);
end