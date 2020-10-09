function out = transferLM(nose,refnose,LM)
         out = clone(LM);
         bar = cart2bary(refnose.Vertices,refnose.Faces,LM.Vertices,LM.Fi);
         out.Vertices = bary2cart(nose.Vertices,nose.Faces,bar,LM.Fi);
end