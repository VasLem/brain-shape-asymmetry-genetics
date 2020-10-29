function out = amplifyIdentity(in,AF)
out = in;
dir = out.PredFace.Vertices-out.BaseFace.Vertices;
out.AmplFace = clone(out.PredFace);
out.AmplFace.Vertices = out.BaseFace.Vertices + AF*dir;
end