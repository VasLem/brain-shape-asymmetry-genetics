function out = selectSNP(in,index, display)
   if nargin<3, display = false; end
% function to make a selection
   out = in;
   if isfield(in,'POS'), out.POS = in.POS(index); end
   if isfield(in,'A1'),out.A1 = in.A1(index);end
   if isfield(in,'A2'),out.A2 = in.A2(index);end
   if isfield(in,'RSID'),out.RSID = in.RSID(index);end
   if isfield(in,'GCount'),out.GCount = in.GCount(index);end
   if isfield(in,'Ohet'),out.Ohet = in.Ohet(index);end
   if isfield(in,'Ehet'),out.Ehet = in.Ehet(index);end
   if isfield(in,'P'),out.P = in.P(index);end
   if isfield(in,'MAFplink'),out.MAFplink = in.MAFplink(index);end
   if isfield(in,'NCHROBS'),out.NCHROBS = in.NCHROBS(index);end
   if isfield(in,'FMISsnp'),out.FMISsnp = in.FMISsnp(index);end
   out.nSNP = length(out.POS);
   if display, disp(['Selected ' num2str(out.nSNP) ' from ' num2str(in.nSNP) ' SNPs']);end
end