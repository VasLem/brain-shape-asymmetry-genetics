function Path = FindLicense

global FastRBF_LICENSE

Path = FastRBF_LICENSE;

if isempty(Path)
  name = 'FastRBF_MEX';
  len = length(name);
  Path = which( 'FindLicense' );
  Path = [Path(1:end-len-2), 'FastRBF.lic'];
end
