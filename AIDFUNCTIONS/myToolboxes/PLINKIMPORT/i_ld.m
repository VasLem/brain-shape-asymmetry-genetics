function [d_raw,d_prime,r2,r,p,H] = i_ld(site1,site2)
	n=length(site1);

	stAa=unique(site1);
	stBb=unique(site2);

	Aid=1; Bid=1;
	fA=sum(site1==stAa(Aid))/n;
	fB=sum(site2==stBb(Bid))/n;

	%D_raw's sign is arbitrary:
	%A common convention is to set A, B to be the
	%common allele and a, b to be the rare allele

	if (fA<0.5), fA=1-fA; Aid=2; end  % Aid tells which one is common allele
	if (fB<0.5), fB=1-fB; Bid=2; end  % Bid tells which one is common allele

		x=0;
		for k=1:n
		if (site1(k)==stAa(Aid) && site2(k)==stBb(Bid)),
			x=x+1;
		end
		end
		d_raw=x./n-fA*fB;

	if (nargout>1),
		fa=1-fA; fb=1-fB;

        %r2=(d_raw*d_raw)./(prod(probMajor)*prod(1-probMajor));
        %r2=min([1 r2]);

		r2=(d_raw*d_raw)./(fA*fa*fB*fb);
		r=sqrt(r2);

		if (d_raw<0),
		      x=min(fA*fB,fa*fb);
		      r=-1*r;
		else
		      x=min(fA*fb,fa*fB);
		end
		d_prime=abs(d_raw)./x;
	end


if nargout>4
	 	M=zeros(2);
	 	for k=1:n
	 		if (site1(k)==stAa(1) && site2(k)==stBb(1)),
	 		      M(1,1)=M(1,1)+1;
	 		elseif (site1(k)==stAa(1) && site2(k)==stBb(2)),
	 		      M(1,2)=M(1,2)+1;
	 		elseif (site1(k)==stAa(2) && site2(k)==stBb(1)),
	 		      M(2,1)=M(2,1)+1;
	 		elseif (site1(k)==stAa(2) && site2(k)==stBb(2)),
	 		      M(2,2)=M(2,2)+1;
	 		end
	 	end
		%[p]=fisherextest(M(1,1),M(1,2),M(2,1),M(2,2));
        [H,p] = fishertest(M);
end