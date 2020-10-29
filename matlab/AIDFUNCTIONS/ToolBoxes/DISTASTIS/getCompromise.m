function [compromise] = getCompromise(S,W)
    compromise=zeros(size(S,1),size(S,1));
    for k=1:size(S,3);
        compromise=compromise+W(k)*S(:,:,k);
    end
end