function stack=mtifread(file,channel)

info = imfinfo(file);
num= numel(info);

stack=zeros(info(1).Height,info(1).Width,num/channel,channel);
k=1;
for j=1:num/channel
    for i=1:channel
        stack(:,:,j,i)=imread(file,k);
        k=k+1;
    end
end