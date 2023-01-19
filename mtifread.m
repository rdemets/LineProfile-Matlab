function stack=mtifread(file)

info = imfinfo(file);
num= numel(info);

stack=zeros(info(1).Height,info(1).Width,num);

for i=1:num
    stack(:,:,i,:)=imread(file,i);
end