function ImageArray = openxlsFigures(FilePath, FileName, pos, head, tail)

impPath = [FilePath];
expPath = [FilePath 'Merge\'];

loadexist = 1;
if(loadexist && exist([expPath FileName pos '.mat'], 'file') )
    load([expPath FileName pos '.mat'], 'ImageArray');
    return;
end

ImageArray = 0;

if(nargin == 3)
    head = 0;
    tail =head;
end

for ii = head : tail
    if(nargin == 3)
        openFile = [impPath FileName '.xls'];
    else
        openFile = [impPath FileName num2str(ii) '.xls'];
    end
    RawImage = fopen(openFile, 'r');
    ImageArray = ImageArray + fread(RawImage, [1392 2080], 'uint16', 'ieee-be')';
    %     ImageArray = ImageArray + (fread(RawImage, [1392 2080], 'uint16', 'ieee-be')') .^ (-2);
    fclose(RawImage);
end
    
if(strcmp(pos, 'up') )
    ImageArray = ImageArray(1 : 1040, :);
else if(strcmp(pos, 'down') )
        ImageArray = ImageArray(1041 : 2080, :);
    else
        error('Did not specify the position');
    end
end

ImageArray =  ImageArray / (tail - head +1);
save([expPath FileName pos '.mat'], 'ImageArray');

% ImageArray = sqrt( (tail - head +1) ./ ImageArray);

%figure;imagesc(ImageArray);
%a11=ImageArray(1:1040,:);
%a12=ImageArray(1041:2080,:);