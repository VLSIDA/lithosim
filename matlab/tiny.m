Resolution = 15;    %nm per pixel in the input image
NA = 0.95;          %numerical aperture
lambda = 193;       %wavelength
threshold = 0.2;    %Resist Threshold
Filter_Size = 310;   %Jinc filter size

I = double(imread('../tests/tiny.pbm'));

H = make_jinc( Filter_Size, (NA/lambda)*Resolution );

Field = imfilter ( I, H, 'replicate', 'same' ) ;
Aerial = abs(Field).^2;

Contour = Aerial > threshold;

f1=figure; imagesc(I); axis image; colormap gray; title('Original');
f2=figure; imagesc(Aerial); axis image; title('Aerial Image');
f3=figure; imagesc(Contour); axis image; colormap gray; title('Contours');
%saveas(f1,'tiny-original.png');
%saveas(f2,'tiny-aerial.png');
%saveas(f3,'tiny-contours.png');
imwrite(I,'tiny-original.tiff','Compression','none')
imwrite(Aerial,'tiny-aerial.tiff','Compression','none')
imwrite(Contour,'tiny-contours.tiff','Compression','none')

