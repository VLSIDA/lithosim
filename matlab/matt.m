Resolution = 10;    %nm per pixel in the input image
NA = 0.95;          %numerical aperture
lambda = 193;       %wavelength
threshold = 0.2;    %Resist Threshold
Filter_Size = 100;   %Jinc filter size

I = double(imread('../tests/2_thick_bar7.pbm'));

H = make_jinc( Filter_Size, (NA/lambda)*Resolution );

Field = imfilter ( I, H, 'replicate', 'same' ) ;
Aerial = abs(Field).^2;

Contour = Aerial > threshold;

f1=figure; imagesc(I); axis image; colormap gray; title('Original');
f2=figure; imagesc(Aerial); axis image; title('Aerial Image');
f3=figure; imagesc(Contour); axis image; colormap gray; title('Contours');
imwrite(I,'original.tiff','Compression','none')
imwrite(Aerial,'aerial.tiff','Compression','none')
imwrite(Contour,'contours.tiff','Compression','none')
