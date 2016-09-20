% Flush out the MATLAB.
close all;

clc;

clear all;


% Input image.
InputImage = imread('building.pnm');


% Converting the input image and getting the height an width of the image.
InputImage = im2double(InputImage);

[H,W] = size(InputImage);

InputImage = InputImage(1:125,1:125);


% Display the original image.
figure,imshow(InputImage);

title(' Original Image : ');


% Get edges from the image.
EdgeImage = edge(InputImage,'canny',0.05,'both',1);


% Display the original image.
figure,imshow(EdgeImage);

title(' Edge Image : ');


% Getting the height an width of the image.
[Height,Width] = size(EdgeImage);


% Implementing the hough transform.
Max = ceil(.70710678*(Height+Width));

Min = -Max;

Parameter = 100;

Hough = zeros(2*Max+1,Parameter);

Temp = linspace(0,pi,Parameter+ 1);

for a = 1:Height,
    
    for b = 1:Width,
        
        if (EdgeImage(a,b) == 1)
            
            for c = 1:Parameter,
                
                Round = round(b*cos(Temp(c))+a*sin(Temp(c)));
                
                Temp2 = Round+Max+1;
                
                Hough(Temp2,c) = Hough(Temp2,c)+ 1;
                
            end
            
        end
        
    end
    
end


% Detecting the lines.
Threshold = 50;

Line = zeros(Height,Width);

for Temp2 = 1:2*Max+1,
    
    for c = 1:Parameter,
        
        if (Hough(Temp2 ,c) > Threshold)
            
            for a = 1:Height,
                
                for b = 1:Width,
                    
                    Round = Temp2-Max-1;
                    
                    if (Round == round(b*cos(Temp(c)) + a *sin(Temp(c))))
                        
                        Line(a,b) = 1;
                        
                    end
                    
                end
                
            end
            
        end
        
    end
    
end


% Display the original image.
figure,imshow(Line);

title(' Final Output : ');