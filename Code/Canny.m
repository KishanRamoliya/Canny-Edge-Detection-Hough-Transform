% Flush out the MATLAB.
close all;

clc;

clear all;


% Input image.
InputImage = imread('pillsetc.pnm'); 


% Converting the input image and getting the height an width of the image.
InputImage = im2double(InputImage);    

[H,W] = size(InputImage);


% Take input of the required pameters. 
Sigma = input(' Enter the Sigma : ');

MinimumThreshold = input(' Enter the Minimum Threshold : ');

MaximumThreshold = input(' Enter the Maximum Threshold : ');


% Getting the derivative in X direction.
XDerivative = zeros(H,W);


% Getting the derivative in Y direction.
YDerivative = zeros(H,W);  
 
 
% Define the gaussian kernal on the basis of X and Y direction.
Kernal = (6*Sigma)+1;  

XDirection = zeros(Kernal,Kernal);

YDirection = zeros(Kernal,Kernal);

for a = 1:Kernal
    
    for b = 1:Kernal
        
        YDirection(a,b) = -((a-((Kernal-1)/2)-1)/(2*pi*Sigma^3))*exp(-((a-((Kernal-1)/2)-1)^2+(b-((Kernal-1)/2)-1)^2)/(2*Sigma^2));
         
    end
    
end
 
for a = 1:Kernal
    
    for b = 1:Kernal
        
        XDirection(a,b) = -((b-((Kernal-1)/2)-1)/(2*pi*Sigma^3))*exp(-((a-((Kernal-1)/2)-1)^2+(b-((Kernal-1)/2)-1)^2)/(2*Sigma^2));
        
    end
    
end

G = zeros(H,W);

NonmaximaSupression = zeros(H,W);     
 
Hysteresis = zeros(H,W); 
 
 
% Calculate the derivative wrt. X direction of image. 
for i = 1+ceil(Kernal/2):H-ceil(Kernal/2)
    
    for j = 1+ceil(Kernal/2):W-ceil(Kernal/2)  
        
        Row =  i-ceil(Kernal/2); 
        
        Column =  j-ceil(Kernal/2); 
        
        for K = 1:Kernal  
            
            for l = 1:Kernal  
                
                XDerivative(i,j) = XDerivative(i,j) + InputImage(Row+K-1, Column+l-1)*XDirection(K,l);
                
            end
            
        end
        
    end
    
end
 

% Calculate the derivative wrt. Y direction of image.
for i = 1+ceil(Kernal/2):H-ceil(Kernal/2)
    
    for j = 1+ceil(Kernal/2):W-ceil(Kernal/2) 
        
        Row = i-ceil(Kernal/2); 
        
        Column = j-ceil(Kernal/2); 
        
        for K = 1:Kernal  
            
            for l = 1:Kernal 
                
                YDerivative(i,j) = YDerivative(i,j) + InputImage(Row+K-1, Column+l-1)*YDirection(K,l);
                
            end
            
        end
        
    end
    
end
 
 
% Calculate magnitude using the above derivatives. 
for i = 1+ceil(Kernal/2):H-ceil(Kernal/2)  
    
    for j = 1+ceil(Kernal/2):W-ceil(Kernal/2)  
        
        G(i,j) = sqrt(XDerivative(i,j)^2 + YDerivative(i,j)^2);
        
    end
    
end
 

% Implement non maxima supression. 
NonmaximaSupression = G;

for i = 1+ceil(Kernal/2):H-ceil(Kernal/2) 
    
    for j = 1+ceil(Kernal/2):W-ceil(Kernal/2)  
       
        if (XDerivative(i,j) == 0)
            
            T = 5;
            
        else
            
            T = (YDerivative(i,j)/XDerivative(i,j));
            
        end
        
        if (-0.4142<T & T<=0.4142)
            
            if(G(i,j)<G(i,j+1) | G(i,j)<G(i,j-1))
                
                NonmaximaSupression(i,j)=0;
                
            end
            
        end
        
        if (0.4142<T & T<=2.4142)
            
            if(G(i,j)<G(i-1,j+1) | G(i,j)<G(i+1,j-1))
                
                NonmaximaSupression(i,j)=0;
                
            end
            
        end
        
        if (abs(T) >2.4142)
            
            if(G(i,j)<G(i-1,j) | G(i,j)<G(i+1,j))
                
                NonmaximaSupression(i,j)=0;
                
            end
            
        end
        
        if (-2.4142<T & T<= -0.4142)
            
            if(G(i,j)<G(i-1,j-1) | G(i,j)<G(i+1,j+1))
                
                NonmaximaSupression(i,j)=0;
                
            end
            
        end
        
    end
    
end
 

% Implement hysteresis. 
Hysteresis = NonmaximaSupression;
 
for i = 1+ceil(Kernal/2):H-ceil(Kernal/2) 
    
    for j = 1+ceil(Kernal/2):W-ceil(Kernal/2)  
        
        if(Hysteresis(i,j) >= MaximumThreshold)
            
            Hysteresis(i,j) = 1;
            
        end
        
        if(Hysteresis(i,j)<MaximumThreshold & Hysteresis(i,j)>=MinimumThreshold)
            
            Hysteresis(i,j) = 2;
            
        end
        
        if(Hysteresis(i,j)<MinimumThreshold)
            
            Hysteresis(i,j) = 0;
        
        end 
        
    end
    
end
 
Temp = 1; 
 
while (Temp == 1)
   
    Temp = 0;
   
    for i = 1+ceil(Kernal/2):H-ceil(Kernal/2)  
        
        for j = 1+ceil(Kernal/2):W-ceil(Kernal/2)  
            
            if (Hysteresis(i,j) > 0)      
                
                if(Hysteresis(i,j) == 2)    
                   
                    if( Hysteresis(i-1,j-1) == 1 | Hysteresis(i-1,j) == 1 | Hysteresis(i-1,j+1) == 1 | Hysteresis(i,j-1) == 1 |  Hysteresis(i,j+1) == 1 | Hysteresis(i+1,j-1) == 1 | Hysteresis(i+1,j) == 1 | Hysteresis(i+1,j+1) == 1 )
                        
                        Hysteresis(i,j) = 1;
                        
                        Temp == 1;
                        
                    end
                    
                end
                
            end
            
        end
        
    end
   
end
 
for i = 1+ceil(Kernal/2):H-ceil(Kernal/2) 
    
    for j = 1+ceil(Kernal/2):W-ceil(Kernal/2)  
        
        if(Hysteresis(i,j) == 2) 
            
            Hysteresis(i,j) == 0;
            
        end 
        
    end
    
end
 
 
% Outputs.
imwrite(InputImage,'OriginalImage.bmp');

imwrite(XDerivative,'XDerivative.bmp');            

imwrite(YDerivative,'YDerivative.bmp');            

imwrite(G,'Gradient.bmp');       

imwrite(NonmaximaSupression,'NonmaximaSupression.bmp');      

imwrite(Hysteresis,'FinalImage.bmp');