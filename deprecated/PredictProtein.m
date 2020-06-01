%Comparison of Boundary position (Xhalf, PXhalf(for different parameters), and NBXhalf) ,
%and Width (Width, PWidth(for diff. parameters), and NBWidth)

%Function that has matrix of Dp and ramma P as input, and gives output as
%matrix(rows: for different time points) of PXhalf and PWidth each.
function [PXhalf,PWidth] = PredictProtein(Pdiffusion,Phalf,NewTime,InterpSmoothmRNA)
    clear protein
    %clear PXhalf
    %PDiffusion and GammaP are matrix which have different parameters.
    %For specific values of PDiffusion and GammaP, I need to simulate
    %protein from AccumulatedmRNA (Interpolated)
    
    AP=[1:41];
    APbin=[0:0.025:1];
    dx=500/40; %delta x = 500um/40(bins) = 12.5um/bin
    APbinpos = (AP-1)*dx;
    
    deltaT=0.01; %min (I can change it later)
    protein=zeros(length(InterpSmoothmRNA),length(AP));
    protein(1,:)=0;  %protein(t=0)=0 at all APbins. Initial Condition
    
    
    %r_p:Translation rate for each AP bin, at nc14. 
    rp=2; %(2 proteins/mRNA, min), I can change it later. from Petkova et al., 2014 Curr. Bio (Gregor lab)
    
    %First, I need to set Dp and gammaP 's units correctly
    Dp=Pdiffusion*60; %um^2/min
    kp=Dp/(dx)^2; %1/Min. unit.
    
    gammaP=log(2)/Phalf; %Use Min. as time unit. I can change this later
    
    for k=2:length(NewTime)
    %Make a Protein reaction-diffusion equation
        protein(k,10)=protein(k-1,10)+rp*InterpSmoothmRNA(k-1,10)*deltaT-gammaP*protein(k-1,10)*deltaT; %+kp*protein(i-1,11)*dt-kp*protein(i-1,10)*dt;
        protein(k,19)=protein(k-1,19)+rp*InterpSmoothmRNA(k-1,19)*deltaT-gammaP*protein(k-1,19)*deltaT; %+kp*protein(i-1,25)*dt-kp*protein(i-1,41)*dt;
        
        for m=11:18%length(AP)-1 %for all APbins
            protein(k,m)=protein(k-1,m)+rp*InterpSmoothmRNA(k-1,m)*deltaT-gammaP*protein(k-1,m)*deltaT+kp*deltaT*(protein(k-1,m-1)+protein(k-1,m+1))-2*kp*protein(k-1,m)*deltaT;
        end
    end
    
    
    %Find the x position where it has half of the Maximum intensity
    for t=1:length(protein)
        MaxIntensity(t)=max(protein(t,:));
    end
    
    for tpoint=3:length(NewTime)
        %pchip fitting
        dxx=0.001;
        xx=0.2:dxx:0.6;
        yy=pchip(APbin,protein(tpoint,:),xx);
    
%         xxx=0.2:dxx:0.6;
%         yyy=pchip(APbin,protein(tpoint,:),xx);


        for i=2:length(xx)
            if (yy(i)<0.5*MaxIntensity(tpoint))%&(yy(i)~=0)
                if yy(i-1)>0.5*MaxIntensity(tpoint)
                    PXhalf(tpoint)=xx(i);
                    HalfMaxY(tpoint)=yy(i);
                    Index(tpoint)=i;
                end  
            end
        end
        %PXhalf(tpoint)
        
        %Get the Width using xx,yy pchip fit
        %First, draw a tangent line at the midpoint.
        %Mid point
        MidX(tpoint)=PXhalf(tpoint);
        MidY(tpoint)=HalfMaxY(tpoint);
    
        %tangent line
        dX=0.001;
        X=0:dX:1;
        grad=diff(yy)./diff(xx);
        %grad(Index) : slope(gradient) at the midpoint.
        MidIndex=Index(tpoint);
        Y=grad(MidIndex)*(X-MidX(tpoint))+MidY(tpoint);
    
        Maximum=max(yy);
        Minimum=min(yy);
    
        for i=1:length(X)
            if Y(i)>Maximum
                if Y(i+1)<Maximum
                    Left=i;
            %else
                %print('sth went wrong')
                end
            end
        end
    
        for j=1:length(X)
            if Y(j)<Minimum
                if Y(j-1)>Minimum
                    Right=j;
                %else
                    %print('sth went wrong')
                end
            end
        end
            
        PWidth(tpoint) = (Right-Left)*dX;
        
    end
    
end