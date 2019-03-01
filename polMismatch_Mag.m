function [z,Pt]= polMismatch_Mag(M,R,tag_identical,tag_isot,delta,E_GS,E_UAV,lambda,N0,B);

for e = 1:1000
    e
    XX = zeros(1,M);
    
%     alphaxk = randsrc(1,1,[-180:.01:180]);%300;
%     alphayk = randsrc(1,1,[-180:.01:180]);%50;
%     alphazk = randsrc(1,1,[0:.01:90]);%150;
    
    alphaxk = rand*360;%randsrc(1,1,[-180:.01:180]);%300;
    alphayk = rand*360;%randsrc(1,1,[-180:.01:180]);%50;
    alphazk = rand*360;%randsrc(1,1,[0:.01:90]);%150;
    
    rvals = 2*rand(1,1)-1;
    el = asin(rvals)+pi/2;
    az = 2*pi*rand(1,1);
    radii = R*(rand(1,1)).^(1/3);
    
%     rvals = 2*rand(1,1)-1;
%     el = 90;asin(rvals)+pi/2;
%     az = 90;2*pi*rand(1,1);
%     radii = R*(rand(1,1)).^(1/3);
    
    if radii>20
        
        x_j = radii.*cos(az).*sin(el);
        y_j = radii.*sin(az).*sin(el);
        z_j = radii.*cos(el);
        dk = sqrt(x_j^2+y_j^2+z_j^2);
        
        
        
        thetak = el; %randsrc(1,1,[0:.1:pi]);
        phik = az; %randsrc(1,1,[0:.1:pi]);
        
        for l=1:M
            % l
            if tag_identical ==1
                alphaxl = 0;
                alphayl = 0;
                alphazl = 0;
            else
                alphaxl = rand*360;%randsrc(1,1,[0:.01:360]);(l/M)*360;
                alphayl = rand*360;%randsrc(1,1,[0:.01:360]);
                alphazl = rand*360;%randsrc(1,1,[0:.01:360]);
            end
            
            theta_array = [0:pi/200:pi];
            phi_array = [0:pi/200:2*pi];
            
            R3xk = rotx(alphaxk);
            R3yk = roty(alphayk);
            R3zk = rotz(alphazk);
            
            R3k = R3xk*R3yk*R3zk;%[1 0 0; 0 1 0; 0 0 1];
            
            R3xl = rotx(alphaxl);
            R3yl = roty(alphayl);
            R3zl = rotz(alphazl);
            
            R3l = R3xl*R3yl*R3zl;%[1 0 0; 0 1 0; 0 0 1];
            
            dkl = dk*(1+(1/(2*dk))*(((l-1)*delta)^2-((l-1)*delta*sin(thetak)*cos(phik))));
            
            thetakl(l) = acos(dk*cos(thetak)/dkl);
            phikl(l) = acot(((dk*sin(thetak)*cos(phik))-((l-1)*delta))/((dk*sin(thetak)*sin(phik))-(0)))-(pi/2)*(sign(dk*sin(thetak)*sin(phik)-0)-1);
            
            thetakl_1 = acos((-sin(thetakl(l))*cos(phikl(l))*R3k(3,1))+(-sin(thetakl(l))*sin(phikl(l))*R3k(3,2))+(-cos(thetakl(l))*R3k(3,3)));
            psikl_1 = acos((-sin(thetakl(l))*cos(phikl(l))*R3k(2,1))+(-sin(thetakl(l))*sin(phikl(l))*R3k(2,2))+(-cos(thetakl(l))*R3k(2,3)));
            
            thetakl_2 = acos((sin(thetakl(l))*cos(phikl(l))*R3l(3,1))+(sin(thetakl(l))*sin(phikl(l))*R3l(3,2))+(cos(thetakl(l))*R3l(3,3)));
            psikl_2 = acos((sin(thetakl(l))*cos(phikl(l))*R3l(2,1))+(sin(thetakl(l))*sin(phikl(l))*R3l(2,2))+(cos(thetakl(l))*R3l(2,3)));
            
            
            F_theta_1 = cos((pi/2)*cos(thetakl_1))/sin(thetakl_1);
            F_psi_1 = cos((pi/2)*cos(psikl_1))/sin(psikl_1);
            
            F_theta_2 = cos((pi/2)*cos(thetakl_2))/sin(thetakl_2);
            F_psi_2 = cos((pi/2)*cos(psikl_2))/sin(psikl_2);
            
            if tag_isot == 0
                F_GS = [F_theta_1 F_psi_1].';
                F_UAV = [F_theta_2 F_psi_2].';
            else
                F_GS = [1 1].';
                F_UAV = [1 1].';
            end
            
            E_GS_R = E_GS.*F_GS;
            E_UAV_R = E_UAV.*F_GS;
            
            
            cosa1 = inv(sin(thetakl_2))*((R3l(3,1)-cos(thetakl_2)*sin(thetakl(l))*cos(phikl(l)))*R3k(3,1)+(R3l(3,2)-cos(thetakl_2)*sin(thetakl(l))*sin(phikl(l)))*R3k(3,2)+(R3l(3,3)-cos(thetakl_2)*cos(thetakl(l)))*R3k(3,3));
            cosa2 = inv(sin(psikl_2))*((R3l(2,1)-cos(psikl_2)*sin(thetakl(l))*cos(phikl(l)))*R3k(3,1)+(R3l(2,2)-cos(psikl_2)*sin(thetakl(l))*sin(phikl(l)))*R3k(3,2)+(R3l(2,3)-cos(psikl_2)*cos(thetakl(l)))*R3k(3,3));
            cosb1 = inv(sin(thetakl_2))*((R3l(3,1)-cos(thetakl_2)*sin(thetakl(l))*cos(phikl(l)))*R3k(2,1)+(R3l(3,2)-cos(thetakl_2)*sin(thetakl(l))*sin(phikl(l)))*R3k(2,2)+(R3l(3,3)-cos(thetakl_2)*cos(thetakl(l)))*R3k(2,3));
            cosb2 = inv(sin(psikl_2))*((R3l(2,1)-cos(psikl_2)*sin(thetakl(l))*cos(phikl(l)))*R3k(2,1)+(R3l(2,2)-cos(psikl_2)*sin(thetakl(l))*sin(phikl(l)))*R3k(2,2)+(R3l(2,3)-cos(psikl_2)*cos(thetakl(l)))*R3k(2,3));
            
            
            R2 = [cosa1 cosa2; cosb1 cosb2];
            
            X = 1.64^2*abs((E_GS_R.'*R2*E_UAV_R))^2;
            if X>.000001
                Y = inv(X);
            else
                Y = NaN;1000;
            end
            
            XX(l) = X;
        end
        z(e) =  1/mean(XX);  %1/sum(XX(~isnan(XX)));        

    end
    Pt=N0*B*(4*pi*dk/lambda)^2*z;
                    XXArray(e,:) = XX;
end
v=1;
