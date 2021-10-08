% settings 
fno=14;   % this is the fault number to read in & used for file name 
L_km=15;  % this is the half length of the fault normal profile centered on trace
S_km=1;   % this is the spacing between fault normal profiles 

%% read in data for a shape file having only on fault trace 
S=shaperead('u14_fault_trace.shp');
lat=S.Y; lat=lat(1:end-1); % name S.Y lat and take all but the last value
lon=S.X; lon=lon(1:end-1); % name S.X lon and take all but the last value

%%  read the MOLA elevation data 
%https://astrogeology.usgs.gov/search/map/Mars/Topography/HRSC_MOLA_Blend/Mars_HRSC_MOLA_BlendDEM_Global_200mp_v2
[Z, R] = readgeoraster('Mars_HRSC_MOLA_BlendDEM_Global_200mp_v2.tif'); 

% trim the DEM to area of interest for speed 
tmp1=range(lat); tmp2=range(lon)/2; 
LATLIM=[min(lat)-km2deg(5+L_km,3396.190), max(lat)+km2deg(5+L_km,3396.190)]; 
LONLIM=[min(lon)-km2deg(5+L_km,3396.190), max(lon)+km2deg(5+L_km,3396.190)]; 
[ZT, RT] = geocrop(Z,R,LATLIM,LONLIM); ZT=double(ZT);
clims=quantile(ZT(:),[0.25,0.75]);    % set color limits 
clear Z R;

%% define positions along fault to profile 
 %calculate the cumulative distance along fault 
         tmpd=nan(1,length(lat)-1); 
         for d=1:length(lat)-1 
             tmpd(d)=distance(lat(d),lon(d),lat(d+1),lon(d+1),[3396.190 0]); 
         end
         tmpd=cat(2,0,tmpd); km_along_fault=cumsum(tmpd);

         ptsS=0:S_km:km_along_fault(end); % distance to points to sample along trace 
         ptlat=interp1(km_along_fault,lat,ptsS); ptlat=cat(2,ptlat,lat(end)); % latitude of points to sample along trace 
         ptlon=interp1(km_along_fault,lon,ptsS); ptlon=cat(2,ptlon,lon(end)); % longitude of points to sample along trace
         
         % now recalculate distance for position of sample points  
         tmpd=nan(1,length(ptlat)-1); 
         for d=1:length(ptlat)-1 
             tmpd(d)=distance(ptlat(d),ptlon(d),ptlat(d+1),ptlon(d+1),[3396.190 0]);
         end
         tmpd=cat(2,0,tmpd); km_along_fault=cumsum(tmpd);


% azimuth of the fault at every point 
    AZ=nan(size(ptlat)); 
    for a=2:length(ptlat)-1  
        AZ(a)=azimuth(ptlat(a-1),ptlon(a-1),ptlat(a+1),ptlon(a+1),[3396.190 0],'degrees'); 
    end
% smooth the azimuth data 
    smfac=floor(length(ptlat)/20); if rem(smfac,2)==0; smfac=smfac+1; end   % make sure its an odd number 
    AZ(1)=AZ(2);  AZ(end)=AZ(end-1);  AZ=smoothazimuth(AZ,smfac); 

    
%% Now extract profile and plot in GUI for selecting top and bottom of scarp 
dataout=nan(length(ptlat),11);   % eleven parameters exported along fault normal profiles. 

%%
figure('Position', [100, 100, 1450, 800]) 
i=1; % start with the first profile 
while i <= length(ptlat) % use a while loop, to make it possible to repeat a pick 

    % map 1 
    subplot(2,1,1); hold off; axesm('MapProjection','Robinson','Geoid', [3396190 0]); h0=meshm(ZT, RT); caxis(clims);   % colored topography  
    h1=plotm(lat,lon,'k','LineWidth',0.5); h2=plotm(ptlat(i),ptlon(i),'.m','MarkerSize',12); hold on;  % trace of fault & current profile 
    h3=plotm(dataout(:,6),dataout(:,7),'w','LineWidth',1);  h4=plotm(dataout(:,9),dataout(:,10),'w','LineWidth',1); 
    title(['Fault #: ' num2str(fno) '  Profile #' num2str(i) '']); 
 
    % calculate the profile and display on map 
     tmp=track1([ptlat(i),ptlat(i)]',[ptlon(i),ptlon(i)]',[wrapTo360(AZ(i)+90),wrapTo360(AZ(i)-90)]',[L_km,L_km]',[3396.190 0],'degrees',1); % find end points 
     tmp=sortrows([tmp(1),tmp(3); tmp(2), tmp(4)],1); % sort them so that the most southern point is first 
     [R_z,R_r,R_lat,R_lon] = mapprofile(ZT,RT,tmp(:,1),tmp(:,2),[3396.190 0],'gc','nearest');   % this will extract the profile between the end pints 
     subplot(2,1,1); hold on; h5=plotm(R_lat,R_lon,'k','LineWidth',1); h6=plotm(R_lat(1),R_lon(1),'ro','MarkerSize',8);  
     
   
     
      % plot the range verses elevation data for the profile and add the
      % local slope 
     subplot(2,1,2); plot(R_r,R_z); hold on; plot(R_r(1),R_z(1),'ro','MarkerSize',12); % profile data 
     scatter(R_r(2:end),R_z(2:end),8,smooth(atand(diff(R_z/1000)./(diff(R_r))),3)); c=colorbar;  c.Label.String='Slope Deg'; c.Label.FontSize=18;    % slope 
     plot(ones(1,2).*range(R_r)/2, [min(R_z),max(R_z)],'m','LineWidth',1);   % draw a line where the trace was picked. 
     grid on; xlim([0,2*L_km]);  % set the plot limits. 
 
     % pick the points and plot them 
     title('CLICK the base and then CLICK the top of scarp and then hit ENTER','FontSize',20)  % give instructions 
     datain=ginput; %  get the point data 
     plot(datain(1,1),datain(1,2),'k+','MarkerSize',12,'LineWidth',2); text(datain(1,1)+2,datain(1,2),'B','FontSize',18);  % plot the bottom point on profile 
     plot(datain(2,1),datain(2,2),'k+','MarkerSize',12,'LineWidth',2); text(datain(2,1)+2,datain(2,2),'T','FontSize',18); hold off  % plot the top point on profile 
     mkB=find(R_r < datain(1,1),1,'last');   mkT=find(R_r > datain(2,1),1,'first');  % get there lat and lons 
     subplot(2,1,1); h7=plotm(R_lat(mkB),R_lon(mkB),'w.','MarkerSize',12); h8=plotm(R_lat(mkT),R_lon(mkT),'w.','MarkerSize',12);  % plot them on the map 
     clims=[min(R_z),max(R_z)];
     
     % now confirm you choices 
     subplot(2,1,2);
     title('Enter 1 to Accept, 2 to repeat, or  0 to skip this profile','FontSize',20)
     resp=input('Enter 1 to Accept, 2 to repeat, or  0 to skip this profile:    '); 

  % now do different things depending on if you want accept, repeat or skip
  % the profile 
    if resp==1
       disp('saving this profile') 
       throw=R_z(mkT)-R_z(mkB); 
       slopeangle=atand((R_z(mkT)-R_z(mkB))/distance(R_lat(mkB), R_lon(mkB), R_lat(mkT), R_lon(mkT),[3396190 0]));  % inverse tan of elevation in m over distance in meters 
       tmpout=[i, km_along_fault(i), AZ(i), throw, slopeangle, R_lat(mkB), R_lon(mkB), R_z(mkB), R_lat(mkT), R_lon(mkT), R_z(mkT)]; 
       dataout(i,:)=tmpout; 
       delete([h0,h1,h2,h3,h4,h5,h6,h7,h8])
       i=i+1; % can use -1 if you want to go the opposite way along the fault
    elseif resp==0
        i=i+1; % can use -1 if you want to go the opposite way along the fault
        disp('not recoding a throw for this profile') 
        delete([h0,h1,h2,h3,h4,h5,h6,h7,h8])
    else 
        disp('repeating the profile') 
        delete([h0,h1,h2,h3,h4,h5,h6,h7,h8])
    end
end


%% now write out the data into a text file 
    fno_s=sprintf('%03.0f',fno);
    fileout=strcat('fault_no', fno_s, '_v', datestr(now,'ddmmyyyy'), '_', datestr(now,'hhmm'), '.csv'); 
    fid=fopen(fileout,'w'); 
    fprintf(fid, 'profileNo, km_along_fault, fault_azimuth, throw_m, slope_deg, lat_base, lon_base, elev_base, lat_top, lon_top, elev_top\n');
    fprintf(fid,'%03.0f,%1.3f, %1.3f,%1.3f, %1.1f, %15.12f, %15.12f,%1.1f,%15.12f, %15.12f,%1.1f\n', dataout'); 
    fclose(fid); 


% General Comments 
% according to this website https://astrogeology.usgs.gov/search/map/Mars/Topography/HRSC_MOLA_Blend/Mars_HRSC_MOLA_BlendDEM_Global_200mp_v2
% the data assume a shereical mars with radius 3396190

%% plot the data throws vs dist
    sm = smooth(dataout(:,4));  % 5 point smoothing 
    figure; plot(dataout(:,2),dataout(:,4),'k',dataout(:,2),sm,'m')
    title('Throw vs. distance along fault u50','FontSize',20)  % give instructions 
    xlabel('Distance along the fault (km)');
    ylabel('Throw (m)');
    legend ('raw data', 'smooth moving avg 5');


%% Functions 

function [smAZ]=smoothazimuth(AZ,smfac)
%USE 
% [smAZ]=smoothazimuth(AZ,smfac)
%
% INPUTS 
% vector of azimuths in degrees 
% smoothing factor, should be an odd integer. e.g., 3,7,9 
%
% OUTPUT
% smAZ = running average vector mean aziuth. 
%
%
% must have function vectormean in the path 
% drbohnen@ncsu.edu 
% 09 July 2019 


azB=buffer(AZ,smfac,smfac-1);  % reshape the input into overlapping columns 
for i=1:smfac-1
azB(1:smfac-i,i)=AZ(i); %deal with first few columns that get zero padded  
end

% now call vector mean for each column and store the output. 
smAZ=nan(size(AZ)); 
for j=1:length(AZ)
smAZ(j)=vectormean(azB(:,j)); 
end

end


function [meanA, Rbar]=vectormean(azdeg)
% USE 
% [meanA, Rbar]=vectormean(azdeg)
% 
% INPUTS
% azdeg is a vector of aziuth (1 x M or M x 1) 
% where 360/0 deg is north 
% 
% OUTPUT 
% meanA is the vector mean direction 
% Rbar is the normalized resultant vector length 
% 
% drbohnen@ncsu.edu 
% 09 July 2019 

cosPhase = sum(cosd(azdeg));
sinPhase = sum(sind(azdeg));
R = sqrt(cosPhase^2 + sinPhase^2); %length of the combined vector squared
Rbar = R/length(azdeg);

%
if (sinPhase >=0 && cosPhase >=0) 
    meanA = atand(sinPhase/cosPhase);  
end
if (cosPhase < 0) 
    meanA = atand(sinPhase/cosPhase) + 180;
end
if (sinPhase <0 && cosPhase >=0) 
    meanA = atand(sinPhase/cosPhase)+360;
end

end


