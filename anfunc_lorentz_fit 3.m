function anfunc_lorentz_fit(dataloc,names,c3,stanbig,backper,lowercut,uppercut,lower,upper,nhood,rsquarelim,binfac)

%anfunc_lorentz_fit is a subfunction of hyper_analysis which processes
%hyperspectral images. This function is used primarily for memory
%management issue and passes no output back to the its parent function. It
%instead writes data to files named as the original TDMS filename followed
%by analysis.mat 1
%function dependencies: convertTDMS.m makergb.m partident.m
%fn_lorentz_fit.m plot_lorentz_fit.m 
%As of 08-21-13, framework for local background subtraction is provided but
%not verified. Some level of assurance must be reached to identify and
%exclude other particles, perhaps 2D gaussian fitting to a sum image, or a
%local version of the existing global subtraction.
    che=strcat(dataloc,'\',names{1,1}(c3+2),'.mat');
    che2=strcat(dataloc,'\',names{1,1}(c3+2),'global','.mat');
    che2_global=strcat(dataloc,'\',names{1,1}(c3+2),'local','.mat'); % A file including global bacsound
    fch=exist(che{1},'file');
    if fch==2
        load(che{1})
    else
        tdmsfile = strcat(names{1,1}(c3+2),'.tdms');
        [ConvertedData,~,~,~,~] = convertTDMS(1,tdmsfile,dataloc);   % TDMS to Matlab converter
    end
    grpnum = ConvertedData.Data.Root.Property(1,7).Value;   %wavelength grouping
    top = ConvertedData.Data.Root.Property(1,9).Value;      %top pixel used for binning
    bottom = ConvertedData.Data.Root.Property(1,10).Value;  %bottom pixel used for binning
    prow = ConvertedData.Data.MeasuredData(1,1).Property(1,10).Value-ConvertedData.Data.MeasuredData(1,1).Property(1,9).Value+1;
    pcol = ConvertedData.Data.MeasuredData(1,1).Property(1,8).Value;
    ncol = length(ConvertedData.Data.MeasuredData(4).Data); % number of columns in future array (CCD spanning wavelengths)
    npdata = zeros(prow*pcol,ncol);                        
    for c2 = 3:prow*pcol+2   % TDMS includes two empty structures for names/properties; data starts in third structure
        npdata(c2-2,:) = (ConvertedData.Data.MeasuredData(c2).Data');     % storing spectral data for one pixel as a row in matrix (turned horizontally)
    end
    wvlths=ConvertedData.Data.MeasuredData(end-2).Data';
    rgbspec=makergb(wvlths(1),wvlths(end),ncol,0.8);
    hypr=rgbspec(1,:);
    hypg=rgbspec(2,:);
    hypb=rgbspec(3,:);
    clear rgbspec
    sumimg = (zeros(prow,pcol));            % sum image matrix
    specim = (zeros(prow,pcol,ncol));       % 3D data cube
    bkavg = zeros(1,ncol);                  % background
    stanim = zeros(prow,pcol,ncol);         % 3D standard
    redim = zeros(prow,pcol,ncol);          % calculate R
    greenim = zeros(prow,pcol,ncol);        % calculate G
    blueim = zeros(prow,pcol,ncol);         % calculate B
    bkgim = zeros(prow,pcol,ncol);          % 3D global background
      if size(stanbig,3) == 1340
        grpnum = 2;
    else
        grpnum = 1;
    end
    for c4 = 1:pcol         % for all columns in the image
        for c5 = 1:prow     % for all rows in the image
            sumimg(c5,c4) = (sum(npdata(c5 + (c4-1)*prow,:)));  % summed data of pixel to produce image
            redim(c5,c4,:)= (hypr);
            greenim(c5,c4,:)= (hypg);
            blueim(c5,c4,:)= (hypb);
            specim(c5,c4,:) = ConvertedData.Data.MeasuredData(1,c5+(c4-1)*prow+2).Data';
        end
        stanim(:,c4,:)=stanbig(top:bottom,1,1:grpnum:end); %selects the local standard based on position on CCD
    end
    clear ConvertedData npdata
    
    smln=ceil(backper*prow*pcol);   % number of pixels to be included in background
%     imnorm=specim./stanim;          % Flatfield correct image first for summing
    imnorm=specim./1;          % Flatfield correct image first for summing
    sumnorm=sum(imnorm,3);          % creates sum image over all wavelengths
    sumnorm=(sumnorm-(min(min(sumnorm))))/(max(max(sumnorm-(min(min(sumnorm))))));  %normalize over this image
    sorted=unique(sumnorm);         % find the lowest index pixels
    nthsmlst=sorted(smln+1);        % find the cutoff value for background
    [a,b]=ind2sub(size(sumnorm), find(sumnorm<nthsmlst,smln)); %Record the backgound positions
    
    clear imnorm
    
    for d1 = 1:smln
        bkavg=bkavg+reshape(specim(a(d1),b(d1),:),1,ncol); % calculate global background sum
    end
    bkavg = bkavg/smln;             %global background average
    
    
    for c4 = 1:pcol                 % for all columns in the image
        for c5 = 1:prow             % for all rows in the image
            bkgim(c5,c4,:) = bkavg; % makes a background cube
        end
    end
    
%    specflat = (specim ./ stanim);              % used if doing local bg     
    specfin = (specim - bkgim) ./ stanim;       % spectral correction
    
    clear bkgim bkgim
    
    specrgb(:,:,1) = (sum(specfin(:,:,lowercut+1:end-uppercut).*redim(:,:,lowercut+1:end-uppercut),3));     % build RGB
    specrgb(:,:,2) = (sum(specfin(:,:,lowercut+1:end-uppercut).*greenim(:,:,lowercut+1:end-uppercut),3));
    specrgb(:,:,3) = (sum(specfin(:,:,lowercut+1:end-uppercut).*blueim(:,:,lowercut+1:end-uppercut),3));
    
    clear redim greenim blueim 
    
    specrgbnorm=(specrgb+abs(min(min(min(specrgb)))))/max(max(max(specrgb))); % normalize the RGB image
    
    %% finding particles

    %nhood = 5;                    % nhood by nhood pixels are grouped for each particle spectrum

    [ptu mm mark2] = partident(specrgbnorm,lower,upper,nhood);    % find particles based on a first comer algorithm
%     mm = 1;
%     mark2 = [148,204];
    %mm is the number of particles found and mark2 is a matrix with their
    %locations
    
    spec=zeros(mm,ncol);        % calculates a spectrum for each identified particle
    %specpad=padarray(specim,[nhood nhood],0,'both');  %pad so as to not run over the edge
    
   % dist = zeros(prow,pcol,mm);
    %extent = zeros(prow,pcol,mm);
    %yprime = zeros(ncol-lowercut-uppercut,mm);
    params = zeros(mm,10);
    
    for k=1:mm 
        for ni=0:nhood-1
            for nj=0:nhood-1
                spec(k,:)=spec(k,:)+reshape(specim(mark2(k,2)-(nhood-1)/2 + ni,mark2(k,1)-(nhood-1)/2 + nj,:),1,ncol);
            end
        end
        
        flats = reshape(stanbig(top+(mark2(k,2)),1,1:grpnum:end),size(bkavg));
%         flats = reshape(stanbig(top-(mark2(k,2)),1,1:grpnum:end),size(bkavg));
        speccor(k,:)=(spec(k,:)-nhood^2*bkavg)./ flats;
        [c,g,xbin]=fn_lorentz_fit_bin(wvlths(lowercut+1:end-uppercut)',speccor(k,lowercut+1:end-uppercut)',1,binfac);
        yprime(:,k)=feval(c,xbin);  %fit to each particle
        
        % params are values for fit as defined in fn_lorentz_fit and have
        % descriptive names
        params(k,1)=c.a1;   
        params(k,2)=c.b1;
        params(k,3)=c.c1;
        params(k,4)=g.sse;
        params(k,5)=g.rsquare;
        params(k,6)=g.dfe;
        params(k,7)=g.adjrsquare;
        params(k,8)=g.rmse;
        params(k,9)=1240/c.b1;
        params(k,10)=1240/(c.b1-0.5*c.c1)-1240/(c.b1+0.5*c.c1);
    end
        i=0;
        l=0;
        
        mark2discard=[];
    for j=1:mm
        if params(j,5)> rsquarelim && params(j,2) > wvlths(lowercut+1) && params(j,2) < wvlths(end-uppercut) 
            i=i+1;
            paramssort(i,:)=params(j,:);
            mmkeep(i)=j;
            yprimesort(i,:)=yprime(:,j)';
            speccorsort(i,:)=speccor(j,:);
            mark2sort(i,:)=mark2(j,:);
        else
            l=l+1;
            mmthrow(l)=j;
            mark2discard(l,:)=mark2(j,:);
        end
    end
    newmm=i;
    badmm=l;
    
        
    
    
        %% local bg section under construction and commented out
        
        %{
        slices(:,:,k) = specfin(:,:,peakloc(k)+lowercut);

        for sc=1:pcol
            for sr = 1:prow
                dist(sr,sc,k) = sqrt((sr-mark2(k,2))^2+(sc-mark2(k,1))^2);
                extent(sr,sc,k) = slices(sr,sc,k) / specfin(mark2(k,2),mark2(k,1),peakloc(k)+lowercut);
            end
        end
    end
    
    
    dist(:,:,mm+1) = max(dist,[],3);
    extent(:,:,mm+1) = max(extent,[],3);
    
    figure; imagesc(dist(:,:,mm+1))
    title('dist')
    figure; imagesc(extent(:,:,mm+1))
    title('extent')
    for kl=1:mm
        weigh(:,:,kl) = (extent(:,:,mm+1)<0.1)./((dist(:,:,kl)/4).^10);
    end
   
    %figure; imagesc(weigh(:,:,1))
    
    %{
    figure;  imagesc(dist{c3}(:,:,mm+1)); title('proximity')
    colorbar
    figure;  imshow(extent{c3}(:,:,mm+1),[0 0.1],'InitialMagnification','fit'); title('extent maxima')
    colorbar; colormap jet
    %}
    
    spfl = zeros(mm,ncol);
    bkloc = zeros(mm,ncol);
    bgn=10;
    locback = zeros(bgn,2,mm);
    peak2 = zeros(mm);
    miny2 = zeros(mm);
    

                
    for l=1:mm
        %[locback{c3}(:,1,l),locback{c3}(:,2,l)]=ind2sub([prow,pcol], find(weigh{c3}(:,:,l)<0,bgn,'last'));
        %[locback{c3}(:,1,l),locback{c3}(:,2,l)]=max(weigh{c3}(:,:,l)<0,bgn,'last'));
        
        ordrd = flipud(unique(weigh(:,:,l)));
        ordrd(isnan(ordrd)) = [];
        l
        [locback(:,1,l),locback(:,2,l)] = ind2sub(size(weigh(:,:,1)), find(weigh(:,:,l)>(ordrd(bgn+1)),bgn));
        
        
        
        %figure; imshow(weigh{c3}(:,:,l),'InitialMagnification','fit')
        %caxis auto; colormap(hot);
        %hold on

        for d2 = 1:bgn
            bkloc(l,:) = bkloc(l,:)+reshape(specflat(locback(d2,1,l),locback(d2,2,l),:),1,ncol);
            %plot(locback{c3}(d2,2,l),locback{c3}(d2,1,l),'gx','MarkerSize',3,'LineWidth',2)
        end
        
       % hold off
        
        bkloc(l,:) = bkloc(l,:)/bgn;
        %bkfl{c3}(l,:) = bkloc{c3}(l,:)./double(standard);
        
        %this section recalculates the spectra with local bg subtraction
        
        for ni=1:nhood
            for nj=1:nhood
                if flatpad(mark2(l,2)-2+ni+nhood,mark2(l,1)-2+nj+nhood,1)==0
                    spfl(l,:)=spfl(l,:)+reshape(flatpad(mark2(l,2)-2+ni+nhood,mark2(l,1)-2+nj+nhood,:),1,ncol);
                else
                    spfl(l,:)=spfl(l,:)+reshape(flatpad(mark2(l,2)-2+ni+nhood,mark2(l,1)-2+nj+nhood,:),1,ncol)-bkloc(l,:);
                end
            end
        end
        %spflst{c3}=spfl./standard;
        
        
        [yprime2(l,:) params2(l,:) resnorm2(l) residual2(l,:)]=lorentzfit(wvlths(lowercut+1:end),spfl(l,lowercut+1:end));
        [peak2(l),~] = max(yprime2(l,:));
        miny2(l) = mean(mean(spfl(l,1:10))+mean(spfl(l,end-10:end)));  %min(yprime2(l,:));
        %half2{c3}(l,:) = abs(yprime2(l,:)-miny2(l)-(peak2(l)-miny2(l))/2);
        %fwhm2{c3}(l) = wvlths{c3}(1,lowercut+find(half2{c3}(l,:)==min(half2{c3}(l,peakloc2{c3}(l):end))))-wvlths{c3}(1,lowercut+find(half2{c3}(l,:)==min(half2{c3}(l,1:peakloc2{c3}(l)))));
        %ratpeak2{c3}(l) = abs(peak2(l)/miny2(l));
    end
    %}
    %{
    figure
    imshow(sumnorm,[min(min(sumnorm)) max(max(sumnorm))],'InitialMagnification','fit')
    colormap(hot)
    hold on
    for ls = 1:mm
        for d3 = 1:bgn
            plot(locback{c3}(d3,2,ls),locback{c3}(d3,1,ls),'gx','MarkerSize',3,'LineWidth',2)
        end
         
    end
    title([names{1,1}(c3+2)])
    hold off
   
    
    figure
  
    imshow(specrgbnorm,'InitialMagnification','fit')
    
    
    figure
    imshow(specrgbnorm,'InitialMagnification','fit')
    hold on
    ColOrd = get(gca,'ColorOrder');
    [m,~] = size(ColOrd);
    for ki = 1:mm
        ColRow = rem(ki,m);
        if ColRow == 0
            ColRow = m;
        end
        Col = ColOrd(ColRow,:);
        plot(mark2(ki,1),mark2(ki,2),'o','MarkerSize',5,'LineWidth',1,'Color',Col)
        hold on
    end
    hold off
    title(num2str(mm))
    %}
    
    %figure
    %plot(results(:,1),results(:,2))
  
    
%    save(che2{1},'specfin','speccor','specrgbnorm','mark2','mm','wvlths','yprime','params','paramssort','newmm','badmm','yprimesort','speccorsort','mark2sort','mark2discard','xbin');
   save(che2{1},'specfin','speccor','wvlths');
   % Save the additional file including global background
   specfin = specim ./ stanim;
   save(che2_global{1}, 'specfin', 'speccor', 'wvlths'); % Keeping original specim before global correction
 %%
   
  

    

end
