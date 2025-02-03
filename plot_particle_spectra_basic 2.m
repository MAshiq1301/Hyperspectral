%Plots the spectra of particles at specified pixels
close all 
clear all
%% Parameters
%Specify path to folder with analysis data
WINSD='C:\CHRISTY\ASHIQ\DFData\250122 spectrum\AuNR1\';
sample = '012225_AuNR_22by77_5analysislocal';
file = 'AuNS803';
cd(strcat(WINSD));
%Specify file name of data you want to plot
load(sample, 'specfin', 'wvlths')
place2save = 'C:\CHRISTY\ASHIQ\DFData\250122 spectrum\AuNR1'
%automatedIUC\data\hyperspectral\240212\2';
mkdir(place2save)

Rawsave = 'C:\CHRISTY\ASHIQ\DFData\250122 spectrum\AuNR1';
Nanosize = '45by85';
%Specify pixel coordinates
colorArea = ''; % Red Orange or Green
Areaval = 'AuNR';

% [part_coord,specfin,wvlths] = particle_search(20,WINSD,file,5);
%manual input
part_coord =[




];
temp_coord =part_coord;
part_coord = circshift(part_coord,[0,1]);
matrix_index_ref = part_coord;
matrix_index_ref = circshift(matrix_index_ref,[0,1]);
% specfin = specfin(:,:,141:1000);
specfin_final = sum(specfin(:,:,141:670),3);

%Testing area for shift code
DriftFind='n';
if DriftFind =='y'
    shifter =3;
    for i= 1:length(part_coord)
        y=(part_coord(i,1));
        x=(part_coord(i,2));
        CheckRange = specfin_final(y-shifter:y+shifter,x-shifter:x+shifter);
        maxValue = max(CheckRange(:));
        [y2, x2] = find(CheckRange == maxValue);
        ny = x+(x2-4);
        nx = y+(y2-4);
        part_coord(i,1)=nx;
        part_coord(i,2)=ny;
    end
    newdriftbase=[part_coord(:,2),part_coord(:,1)];
end

high = max(max(specfin_final));
figure1 = figure;
imshow(specfin_final,[0 high/1])
hold all
for n = 1:size(part_coord,1)
    text(part_coord(n,2),part_coord(n,1),num2str(n),'Color','Green')
end
hgsave('image','-v6')


%Specify cutoff bounds in units of pixels, not nm. 200:1050 gives about 
% 500-850 nm
%NORMAL IS 100 to 670
lower_bound =100;
upper_bound = 670;
rawwvlths = wvlths(100:670);
wvlths = wvlths(lower_bound:upper_bound);
       
%Describes the number of pixels you want to integrate/particle. Typically 3 or 5.
int_size = 3;                                              
int_var_low = -1*(int_size-1)./2;
int_var_high = (int_size-1)./2;

% n = size(matrix_index_ref,1);
% x_part = matrix_index_ref(n+size(matrix_index_ref,1));
% y_part = matrix_index_ref(n);
% part_spec = zeros(1,1,670);
% for m = int_var_low:int_var_high
%     for l = int_var_low:int_var_high
%         a = x_part + m;
%         b = y_part + l;
%         part_spec = part_spec + specfin(a,b,:);
%     end
% end
% part_spec = squeeze(part_spec);
% background = part_spec(lower_bound:upper_bound);

% % % Remove specfic matrix
% temp = [];
% for i = [1,2,6:15];
%     temp = [temp;matrix_index_ref(i,:)];
% end
% matrix_index_ref = temp;


%% Integrates over the specified square and then plots each specified pixel.
diff = [];
cd(place2save)
all_spec = [];
AllSigyNoi =[];
rawdata = [];
fits =[]; 
% Extract the last part of the filename
if endsWith(sample, 'global')
    use_local_correction = false; % Deactivate background subtraction
elseif endsWith(sample, 'local')
    use_local_correction = true;  % Activate background subtraction
else
    error('Sample name must end with either "global" or "local".');
end
for n = 1:size(matrix_index_ref,1)
        % full data background correcting
    x_part = matrix_index_ref(n+size(matrix_index_ref,1));
    y_part = matrix_index_ref(n)+7;
    part_spec = zeros(1,1,670);
    for m = int_var_low:int_var_high
        for l = int_var_low:int_var_high
            a = x_part + m;
            b = y_part + l;
            part_spec = part_spec + specfin(a,b,:);
        end
    end
    part_spec = squeeze(part_spec);
    background = part_spec(1:670);
    
    x_part = matrix_index_ref(n+size(matrix_index_ref,1));
    y_part = matrix_index_ref(n);
    part_spec = zeros(1,1,670);
    for m = int_var_low:int_var_high
        for l = int_var_low:int_var_high
            a = x_part + m;
            b = y_part + l;
            part_spec = part_spec + specfin(a,b,:);
         end
    end
    part_spec = squeeze(part_spec);
    part_specr = part_spec(100:670);
    rawdata = [rawdata,part_specr];
    % fit area
    x_part = matrix_index_ref(n+size(matrix_index_ref,1));
    y_part = matrix_index_ref(n)+7;
    part_spec = zeros(1,1,670);
    for m = int_var_low:int_var_high
        for l = int_var_low:int_var_high
            a = x_part + m;
            b = y_part + l;
            part_spec = part_spec + specfin(a,b,:);
        end
    end
    part_spec = squeeze(part_spec);
    background = part_spec(lower_bound:upper_bound);
    
    x_part = matrix_index_ref(n+size(matrix_index_ref,1));
    y_part = matrix_index_ref(n);
    part_spec = zeros(1,1,670);
    for m = int_var_low:int_var_high
        for l = int_var_low:int_var_high
            a = x_part + m;
            b = y_part + l;
            part_spec = part_spec + specfin(a,b,:);
         end
    end
    part_spec = squeeze(part_spec);
    part_spec = part_spec(lower_bound:upper_bound);
     % Conditionally subtract background if sample name ends with 'local'
    if use_local_correction
        part_spec = part_spec - background;
    end
    all_spec = [all_spec,part_spec];
    
    [param_1,param_2]=fn_lorentz_fit(wvlths',part_spec,1,1);
    a1 = param_1.a1;
    b1 = param_1.b1;
    c1 = param_1.c1;
    resonance = b1;
    FWHM = c1;
    r_list = param_2.rsquare;
    lorentz_fit =(2*a1/pi).*(c1./(4*(wvlths'-b1).^2+c1.^2));
    diff = part_spec-lorentz_fit;
    Noi = std(diff);
    [Notneeded, IndiMax]=min(abs(wvlths-resonance));
    Sigy = part_specr(IndiMax);
    SnN= Sigy/Noi;
    AllSigyNoi = [AllSigyNoi,[Noi;Sigy;SnN]];
    figure1 = figure;
    hold all
    plot(rawwvlths,part_spec,'b','linewidth',3)
    plot(wvlths,lorentz_fit,'k--','linewidth',3)
    xlabel('Wavelength (nm)','fontsize',32)
    ylabel('Scattering','fontsize',32)
    set(gca,'FontSize',22,'box','on')
 Energy(n)= round(resonance);
     Width(n)= round(FWHM);
     SNR(n)= round(SnN);
    text(0.55,0.9,['\lambda_m_a_x = ',num2str(round(resonance)), ' nm'],'fontsize',20,'Units','normalized')
    text(0.55,0.78,['\Gamma = ',num2str(round(FWHM)), ' nm'],'fontsize',20,'Units','normalized')
    text(0.55,0.66,['S/N = ',num2str(round(SnN))],'fontsize',20,'Units','normalized')
    xlim([450 950])
    ylim([0 max(part_spec)+0.02])
    title(num2str([a1 b1 c1]))
%     ylim([0 0.7])
    filename1 = [Nanosize,colorArea,Areaval,num2str(n)];
    saveas(figure1,[filename1,'.tif'])
    hgsave(filename1,'-v6')
    close;
        figure2 = figure;
    hold all
    plot(rawwvlths,part_spec,'b','linewidth',3)
    %plot(wvlths,lorentz_fit,'k--','linewidth',3)
    xlabel('Wavelength (nm)','fontsize',32)
    ylabel('Scattering','fontsize',32)
    set(gca,'FontSize',22,'box','on')
%     text(0.55,0.9,['\lambda_m_a_x = ',num2str(round(resonance)), ' nm'],'fontsize',20,'Units','normalized')
    %text(0.05,0.78,['\Gamma = ',num2str(round(FWHM)), ' nm'],'fontsize',20,'Units','normalized')
    xlim([450 950])
    ylim([-.02 max(part_spec)+0.02])
%     ylim([0 0.7])
    filename2 = [Nanosize,colorArea,Areaval,num2str(n),'r'];
    hgsave(filename2,'-v6')
    temp1 = [wvlths.',part_spec];
    filename3 =[sample,'_','spectra','_',num2str(n),'r'];
    if Rawsave == 'y'
        saveas(figure2,[filename2,'.tif'])
    end
    close;
end
save('spectra','wvlths','all_spec','Energy','Width','SNR','matrix_index_ref')



