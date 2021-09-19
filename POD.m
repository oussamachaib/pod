clc;
close all;
clear;

%% Image location

% Folder
directory='E:\Projet de synthèse\Manips\28.01.2020\20200128_Moog40Hz-5Off-2,5Amp_reactif_C12_2\';

% Image title
root_name='20200128_Moog40Hz-5Off-2,5Amp_reactif_C12';

% Saved .mat file

mkdir '20200128_Moog40Hz-5Off-2,5Amp_reactif_C12_2_derniers';
name_tosave='20200128_Moog40Hz-5Off-2,5Amp_reactif_C12_2';

%Display  1
    k=1;
    image_name=num2str(k,'%0.6d');
    image=imread([directory,root_name,image_name,'.tif']);
    figure(2);
    imagesc(image);


%% Initial parameters

% Nombre d'images
N=2000;

% Fréquence d'acquisition
freq=4000;

% Cropping
    % (x1,y1): top left
    % (x2,y2): bottom right
    
% %No binning
% x1=274;
% y1=5;
% x2=844;
% y2=823;

%Binning 1/4 ----------------------------------------------- Find cropping
%limits on the binned image
b=1;

x1=2;
y1=50;
x2=1023;
y2=360;

    % Cropping limits
cropping=[x1 y1 x2-x1+1 y2-y1+1];


%% Initialization of vectors/matrices

% Image matrix
    % Pixel values ordered in columns
    % Rows = images from 1 to N
image_reshape=zeros((cropping(3)+1)*(cropping(4)+1),N);


%% POD Computation: Part 1

%       Loading chemiluminescence data / filling the reshaped image matrix
disp('Loading chemiluminescence data...');
tic;
for k=748:N
    image_name=num2str(k,'%0.6d');
    image=imread([directory,root_name,image_name,'.tif']);
%Binning
%    image=imresize(image,b);
    image=imcrop(image,cropping);    
    image_reshape(:,k)=reshape(image,[(cropping(3)+1)*(cropping(4)+1),1]);
    k
    imagesc(image);
end

disp('Chemiluminesence data loaded successfully.');

%       Computing the correlation matrix
disp('Computing correlation matrix R...');

R=image_reshape'*image_reshape;                                                 % x.x' memory-hungry. We use a time autocorrelation instead (x'.x) (PPT Richecoeur)

%       Computing eigenvalues and eigenvectors
disp('Computing eigenvalues and eigenvectors...');

[eigvec,eigval]=eig(R);

%       Computing eigenvalue size
[L,c]=size(eigval);
[L1,c1]=size(eigvec);

%       Injecting eigenvalues (diagonal of the matrix) in a new vector
eigval_vector=zeros(L,c);

disp('Sorting eigenvalues...');

%       Flipping the matrix horizontally

for k=1:L
    eigval_vector(k)=eigval(k,k);
end

%       Sorting eigenvalues in descending order
[eigval_vector,order]=sort(eigval_vector,'descend');

%       Sorting eigenvectors
disp('Sorting eigenvectors...');


eigvec_ordered=zeros(L1,N);

for mode=1:N
    p=order(mode);
    eigvec_ordered(:,mode)=eigvec(:,p);
end

toc;
%% POD Computation: Part 2 (Computing modes and coefficients)
tic;
%       Computing modes and coefficients
disp('Computing modes and coefficients...');

Phi=image_reshape*eigvec_ordered;
a=Phi'*image_reshape;

%       Reordering modes
disp('Sorting modes...');

for mode=1:6
    Phi(:,mode)=Phi(:,mode)/eigval_vector(mode);
end
toc;
%% Display
tic;
%       Displaying modes  
disp('Displaying modes...');
Disp_Phi=Phi(:,1:6);

%       Reshaping modes into images
Disp_Phi=reshape(Disp_Phi,(cropping(4)+1),(cropping(3)+1),6);

%       Displaing modes

for mode=1:6
    image=Disp_Phi(:,:,mode);
    image=(image-min(min(image)));
    image=image/max(max(image));
    figure(mode);
    imagesc(image);
    title(['Mode ',num2str(mode)]);
    xlabel('x');
    ylabel('y');
    mode
saveas(gcf,[pwd '/',name_tosave,'/',name_tosave,'_mode',num2str(mode),'.png']);
saveas(gcf,[pwd '/',name_tosave,'/',name_tosave,'_mode',num2str(mode),'.fig']);
end

%       Computing energy content of each mode
%Energy = Sum of eigenvalues

% Etot=0;
% for i=1:L
%     Etot=Etot+eigval_vector(i);
%     
% end
% Sum=0;
% for i=1:L
%     energy(i-1)=eigval_vector(i)*100/Etot;
%     Sum=Sum+eigval_vector(i);
%     cont(i)=Sum*100/Etot;
% end
 
Etot=0;
for i=1:L
    Etot=Etot+eigval_vector(i);
    
end
Sum=0;
for i=1:L
    energy(i)=eigval_vector(i)*100/Etot;
    Sum=Sum+eigval_vector(i);
    cont(i)=Sum*100/Etot;
end

figure(20);
bar(2:21,energy(2:21)*100/sum(energy(2:21)))
grid;
xlabel('Modes');
ylabel('Energy content (%)');
title('Energy content (%)');
saveas(gcf,[pwd '/',name_tosave,'/',name_tosave,'_energy.png']);
saveas(gcf,[pwd '/',name_tosave,'/',name_tosave,'_energy.fig']);
%       Log scale
% figure;
% h=loglog(energy,'b');
% xlabel('Modes');
% ylabel('Energy content (%)');

%       Energy content (%)
% h=sum(energy(2:end));
% en_percent=energy(2:end)*100/h;
% figure;
% bar(en_percent,'b')
% xlabel('Modes');
% ylabel('Energy content (%)');
% grid;
% %ylim([0 30]);

toc;
%% Spectra
tic;
samples=N;
dt=1/freq;
warning off;
window=4; 
for mode=1:6
    figure(10+mode);
    [PxxU,f] = pwelch(detrend(a(mode,:)),hanning(round(samples/window)),0.5,round(samples/window),1/dt);
    loglog(f,PxxU/max(PxxU),'k-');
    axis([0 1000 -Inf Inf])
    xlabel('Frequency [Hz]');
    ylabel('Power spectral density [-]');
    set(gca,'FontName','TimesNewRoman','FontSize',20);
saveas(gcf,[pwd '/',name_tosave,'/',name_tosave,'_PSDmode',num2str(mode),'.png']);
saveas(gcf,[pwd '/',name_tosave,'/',name_tosave,'_PSDmode',num2str(mode),'.fig']);
end
save(name_tosave);
 
% figure(7);plot(a(2,:),a(1,:),'ko');hold on
% figure(7);plot(a(2,1:500),a(1,1:500),'ko','MarkerFaceColor','g');
% figure(7);plot(a(2,end-500:end),a(1,end-500:end),'ko','MarkerFaceColor','r');
% xlabel('Mode 1');ylabel('Mode 2');
% set(gca,'FontName','TimesNewRoman','FontSize',20);

modex=3;
modey=2;

figure(8);plot(a(modex,:),a(modey,:),'ko');hold on
figure(8);plot(a(modex,1:500),a(modey,1:500),'ko','MarkerFaceColor','g');
figure(8);plot(a(modex,end-500:end),a(modey,end-500:end),'ko','MarkerFaceColor','r');
xlabel(['Mode ' num2str(modex)]);ylabel(['Mode ' num2str(modey)]);
set(gca,'FontName','TimesNewRoman','FontSize',20);
saveas(gcf,[pwd '/',name_tosave,'/',name_tosave,'.png']);
saveas(gcf,[pwd '/',name_tosave,'/',name_tosave,'.fig']);
close all;

toc;



