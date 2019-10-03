%computes CSD of toy data (PP stimulation recorded by Ain Chung 2015)
%used CSDplotter library
%Dino Dvorak (dino@indus3.net)
%assumes FS = 2000;
clear; close all; clc;

addpath('CSDplotter');
addpath('CSDplotter/methods/');

%select electrode spacing (in um)
electrodeSpacing = 50;
%electrodeSpacing = 100;

Nchannels = 16;

methodSpline = 1; %0-standard,1-spline(smooth)

%example data
load('dataSTIM.mat')
%%
%CSD properties
if electrodeSpacing == 50
    h = 50e-6; %electrode distance m
    gauss_sigma = 0.05e-3; %mm
    el_pos = 0.05:0.05:(Nchannels-1)*0.05+0.05; %distances
elseif electrodeSpacing == 100
    h = 100e-6; %electrode distance m
    gauss_sigma = 0.1e-3; %mm
    el_pos = 0.1:0.1:(Nchannels-1)*0.1+0.1; %distances
end
ex_cond = 0.3; %external conductivity S/m
top_cond = 0.3; %S/m
diam = 0.5e-3; %mm
filter_range = 5*gauss_sigma; % numeric filter must be finite in extent
dt = 0.5; %sampling time s (1000/eegFS)
scale_plot = 1; %focus
max_plot = 0;
el_pos = el_pos*1e-3; %mm


if methodSpline == 0
    %compute CSD - standard method
    %adding vaknin electrodes not to loose first and last channel
    %duplicate first and last channel
    eegAver = cat(1,eegAver(1,:),eegAver,eegAver(end,:));

    Nch = size(eegAver,1);
    CSD = -ex_cond*D1(Nch,h)*eegAver;

    %filter iCSD
    b0 = 0.54; %center
    b1 = 0.23; %neighbor
    [n1,n2]=size(CSD);            
    CSD_add(1,:) = zeros(1,n2);   %add top and buttom row with zeros
    CSD_add(n1+2,:)=zeros(1,n2);
    CSD_add(2:n1+1,:)=CSD;        %CSD_add has n1+2 rows
    CSD = S_general(n1+2,b0,b1)*CSD_add; % CSD has n1 rows

    %plot CSD
    plot_CSD(CSD,el_pos,dt,scale_plot,max_plot)
else

    % compute spline iCSD:
    Fcs = F_cubic_spline(el_pos,diam,ex_cond,top_cond);
    [zs,CSD_cs] = make_cubic_splines(el_pos,eegAver,Fcs);
    if gauss_sigma~=0 %filter iCSD
      [~,CSD_cs]=gaussian_filtering(zs,CSD_cs,gauss_sigma,filter_range);
    end

    %plot CSD
    plot_CSD(CSD_cs,zs,dt,1,0) %length(el_pos) must equal rows of CSD!
end

axis ij tight
colormap(jet)

