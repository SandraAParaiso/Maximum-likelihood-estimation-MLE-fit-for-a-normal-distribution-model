% Example analyzing psychophysical data, fitting a model, and visualizing the results for one subject. 
% It involves statistical estimation and graphical representation of the relationship between stimulus duration and perceptual responses.

close all
clear all

addpath('./DATA/');
files=dir('./DATA/*.mat');

for kk=1:length(files);

addpath('./DATA/');
files=dir('./DATA/*.mat');
    
    
load(files(kk).name);

if length([data.trials])==270;

%%%%%%%%%%%

stairs=[ 1 2 3 4 5 6 7 8 9];
A=(find([data.typestairs]==stairs(1)));
B=(find([data.typestairs]==stairs(2)));
C=(find([data.typestairs]==stairs(3)));
D=(find([data.typestairs]==stairs(4)));
E=(find([data.typestairs]==stairs(5)));
F=(find([data.typestairs]==stairs(6)));
G=(find([data.typestairs]==stairs(7)));
H=(find([data.typestairs]==stairs(8)));
I=(find([data.typestairs]==stairs(9)));

Dur=[A B C D E F G H I];
duration=unique(log10([data(Dur).duration]));
duration_presented=log10([data(Dur).duration]);
responses=[data(Dur).correcterror];

for h=1:length(duration);
   totaldata=find(duration_presented<duration(h)+0.00001&duration_presented>duration(h)-0.00001);
   correct=responses(totaldata);
   numbertotal(h)=length(totaldata);
   nYes(h)=sum(correct);
   nNo(h)=numbertotal(h)-nYes(h);
   probcorrect(h)=sum(correct)./length(totaldata);
   
end

addpath('./DATA');

global simulation

simulation=[duration' probcorrect' nYes' nNo' numbertotal'];
H=find(simulation(:,5)<30); %<30
simulation(H,:)=[];

start = [(mean((simulation(:,1)))) 0.5];
options = optimset('TolX',1e-8,'MaxFunEvals',4000, 'MaxIter',2000);
[parametros_estimados,fval,exitflag,output] = fminsearch('adjustML_norm',start,options);

mu1=parametros_estimados(1);
sigma1=parametros_estimados(2);
pg1=0.5;;
umbralT = 10^mu1;

[II]=find([data.typestairs]==stairs(1));
[JJ]=find([data.typestairs]==stairs(2));
[KK]=find([data.typestairs]==stairs(3));
[LL]=find([data.typestairs]==stairs(4));
[MM]=find([data.typestairs]==stairs(5));
[NN]=find([data.typestairs]==stairs(6));
[OO]=find([data.typestairs]==stairs(7));
[PP]=find([data.typestairs]==stairs(8));
[QQ]=find([data.typestairs]==stairs(9));

aciertos=[data.correcterror];

datosA=[[1:length(aciertos(II))];aciertos(II)]';
datosB=[[1:length(aciertos(JJ))];aciertos(JJ)]';
datosC=[[1:length(aciertos(KK))];aciertos(KK)]';
datosD=[[1:length(aciertos(LL))];aciertos(LL)]';
datosE=[[1:length(aciertos(MM))];aciertos(MM)]';
datosF=[[1:length(aciertos(NN))];aciertos(NN)]';
datosG=[[1:length(aciertos(OO))];aciertos(OO)]';
datosH=[[1:length(aciertos(PP))];aciertos(PP)]';
datosI=[[1:length(aciertos(QQ))];aciertos(QQ)]';

%%%Results of the simulation
pl1=0;
x=[0:0.001:3];
pYes=pg1+((1-pg1-pl1).*normcdf(x,mu1,sigma1));   

threshold_prob=0.8200;
threshold=norminv((threshold_prob-pg1)/(1-pg1-pl1),mu1,sigma1)


figure('NumberTitle', 'off', 'Name', files(kk).name);

plot(duration,probcorrect,'bo','MarkerSize',6,'MarkerFaceColor',[0 0 1],'Color',[0 0 1],'LineWidth',3);
hold on
plot(x,pYes,'k-','MarkerSize',6,'MarkerFaceColor',[0/255 206/255 209/255],'Color',[0/255 206/255 209/255],'LineWidth',3);
hold on
plot([mu1 mu1],[0 threshold_prob],'g--','Color',[205/255 205/255 205/255],'LineWidth',2);
hold on
plot([0 mu1],[threshold_prob threshold_prob],'g--','Color',[205/255 205/255 205/255],'LineWidth',2);
hold on
% plot(0*ones(11,1),[0:0.1:1],'k--'); %
xlabel('log_1_0(Duration (mseg))');
ylabel('Proportion correct');
xlim([0 2.6]);

%Plot
set(gca,'Fontsize', 12,'FontWeight','bold','Fontname', 'Times New Roman');
ax = gca;
ax.LineWidth = 2;

for k=1:length(numbertotal);
text(duration(k)-0.1,probcorrect(k)-0.03,sprintf(' %d',numbertotal(k)),'col',[0 0 0],'Fontsize', 10,'Fontname', 'Times New Roman');%Represents the number of samples per point
end

izq=0.1+mu1;

date=[files(kk).name(20:length(files(2).name))];

subject = 'Subject: %s';
subjectname=[files(kk).name(16:18)];
text(izq,0.4,sprintf(subject,subjectname),'col',[0 0 0],'Fontsize', 10,'Fontname', 'Times New Roman')

age = 'Age: %s y';
age_number = expt.age;
text(izq,0.35,sprintf(age,age_number),'col',[0 0 0],'Fontsize', 10,'Fontname', 'Times New Roman')

if expt.window==1
    gauss='Spatial deviation Gaussian: %u deg';
    gauss_n=expt.spatialdeviationgaussian;
    text(izq,0.3,sprintf(gauss,gauss_n),'col',[0 0 0],'Fontsize', 10,'Fontname', 'Times New Roman')
end
if expt.window==2
    butter='Butterworth diameter: %u deg';
    butter_n=expt.butterworth_diameter;
    text(izq,0.3,sprintf(butter,butter_n),'col',[0 0 0],'Fontsize', 10,'Fontname', 'Times New Roman')
end

freq = 'Frequency: %u c/deg';
freq_number = expt.frequency;
text(izq,0.25,sprintf(freq,freq_number),'col',[0 0 0],'Fontsize', 10,'Fontname', 'Times New Roman')

speed = 'Speed: %u deg/sec';
Speed_number = expt.speed;
text(izq,0.2,sprintf(speed,Speed_number),'col',[0 0 0],'Fontsize', 10,'Fontname', 'Times New Roman')

Contrast = 'Michelson contrast of %d%%';
Contrast_number = expt.contrast.*100;
text(izq,0.15,sprintf(Contrast,Contrast_number),'col',[0 0 0],'Fontsize', 10,'Fontname', 'Times New Roman')

Thres = 'Threshold (%d%%): %1.4f mseg';
cresponses=(threshold_prob*100);
th = 10^threshold;
text(izq,0.1,sprintf(Thres,cresponses,th),'col',[0 0 0],'Fontsize', 10,'Fontname', 'Times New Roman')

sigma_name = 'Sigma: %1.4f';
sigm = sigma1;
text(izq,0.05,sprintf(sigma_name,sigm),'col',[0 0 0],'Fontsize', 10,'Fontname', 'Times New Roman')

end 
    
end

