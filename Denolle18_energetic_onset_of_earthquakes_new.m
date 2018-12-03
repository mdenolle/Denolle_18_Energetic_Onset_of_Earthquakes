clear all
close all
clc
DIR='~/Dropbox/SOURCE/SCARDEC/';
% DIR='YOU_PATH_TO_SCARDEC/';
F = importdata('./scardec_filenames.dat');
dt1=0.0703;
tt1=(-10:dt1:1000);
dt2=0.01;tt2=(-0.05:dt2:2);

stf=zeros(length(F),length(tt1));er1=stf;erstf=stf;

stf_sd=zeros(length(F),length(tt2));er_sd=stf_sd;erstf_sd=stf_sd;
stf_sm=stf_sd;er_sm=stf_sm;erstf_sm=stf_sd;

mo=zeros(length(F),1);AMtype=mo;mw=mo;
E=mo;E0=mo;TD5=mo;Es=mo;fbest1=mo;fbest2=mo;
mo1=mo;TD=mo;TD2=mo;TD3=mo;TM=mo;Esum=mo;TD4=mo;Ecorr=mo;
ER=mo;Tpeak_stf=mo;Tpeak_er=mo;Tpeak_erstf=mo;
A_peak_er=mo;A_peak_stf=mo;A_peak_erstf=mo;
depth=mo;r_p=mo;a_p=mo;mu_p=mo;b_p=mo;ERfac_p=mo;
A_peak_erstf_sd=mo;A_peak_er_sd=mo;A_peak_stf_sd=mo;
sKd=mo;sKm=mo;Kd=mo;Km=mo;tgrowth=mo;
Es_final_test=mo;
fit_skG=mo;

% FF1=logspace(log(2),log10(50),500);
% G1 = [ones(length(FF1),1) -log10(FF1)'];
FF0=logspace(-3,1,500);
FF=logspace(-2,1.5,100);
% G = [ones(length(FF),1) -log10(FF)'];
% G10 = [ones(length(FF0),1) -log10(FF0)'];
ff1=logspace(-3,0,50);
hffr=linspace(1.5,4,50);
xspec2=zeros(length(F),length(FF));
xspec3=zeros(length(F),length(FF));


%% read PREM
A = csvread('./PREM_1s.csv');
[depth_prem,aa]=unique(A(:,2));
mu_prem=A(aa,3)*1E3.*(A(aa,6)*1E3).^2;
rho_prem=A(aa,3)*1E3;
beta_prem=A(aa,6)*1E3;
alpha_prem=A(aa,4)*1E3;
ik=find(depth_prem<=50);
rho_prem(ik)=2900;
beta_prem(ik)=3300;
ik=find(depth_prem<=25);
rho_prem(ik)=2700;
beta_prem(ik)=3000;
alpha_prem(ik) = sqrt(3)*beta_prem(ik);


rho=2800;
alpha=6900;
beta=alpha/sqrt(3);
ERfac = (1/(15*pi*rho*alpha^5) +1/(10*pi*rho*beta^5) );
mu=rho*beta^2;


% depth vectors
HH=logspace(log10(5),log10(1200),11);
M0=logspace(18,23.5,11);
HH2=[0 50 100 1000];
NN=linspace(0,4,1000);
CC3=[0    0.4470    0.7410; ...
    0.8500    0.3250    0.0980; ...
    0.9290    0.6940    0.1250];


gaussian = @(x) (1/sqrt((2*pi))*exp(-x.^2/2));
skewedgaussian = @(x,alpha) 2*gaussian(x).*normcdf(alpha*x);

ske_best=mo;
agauss_best=mo;
resmin=mo;
agauss=linspace(1,10,100);
ske=linspace(-6,6,100);
tt3=(-3:dt2:3);
for ig=1:length(agauss)
   for isk=1:length(ske)
       Y{ig,isk}(1:length(tt2))=0;
        crap=skewedgaussian(tt3*agauss(ig), ske(isk));
        ik=find(crap>=0.01*max(crap));
        Y{ig,isk}(1:min(length(ik),length(tt2)))=crap(ik(1:min(length(ik),length(tt2))));
        sacc{ig,isk}=(diff(Y{ig,isk})/dt2).^2;sacc{ig,isk}(end+1)=0;
        es{ig,isk}=cumtrapz(sacc{ig,isk}.*dt2)./cumtrapz(Y{ig,isk}*dt2);es{ig,isk}=es{ig,isk}/es{ig,isk}(end);
        sacc{ig,isk}=sacc{ig,isk}/sum(sacc{ig,isk}(tt2>=0&tt2<=1)*dt2);
   end
end
      

% 
% for i=1:length(F)
% 
%     [a,x]=system(['ls ' char(DIR) '/FCTS/' char(F(i)) '/fctmoysource*']);
%     if a>0;continue;end
%     fid=fopen(char(x(1:end-1)),'r'); % read average STF file
%     fgetl(fid);
%     format='%f %f %f %f %f %f %f %f %f\n';
%     aa=textscan(fid,format ,'delimiter', '\n');
%     fclose('all');
%     xx=importdata(x(1:end-1));
%     
%     if 2/3*log10(xx(2,2))-6<6;continue;end % only keep M6+
%     dt=xx(5,1)-xx(4,1); % sampling rate of scardec
%     depth(i)=xx(2,1);   % extract depth 
%     tt=xx(4:end,1);     % extract time vector 
%     iik=find(xx(4:end,2)>0&tt<=0);% find indices in a standard time vector
%     ik1=find(tt1>=tt(1));
%     stf(i,ik1(1)+(1:length(xx(:,2))-3))=xx(4:end,2);% add STF into matrix of STFs
%    
%     
%     rake(i,1:2)=[aa{6} aa{9}];  % extract focal mechanism type
% %     use Shearer et al.(2006) vector representation of focal mechanism
% %     type
%     r1=rake(i,1); r2=rake(i,2);
%     if abs(r1)>90;r1=(180-abs(r1))*(r1/abs(r1));end
%     if abs(r2)>90;r2=(180-abs(r2))*(r2/abs(r2));end
%     if abs(r1)<abs(r2)
%       AMtype(i)=r1/90;
%     else
%         AMtype(i)=r2/90;
%     end
%     
% %     extract elastic properties from prem at earthquake depth
%     mu_p(i) = interp1(depth_prem,mu_prem,depth(i),'linear');
%     r_p(i) = interp1(depth_prem,rho_prem,depth(i),'linear');
%     b_p(i) = interp1(depth_prem,beta_prem,depth(i),'linear');    
%     a_p(i) = interp1(depth_prem,alpha_prem,depth(i),'linear');
% %     calculate energy coefficient at depth
%     ERfac_p(i) = (1/(15*pi*r_p(i)*a_p(i)^5) +1/(10*pi*r_p(i)*b_p(i)^5) );
%     
%     
%     %%  moments
%     mo(i)=xx(2,2);              % reported moment
%     mo1(i)=sum(stf(i,1:length(tt1))*dt1); % moment from integral of STF
%     mw(i)=2/3*log10(mo(i))-6;   % moment magnitude
%     
%     %% durations  
%     [~,icrap1]=find(stf(i,find(tt1>=0))>0);                 % take indices of positive time values
%     [~,icrap]=find(stf(i,find(tt1>=0))>0.1*max(stf(i,:)));  % pick index for threshold
%     TD(i) = length(icrap1(1):icrap(end))*dt ;               % duration of positive times until threshold
%     TM(i) = 2*sum(stf(i,1:length(tt1)).*tt1*dt1)/mo1(i);    % 2x centroid time
%     
%     
%     
%     % time to peak moment rate
%     [~,icrap]=max(abs(stf(i,:)));
%     Tpeak_stf(i) = icrap*dt1+tt1(1);                        % calculate the time-to-peak value
%     A_peak_stf(i) = abs(stf(i,icrap));                      % calculate peak values
%     
%     %% accelerations IMPORTANT CLARIFICATION:
%     % here, I remove spurious spikes that are due to SCARDEC discretization
%     % and run a running smoothing at those spikes.
%     acc=diff(stf(i,:))/dt1;                         % 1st order time derivative
%     er1(i,1:end-1) = ERfac*(acc).^2;er1(i,end)=0;   % seismic power
%     E0(i) = sum(er1(i,find(tt1>=0)))*dt1;            % energy measured from integrating the "rougher" functions.
%     for it=3:length(tt1)-3
%        if er1(i,it)>1.5*((er1(i,it-2)+er1(i,it-1)+er1(i,it+1)+er1(i,it+2))/4)
%            er1(i,it)=((er1(i,it-2)+er1(i,it-1)+er1(i,it+1)+er1(i,it+2))/4);
%        end
%     end
%     E(i) = sum(er1(i,find(tt1>=0)))*dt1;            % energy measured from integrating the smoothed functions.
%     [~,icrap]=find(er1(i,find(tt1>=0))>0.01*max(er1(i,:)));% measure duration of seismic power
%     TD3(i) = length(icrap1(1):icrap(end))*dt1 ;     % threshold for seismic power
%     II=find(tt1>=0);                                % take positive times
%     [~,icrap]=max(abs(er1(i,II)));                  % find max of seismic power
%     Tpeak_er(i) = icrap*dt1;                        % calculate the time-to-peak value
%     A_peak_er(i)=er1(i,II(icrap));                  % calculate peak values
%     
%     
%     %% choice of earthquake duration
%     TD5(i) = (TD(i)+TD3(i))/2;                      % average between STF and ER durations
%     
%     
%     %% Spectra    
%     crap=fft(stf(i,:))*dt1;freq2=1/2/dt1*linspace(0,1,floor(length(crap)/2));
%     xspec(i,:)=interp1(freq2,abs(crap(1:floor(length(crap)/2))),FF0,'linear'); % fft of stf
%     icrap2=find(isnan(abs(crap(1:floor(length(crap)/2))))==0&abs(crap(1:floor(length(crap)/2)))~=0);
%     ii=find(freq2(icrap2)<1);ii2=find(FF0<=1&isnan(xspec(i,:))==0);
%     Es(i) = 8*pi^2*ERfac*sum( abs(FF0(ii2).*xspec(i,ii2)).^2.*diff(FF0(ii2(1):ii2(end)+1)))/mo1(i); % energy under 1 Hz
%     
%     
%      if i==3165
%         figure
%         subplot(311)
%         yyaxis left
%         plot(tt1/TD5(i)*100,stf(i,:),'linewidth',2);hold on
% %         plot([1 1],[0 A_peak_stf(i)],'k');hold on
%         plot(Tpeak_stf(i)*[1 1]/TD5(i)*100,[0 A_peak_stf(i)],'k--','color',CC3(1,1:3),'linewidth',2);
%         plot(Tpeak_stf(i)*[0 1]/TD5(i)*100,A_peak_stf(i)*[1 1],'k--','color',CC3(1,1:3),'linewidth',2);
%         text(Tpeak_stf(i)/TD5(i)*100,1.05*A_peak_stf(i),'$A_M$','fontsize',14,'Interpreter','latex','color',CC3(1,1:3),'fontsize',14)
%         text(Tpeak_stf(i)/TD5(i)*100,0.1*A_peak_stf(i),'$T_M$','fontsize',14,'Interpreter','latex','color',CC3(1,1:3),'fontsize',14)
%         title('a)','fontsize',20);
%         xlabel('$t/T_D$ ','Interpreter','latex');
%         yyaxis right
%         plot(tt1/TD5(i)*100,cumtrapz(stf(i,:)*dt),'linewidth',2);xlim([0 110]);grid on;hold on
%         text(100,0.85*mo(i),'$M_0$','fontsize',14,'Interpreter','latex','color',CC3(2,1:3),'fontsize',14)
%         set(gca,'fontsize',14,'TickLabelInterpreter', 'latex');
% %         xlabel('time (s) ');
%         ylabel('$M_0(t)$','Interpreter','Latex')
%         yyaxis left; ylabel('$\dot{M}_0(t)$','Interpreter','Latex')
%         
%         subplot(312)
%         yyaxis left
%         plot(tt1/TD5(i)*100,er1(i,:),'linewidth',2);hold on;
% %         plot([1 1],[0 A_peak_er(i)],'k');hold on
%         plot(Tpeak_er(i)*[1 1]/TD5(i)*100,[0 A_peak_er(i)],'k--','color',CC3(1,1:3),'linewidth',2);grid on
%         plot(Tpeak_er(i)*[0 .5]/TD5(i)*100,A_peak_er(i)*[1 1],'k--','color',CC3(1,1:3),'linewidth',2);
%         text(Tpeak_er(i)/TD5(i)*100,A_peak_er(i),'$A_R$','fontsize',14,'Interpreter','latex','color',CC3(1,1:3),'fontsize',14)
%         text(Tpeak_er(i)/TD5(i)*100,0.1*A_peak_er(i),'$T_R$','fontsize',14,'Interpreter','latex','color',CC3(1,1:3),'fontsize',14)
%         title('b)','fontsize',20);xlim([0 110])
%         yyaxis right
%         plot(tt1/TD5(i)*100,cumtrapz(er1(i,:)*dt),'linewidth',2);grid on;hold on
%         text(100,0.85*E(i),'$E$','fontsize',14,'Interpreter','latex','color',CC3(2,1:3),'fontsize',14)
%         title('b)','fontsize',20);hold on
%         set(gca,'fontsize',14,'TickLabelInterpreter', 'latex');
%         xlabel('$t/T_D$ ','Interpreter','latex');
%         ylabel('$E_R(t)$','Interpreter','Latex')
%         yyaxis left; ylabel('$\dot{E}_R(t)$','Interpreter','Latex')
%         
%         
%         subplot(313)
%         yyaxis left
%         plot(tt1/TD5(i)*100,erstf(i,:),'linewidth',2);xlim([0 110]);grid on;hold on
% %         plot(TD5(i)/TD5(i)*[1 1],[0 A_peak_erstf(i)],'k');hold on
%         plot(Tpeak_erstf(i)*[1 1]/TD5(i)*100,[0 A_peak_erstf(i)],'k--','color',CC3(1,1:3),'linewidth',2);grid on
%         plot(Tpeak_erstf(i)*[0 .5]/TD5(i)*100,A_peak_erstf(i)*[1 1],'k--','color',CC3(1,1:3),'linewidth',2);
%         text(Tpeak_erstf(i)/TD5(i)*100,A_peak_erstf(i),'$A_S$','fontsize',14,'Interpreter','latex','color',CC3(1,1:3),'fontsize',14)
%         text(Tpeak_erstf(i)/TD5(i)*100,0.1*A_peak_erstf(i),'$T_S$','fontsize',14,'Interpreter','latex','color',CC3(1,1:3),'fontsize',14)
%         title('c)','fontsize',20);
%         set(gca,'fontsize',14,'TickLabelInterpreter', 'latex');
%         xlabel('$t/T_D $','Interpreter','latex');ylabel('$E_S(t)$','Interpreter','Latex')
%         text(100,E(i),'$E/M_0$','fontsize',14,'Interpreter','latex','color',CC3(2,1:3),'fontsize',14)
%         yyaxis right
%         set(gca,'YTick',[],'YTickLabel',[])
%         
%         print('-dpdf','-fillpage','./NEWFIGS/Fig1.pdf')
%         
%      end
%     
%      
%      
%     % spectral fit with 2 fc.
% %     resmin=9E19;
% %     for in=1:length(hffr)
% %     for ifreq1=1:length(ff1)-1
% %         for ifreq2=ifreq1:length(ff1)
% %             Ghat=mo1(i)./ sqrt(1 + (FF0(ii2)/ff1(ifreq1)).^hffr(in)) ./ sqrt(1 + (FF0(ii2)/ff1(ifreq2)).^hffr(in));
% %             res1=sum((log10( xspec(i,ii2) )  - log10(Ghat)).^2);
% %             if res1<=resmin
% %                 resmin=res1;
% %                 fbest1(i)=ff1(ifreq1);
% %                 NNbest(i)=hffr(in);
% %                 fbest2(i)=ff1(ifreq2);
% %             end
% %         end
% %     end
% %     end
% %     Ghat = mo1(i)./ sqrt(1 + (FF0/fbest1(i)).^NNbest(i)) ./ sqrt(1 + (FF0/fbest2(i)).^NNbest(i));
% %     figure(12)
% %     loglog(FF0(ii2),xspec(i,ii2),FF0,Ghat);grid on;
%     %% energy from the spectra
%     % amount of energy missed in the finite frequency cut off.
% %     ii3=find(FF0>1);
% %     Ecorr(i) =  8*pi^2*ERfac*sum( abs(FF0(ii3).*Ghat(ii3)).^2.*diff(FF0(ii3(1)-1:ii3(end))))/mo1(i); % correction term for energy
% %     disp('Error in unreliable energy after 1Hz.')
% %     disp([Es(i)/E(i)*mo1(i) Ecorr(i)/Es(i)])
% %     Es(i) = Es(i)+Ecorr(i);
%     
% 
%     %% NOTE HERE: OBVIOUS LOSS OF ENERGY FOR SMALL EVENTS. I CHOSE TO NOT ADD A CORRECTION BECAUSE
%     %% I ONLY SHOW FR / ER IN THE PAPER (INDEPENDENT OF CORRECTION). ADDING A CORRECTION DOES NOT CHANGE THE RESULTS, CONTACT ME
%     %% IF YOU NEED TO DISCUSS THIS: mdenolle@fas.harvard.edu
% 
%     %% time-lapse scaled energy     
%     erstf(i,:)=cumtrapz(er1(i,:)*dt1)./cumtrapz(stf(i,:)*dt1);
% %     erstf2(i,:)=cumsum(er1(i,:)*dt1)./cumsum(stf(i,:)*dt1);
%     
% %     plot(tt1,erstf(i,:),tt1,erstf2(i,:),'linewidth',2);grid on;xlim([0 TD5(i)]);pause
%     II=find(isnan(erstf(i,:))==0&isinf(erstf(i,:))==0&tt1>=-0.05);
%     [~,icrap]=max(erstf(i,II));                     % find max of time-lapse scaled energy
%     Tpeak_erstf(i) = icrap*dt1+tt1(II(1));          % calculate time-to-peak value
%     A_peak_erstf(i)=max(erstf(i,II));               % calculate peak value
%     Es_final_test(i)=erstf(i,II(end));               % calculate peak value
%     
%     
%     
%     %% Normalizing the functions
%     % with duration
%     stf_sd(i,1:length(tt2)) = interp1(tt1/TD5(i),stf(i,:)*TD5(i)/mo1(i),tt2,'linear'); % STF
%     er_sd(i,:) = interp1(tt1/TD5(i),er1(i,:),tt2,'linear');                             % energy
%     II1=find(isnan(er_sd(i,:))==0&isinf(er_sd(i,:))==0);
%     Esum=sum(er_sd(i,II1)*dt2); 
%     er_sd(i,II1)=er_sd(i,II1)/Esum;                   % normalize amplitude
%     
%     erstf_sd(i,:) = interp1(tt1(II)/TD5(i),erstf(i,II),tt2,'linear');               % scaled energy
%     II=find(isnan(erstf_sd(i,:))==0&isinf(erstf_sd(i,:))==0);
%     Esum=sum(erstf_sd(i,II)*dt2);   
%     erstf_sd(i,II)=erstf_sd(i,II)/Es(i) ;             % normalize amplitude
%     
%     
%     % with 2x centroid moment
%     stf_sm(i,1:length(tt2)) = interp1(tt1/TM(i),stf(i,:)*TM(i)/mo1(i),tt2,'linear'); % STF
%     er_sm(i,:) = interp1(tt1/TM(i),er1(i,:),tt2,'linear');                             % energy
%     II1=find(isnan(er_sm(i,:))==0&isinf(er_sm(i,:))==0);
%     Esum=sum(er_sm(i,II1)*dt2); 
%     er_sm(i,II1)=er_sm(i,II1)/Esum;                   % normalize amplitude
%     
%     erstf_sm(i,:) = interp1(tt1(II)/TM(i),erstf(i,II),tt2,'linear');               % scaled energy
%     II=find(isnan(erstf_sm(i,:))==0&isinf(erstf_sm(i,:))==0);
%     Esum=sum(erstf_sm(i,II)*dt2);   
%     erstf_sm(i,II)=erstf_sm(i,II)/Es(i) ;             % normalize amplitude to total energy
%     
%     
%     
%     % fit the gaussian for simple models.
%     resmin(i)=9E9; 
%     II=find(tt2>=0&tt2<=1);   
%     for ig=1:length(agauss)
%     for isk=1:length(ske)
%         res1 = norm([erstf_sd(i,II)-es{ig,isk}(II) ...
%             er_sd(i,II)-sacc{ig,isk}(II) stf_sd(i,II)-Y{ig,isk}(II)],2);
%         if res1<=resmin(i)
%             resmin(i)=res1;
%             agauss_best(i)=agauss(ig);
%             ske_best(i)=ske(isk);
%         end
%     end
%     end
%     if abs(ske_best(i))==4||agauss_best(i)==min(agauss)||agauss_best(i)==min(agauss)
%         fit_skG(i)=0;
%     else
%         fit_skG(i)=1;
%         disp([agauss_best(i) ske_best(i)])
%         
%         
% %         crap=skewedgaussian(tt3*agauss_best(i), ske_best(i));
% %         ik=find(crap>=0.01*max(crap));
% %         Y1=zeros(size(tt2));
% %         Y1(1:min(length(ik),length(tt2)))=crap(ik(1:min(length(ik),length(tt2))));
% %         Y1=Y1/sum(Y1*dt2);
% %         sacc1=(diff(Y1)/dt2).^2;sacc1(end+1)=0;
% %         es1=cumtrapz(sacc1.*dt2)./cumtrapz(Y1*dt2);es1=es1/es1(end);
% %         sacc1=sacc1/sum(sacc1(tt2>=0&tt2<=1)*dt2);
% %         crap=fft(Y);n2=floor(length(crap)/2);freq1=1/2/dt2*linspace(0,1,n2);
% %         xmod(i,:) = interp1(freq1,abs(crap(1:n2)),FF0,'linear');
% 
% %         figure(1)
% %         subplot(311)
% %         plot(tt2,stf_sd(i,:),tt2,Y1,'linewidth',2);grid on;xlim([0 1.1])
% %         title(['mw ' num2str(mw(i)) 'skewness ' num2str(ske_best(i))])
% %         subplot(312)
% %         plot(tt2,er_sd(i,:),tt2,sacc1,'linewidth',2);grid on;xlim([0 1.1])
% %         subplot(313)
% %         plot(tt2,erstf_sd(i,:),tt2,es1,'linewidth',2);grid on;xlim([0 1.1]);legend('data','model')
% % %         subplot(414)
% % %         loglog(FF0,xmod(i,:),'linewidth',2);grid on;xlim([0.5 10])
% %         pause(0.01)
%     end
%     
% %     % skewness until the 2 metrics of duration
% %     sKd(i) = skewness(stf(i,tt1>=0&tt1<=TD5(i)));
% %     sKm(i) = skewness(stf(i,tt1>=0&tt1<=TM(i)));
% %     % kurtosis: impulsivity
% %     Kd(i) = kurtosis(stf(i,tt1>=0&tt1<=TD5(i)));
% %     Km(i) = kurtosis(stf(i,tt1>=0&tt1<=TM(i)));
%     
%    
%     %% show normalize function spectra
%     crap=fft(stf_sd(i,:))*dt2;freq2=1/2/dt2*linspace(0,1,floor(length(crap)/2));
%     xspec2(i,:)=interp1(freq2,abs(crap(1:floor(length(crap)/2))),FF,'linear'); % fft of stf
%     
%     crap=fft(stf_sm(i,:))*dt2;freq2=1/2/dt2*linspace(0,1,floor(length(crap)/2));
%     xspec3(i,:)=interp1(freq2,abs(crap(1:floor(length(crap)/2))),FF,'linear'); % fft of stf
%     
%     
% end
% save scardec_moy_M6_nofit.mat
load scardec_moy_M6_nofit.mat

HH2=[0 50 350 1000];
tgrowth = min([Tpeak_stf   Tpeak_er Tpeak_erstf ],[],2);

[nn,~]=find(mw~=0&tgrowth>=1);  % select earthquakes where tgrowth>=1)
%  [nn,~]=find(mw~=0);  % select earthquakes where T_S>1)


% a bunch of colormaps 
CC=colormap(parula(100));   
CC2=colormap(jet(100));
CC4=colormap(summer(100));
Ch=linspace(0.5,3,100);
Cm=linspace(6.5,9,100);
C1=linspace(0,1,100);
C2=linspace(-1,1,100);
Cg=linspace(0,0.75,100);

% sort depth, magnitude, time-to-peak value
[~,aa]=sort(depth(nn));
[~,aa1]=sort(mo(nn));
[~,aa4]=sort(Tpeak_erstf(nn)./TD5(nn));


% 
for i=856%1:3165;
if mw(i)<7||ske_best(i)==0;continue;end
crap=skewedgaussian(tt3*agauss_best(i), ske_best(i));
ik=find(crap>=0.01*max(crap));
Y1=zeros(size(tt2));
Y1(1:min(length(ik),length(tt2)))=crap(ik(1:min(length(ik),length(tt2))));
Y1=Y1/sum(Y1*dt2);
sacc1=(diff(Y1)/dt2).^2;sacc1(end+1)=0;
es1=cumtrapz(sacc1.*dt2)./cumtrapz(Y1*dt2);es1=es1/es1(end);
sacc1=sacc1/sum(sacc1(tt2>=0&tt2<=1)*dt2);

figure
% subplot(311)
% plot(tt2,stf_sd(i,:),tt2,Y1,'linewidth',2);grid on;xlim([0 1.1])
% title(['mw ' num2str(mw(i)) 'skewness ' num2str(ske_best(i))])
% subplot(312)
% plot(tt2,er_sd(i,:),tt2,sacc1,'linewidth',2);grid on;xlim([0 1.1])
% subplot(313)
% plot(tt2,erstf_sd(i,:),tt2,es1,'linewidth',2);grid on;xlim([0 1.1]);legend('data','model')
% %         subplot(414)
% %         loglog(FF0,xmod(i,:),'linewidth',2);grid on;xlim([0.5 10])
% pause(0.01)
% figure
subplot(311)
yyaxis left
plot(tt1/TD5(i)*100,stf(i,:),tt2*100,Y1*mo1(i)/TD5(i),'linewidth',2);hold on
plot(Tpeak_stf(i)*[1 1]/TD5(i)*100,[0 A_peak_stf(i)],'k-.','color',[128,0,128]/255,'linewidth',2);
plot(Tpeak_stf(i)*[0 1]/TD5(i)*100,A_peak_stf(i)*[1 1],'k-.','color',[128,0,128]/255,'linewidth',2);
text(Tpeak_stf(i)/TD5(i)*100,1.05*A_peak_stf(i),'$A_M$','fontweight','bold','fontsize',14,'Interpreter','latex','color',[128,0,128]/255,'fontsize',14)
text(Tpeak_stf(i)/TD5(i)*100,0.1*A_peak_stf(i),'$T_M$','fontweight','bold','fontsize',14,'Interpreter','latex','color',[128,0,128]/255,'fontsize',14)
title('a)','fontsize',20);ylim([0 5.5E19])
xlabel('$t/T (\%)$ ','Interpreter','latex');
yyaxis right;box on
plot(tt1/TD5(i)*100,cumtrapz(stf(i,:)*dt),'linewidth',2);xlim([0 110]);grid on;hold on
text(100,0.85*mo(i),'$M_0$','fontsize',14,'Interpreter','latex','color',CC3(2,1:3),'fontsize',14)
set(gca,'fontsize',14,'TickLabelInterpreter', 'latex');
ylabel('$M_0(t)$','Interpreter','Latex')
yyaxis left; ylabel('$\dot{M}_0(t)$','Interpreter','Latex')



subplot(312)
yyaxis left
plot(tt1/TD5(i)*100,er1(i,:),tt2*100,sacc1*E(i)/TD5(i),'linewidth',2);hold on;
%         plot([1 1],[0 A_peak_er(i)],'k');hold on
plot(Tpeak_er(i)*[1 1]/TD5(i)*100,[0 A_peak_er(i)],'k-.','color',[128,0,128]/255,'linewidth',2);grid on
plot(Tpeak_er(i)*[0 1]/TD5(i)*100,A_peak_er(i)*[1 1],'k-.','color',[128,0,128]/255,'linewidth',2);
text(Tpeak_er(i)/TD5(i)*100,A_peak_er(i),'$A_R$','fontweight','bold','fontsize',14,'Interpreter','latex','color',[128,0,128]/255,'fontsize',14)
text(Tpeak_er(i)/TD5(i)*100,0.1*A_peak_er(i),'$T_R$','fontweight','bold','fontsize',14,'Interpreter','latex','color',[128,0,128]/255,'fontsize',14)
title('b)','fontsize',20);xlim([0 110])
yyaxis right;box on
plot(tt1/TD5(i)*100,cumtrapz(er1(i,:)*dt),'linewidth',2);grid on;hold on
text(100,0.85*E(i),'$E$','fontsize',14,'Interpreter','latex','color',CC3(2,1:3),'fontsize',14)
title('b)','fontsize',20);hold on
set(gca,'fontsize',14,'TickLabelInterpreter', 'latex');
xlabel('$t/T (\%)$ ','Interpreter','latex');
ylabel('$E_R(t)$','Interpreter','Latex')
yyaxis left; ylabel('$\dot{E}_R(t)$','Interpreter','Latex')


subplot(313)
yyaxis left
plot(tt1/TD5(i)*100,erstf(i,:),tt2*100,es1*E(i)/mo1(i),'linewidth',2);xlim([0 110]);grid on;hold on
%         plot(TD5(i)/TD5(i)*[1 1],[0 A_peak_erstf(i)],'k');hold on
plot(Tpeak_erstf(i)*[1 1]/TD5(i)*100,[0 A_peak_erstf(i)],'k-.','color',[128,0,128]/255,'linewidth',2);grid on
plot(Tpeak_erstf(i)*[0 .5]/TD5(i)*100,A_peak_erstf(i)*[1 1],'k-.','color',[128,0,128]/255,'linewidth',2);
text(Tpeak_erstf(i)/TD5(i)*100,A_peak_erstf(i),'$A_S$','fontweight','bold','fontsize',14,'Interpreter','latex','color',[128,0,128]/255,'fontsize',14)
text(Tpeak_erstf(i)/TD5(i)*100,0.1*A_peak_erstf(i),'$T_S$','fontweight','bold','fontsize',14,'Interpreter','latex','color',[128,0,128]/255,'fontsize',14)
title('c)','fontsize',20);
set(gca,'fontsize',14,'TickLabelInterpreter', 'latex');
xlabel('$t/T (\%) $','Interpreter','latex');ylabel('$E_S(t)$','Interpreter','Latex')
text(100,E(i),'$E/M_0$','fontsize',14,'Interpreter','latex','color',CC3(2,1:3),'fontsize',14)
yyaxis right;box on
set(gca,'YTick',[],'YTickLabel',[])

disp(i)
pause
print('-dpdf','-fillpage','./NEWFIGS/Fig1.pdf')


end
% 14
% 54 (ya)
% 223
% 562
% 610 (small quake?)
% 641
% 724
% 780
% 795
% 821 
% 
% 856
%% correct tohoku: ONE OBVIOUS OUTLIER IN THIS METHOD IS TOHOKU. CORRECT FOR TRUE DURATION ~ 150 s.
[~,amax]=max(mw(nn));
[~,icrap]=find(stf(nn(amax),find(tt1>=0))>0);
TD5(nn(amax))=length(icrap)*dt1;

%% make regression to all parameters
AM=A_peak_stf(nn)./(r_p(nn).^(1/3).*b_p(nn).^(5/3))...
    .*rho^(1/3).*beta^(5/3);    % you can check that this confirms Vallee 2013 .
AR=A_peak_er(nn)./(mo1(nn).*Es(nn));     % peak amplitude of seismic power over radiated energy
AR1=(A_peak_er(nn) + Ecorr(nn)./TD5(nn))./(mo1(nn).*(Es(nn)+Ecorr(nn)));     % peak amplitude of seismic power over radiated energy
AS=((A_peak_erstf(nn)))./(E(nn)./mo1(nn));       % peak scaled energy over total scaled energy. 
                    %effectivelty this is the peak amplitude over the last amplitude of the function
    
%% make regression to see what variable best reduces the variance.
AR_MO=fitlm(log10(mo(nn)),log10(AR),'linear');
T_MO=fitlm(log10(mo(nn)),log10(TD5(nn)),'linear');
AR_MO_DP=fitlm(log10(depth(nn)),log10(AR.*mo(nn).^(-AR_MO.Coefficients.Estimate(2))),'linear');
AS_TS=fitlm(log10(Tpeak_stf(nn)./TD5(nn)),log10(AS),'linear');


close all
%% Fig S1: variations of duration with H, Mo, FM
figure
crap=TD5(nn)./mo(nn).^(1/3).*r_p(nn).^(1/3).*b_p(nn).^(5/3)/(rho^(1/3)*beta^(5/3));
set(0,'defaulttextinterpreter','latex')
subplot(311)
s=scatter(log10(depth(nn)), log10(crap),'filled');grid on;hold on;
    s.MarkerFaceAlpha=0.4;
    s.MarkerEdgeColor=0.5*[1 1 1];
    s.MarkerFaceColor=[173,216,230]/255;
set(gca,'fontsize',14);box on
set(gca,'TickLabelInterpreter', 'latex');
title('a)','position',[0.6 -4.8],'fontsize',20);
xlabel('$log_{10}$ H (km) ');ylabel('$\bar{T}$','Interpreter','Latex')

subplot(312)
s=scatter(log10(mo(nn)), log10(crap),'filled');grid on;hold on;
    s.MarkerFaceAlpha=0.4;
    s.MarkerEdgeColor=0.5*[1 1 1];
    s.MarkerFaceColor=[173,216,230]/255;
set(gca,'fontsize',14);box on
title('b)','position',[18.1 2E-5],'fontsize',20)
set(gca,'TickLabelInterpreter', 'latex');
title('b)','position',[18.2 -4.8],'fontsize',20);
xlabel('$log_{10} M_0$ (Nm)');ylabel('$\bar{T}$','Interpreter','Latex')
subplot(313)
s=scatter(AMtype(nn), log10(crap),'filled');grid on;hold on;
    s.MarkerFaceAlpha=0.4;
    s.MarkerEdgeColor=0.5*[1 1 1];
    s.MarkerFaceColor=[173,216,230]/255;
set(gca,'fontsize',14);box on
set(gca,'TickLabelInterpreter', 'latex');
title('c)','position',[-0.95 -4.8],'fontsize',20);
xlabel('Focal mechanism type');ylabel('$\bar{T}$','Interpreter','Latex')
xlim([-1 1])
print('-dpdf','-fillpage','./NEWFIGS/FigS1.pdf')




%% Fig S2 = moment dependence of TM, TR, TS
figure
set(0,'defaulttextinterpreter','latex')
subplot(311)
crap=Tpeak_stf(nn)./TD5(nn)*100;
for i=1:length(nn)
    [~,ik]=min(abs(log10(depth(nn(aa(i))))-Ch));
    plot(log10(mo(nn(aa(i)))),crap(aa(i)),'bo','markeredgecolor',0.5*[1 1 1],'markerfacecolor',CC(ik,1:3),...
        'markersize',8);
    hold on; grid on;
end
set(gca,'fontsize',14);box on
title('a)','position',[18.1 105],'fontsize',20);
set(gca,'TickLabelInterpreter', 'latex');ylim([0 100])
xlabel('$log_{10} M_0 $ (Nm) ');ylabel('$T_M (\%)$','Interpreter','Latex')
subplot(312)
crap=Tpeak_er(nn)./TD5(nn)*100;
for i=1:length(nn)
    [~,ik]=min(abs(log10(depth(nn(aa(i))))-Ch));
    plot(log10(mo(nn(aa(i)))),crap(aa(i)),'bo','markeredgecolor',0.5*[1 1 1],'markerfacecolor',CC(ik,1:3),...
        'markersize',8);
    hold on; grid on;
end
set(gca,'fontsize',14);box on
title('b)','position',[18.1 105],'fontsize',20);
set(gca,'TickLabelInterpreter', 'latex');ylim([0 100])
xlabel('$log_{10} M_0$ (Nm)');ylabel('$T_R  (\%)$','Interpreter','Latex')
subplot(313)
crap=Tpeak_erstf(nn)./TD5(nn)*100;
for i=1:length(nn)
    [~,ik]=min(abs(log10(depth(nn(aa(i))))-Ch));
    plot(log10(mo(nn(aa(i)))),crap(aa(i)),'bo','markeredgecolor',0.5*[1 1 1],'markerfacecolor',CC(ik,1:3),...
        'markersize',8);
    hold on; grid on;
end
set(gca,'fontsize',14);ylim([0 100]);box on
title('c)','position',[18.1 105],'fontsize',20);
set(gca,'TickLabelInterpreter', 'latex');
xlabel('$log_{10} M_0$ (Nm)');ylabel('$T_S  (\%)$','Interpreter','Latex')
print('-dpdf','-fillpage','./NEWFIGS/FigS2.pdf')

close all
%% figure S3: Tr, Tm, Ts as a function of focal mechanism
AAtype=(-1:0.2:1);
figure
subplot(311)
s=scatter(AMtype(nn), Tpeak_stf(nn)./TD5(nn)*100,'filled');grid on;hold on;
    s.MarkerFaceAlpha=0.4;
    s.MarkerEdgeColor=0.5*[1 1 1];
    s.MarkerFaceColor=[173,216,230]/255;
crap=Tpeak_stf./TD5*100;
for ih=1:length(AAtype)-1
    [ik,~]=find(AMtype(nn)>=AAtype(ih)&AMtype(nn)<AAtype(ih+1));
    if isempty(ik);continue;end
    hh=AAtype(ih)+0.1;
    % bootstrap median
    medT=zeros(100,1);
    for k=1:1000
        kk=datasample((1:length(ik)),length(ik));
        medT(k)=median(crap(nn(ik(kk))));
    end
    medB=mean(medT);stdB=std(medT);
    stdT=std(crap(nn(ik)));
    plot(hh+[-0.1 0.1]/2,(medB-stdT)*[1 1],'k','linewidth',1.5);
    plot(hh+[-0.1 0.1]/2,(medB+stdT)*[1 1],'k','linewidth',1.5);
    plot(hh*[1 1],medB+[-stdT stdT],'k','linewidth',1.5);
    plot(hh+[-0.1 0.1]/2,(medB-stdB)*[1 1],'r','linewidth',1.5);
    plot(hh+[-0.1 0.1]/2,(medB+stdB)*[1 1],'r','linewidth',1.5);
    plot(hh*[1 1],medB+[-stdB stdB],'k','linewidth',1.5);
    plot(hh,medB,'ks','MarkerFaceColor',[255,99,71]/255,'markersize',14);
end
set(gca,'fontsize',14);box on
title('a)','position',[105 -0.75],'fontsize',20);ylim([0 100])
set(gca,'TickLabelInterpreter', 'latex');
xlabel('fm type ');ylabel('$T_M (\%)$','Interpreter','Latex')

subplot(312)
s=scatter(AMtype(nn), Tpeak_er(nn)./TD5(nn)*100,'filled');grid on;hold on;
    s.MarkerFaceAlpha=0.4;
    s.MarkerEdgeColor=0.5*[1 1 1];
    s.MarkerFaceColor=[173,216,230]/255;
crap=Tpeak_er./TD5*100;
for ih=1:length(AAtype)-1
    [ik,~]=find(AMtype(nn)>=AAtype(ih)&AMtype(nn)<AAtype(ih+1));
    if isempty(ik);continue;end
    hh=AAtype(ih)+0.1;
    % bootstrap median
    medT=zeros(100,1);
    for k=1:1000
        kk=datasample((1:length(ik)),length(ik));
        medT(k)=median(crap(nn(ik(kk))));
    end
    medB=mean(medT);stdB=std(medT);
    stdT=std(crap(nn(ik)));
    plot(hh+[-0.1 0.1]/2,(medB-stdT)*[1 1],'k','linewidth',1.5);
    plot(hh+[-0.1 0.1]/2,(medB+stdT)*[1 1],'k','linewidth',1.5);
    plot(hh*[1 1],medB+[-stdT stdT],'k','linewidth',1.5);
    plot(hh+[-0.1 0.1]/2,(medB-stdB)*[1 1],'r','linewidth',1.5);
    plot(hh+[-0.1 0.1]/2,(medB+stdB)*[1 1],'r','linewidth',1.5);
    plot(hh*[1 1],medB+[-stdB stdB],'k','linewidth',1.5);
    plot(hh,medB,'ks','MarkerFaceColor',[255,99,71]/255,'markersize',14);
end
set(gca,'fontsize',14);box on
title('b)','position',[105 -0.75],'fontsize',20);ylim([0 100])
set(gca,'TickLabelInterpreter', 'latex');
xlabel('fm type ');ylabel('$T_R (\%)$','Interpreter','Latex')
subplot(313)
crap=Tpeak_erstf./TD5*100;
s=scatter(AMtype(nn), Tpeak_erstf(nn)./TD5(nn)*100,'filled');grid on;hold on;
    s.MarkerFaceAlpha=0.4;
    s.MarkerEdgeColor=0.5*[1 1 1];
    s.MarkerFaceColor=[173,216,230]/255;
for ih=1:length(AAtype)-1
    [ik,~]=find(AMtype(nn)>=AAtype(ih)&AMtype(nn)<AAtype(ih+1));
    if isempty(ik);continue;end
    hh=AAtype(ih)+0.1;
    % bootstrap median
    medT=zeros(100,1);
    for k=1:1000
        kk=datasample((1:length(ik)),length(ik));
        medT(k)=median(crap(nn(ik(kk))));
    end
    medB=mean(medT);stdB=std(medT);
    stdT=std(crap(nn(ik)));
    plot(hh+[-0.1 0.1]/2,(medB-stdT)*[1 1],'k','linewidth',1.5);
    plot(hh+[-0.1 0.1]/2,(medB+stdT)*[1 1],'k','linewidth',1.5);
    plot(hh*[1 1],medB+[-stdT stdT],'k','linewidth',1.5);
    plot(hh+[-0.1 0.1]/2,(medB-stdB)*[1 1],'r','linewidth',1.5);
    plot(hh+[-0.1 0.1]/2,(medB+stdB)*[1 1],'r','linewidth',1.5);
    plot(hh*[1 1],medB+[-stdB stdB],'k','linewidth',1.5);
    plot(hh,medB,'ks','MarkerFaceColor',[255,99,71]/255,'markersize',14);
end
set(gca,'fontsize',14);box on
title('c)','position',[105 -0.75],'fontsize',20);ylim([0 100])
set(gca,'TickLabelInterpreter', 'latex');
xlabel('fm type ');ylabel('$T_S (\%)$','Interpreter','Latex')
print('-dpdf','-fillpage','./NEWFIGS/FigS3.pdf')



%% Fig S4
% is energy coming from acceleration or deceleration?
figure                               
 subplot(211)
 hist(log10(Tpeak_er(nn)./Tpeak_stf(nn)),100);set(gca,'TickLabelInterpreter', 'latex');
 xlabel('$\log_{10}\left(T_R/T_M\right)$','Interpreter','latex');set(gca,'Fontsize',14);       
 xlim([-2 2])
 grid on
 subplot(212)
 hist(log10(Tpeak_erstf(nn)./Tpeak_stf(nn)),100);set(gca,'TickLabelInterpreter', 'latex');
 xlabel('$\log_{10}\left(T_S/T_M\right)$','Interpreter','latex');set(gca,'Fontsize',14); 
 xlim([-2 2])
 grid on
set(gcf, 'renderer', 'painters')
 print('-dpdf','-fillpage','./NEWFIGS/FigS4.pdf')

close all


figure
s=scatter(AMtype(nn), AS,'filled');grid on;hold on;
    s.MarkerFaceAlpha=0.4;
    s.MarkerEdgeColor=0.5*[1 1 1];
    s.MarkerFaceColor=[173,216,230]/255;
for ih=1:length(AAtype)-1
    [ik,~]=find(AMtype(nn)>=AAtype(ih)&AMtype(nn)<AAtype(ih+1));
    if isempty(ik);continue;end
    hh=AAtype(ih)+0.1;
    % bootstrap median
    medT=zeros(100,1);
    for k=1:1000
        kk=datasample((1:length(ik)),length(ik));
        medT(k)=median(AS(ik(kk)));
    end
    medB=mean(medT);stdB=std(medT);
    stdT=std(AS(ik));
%     plot(hh+[-0.1 0.1]/2,(medB-stdT)*[1 1],'k','linewidth',1.5);
%     plot(hh+[-0.1 0.1]/2,(medB+stdT)*[1 1],'k','linewidth',1.5);
    plot(hh*[1 1],medB+[-stdB stdB],'k','linewidth',1.5);
    plot(hh+[-0.1 0.1]/2,(medB-stdB)*[1 1],'r','linewidth',1.5);
    plot(hh+[-0.1 0.1]/2,(medB+stdB)*[1 1],'r','linewidth',1.5);
    plot(hh*[1 1],medB+[-stdB stdB],'k','linewidth',1.5);
    plot(hh,medB,'ks','MarkerFaceColor',[255,99,71]/255,'markersize',14);
end
set(gca,'fontsize',14);box on
ylim([0 15])
set(gca,'TickLabelInterpreter', 'latex');
xlabel('fm type ');ylabel('$A_S $','Interpreter','Latex')
print('-dpdf','-fillpage','./NEWFIGS/FigS5.pdf')


%% figure S4: TD vd Tm
figure
for i=1:length(nn)
    [~,ik]=min(abs(log10(depth(nn(aa(i))))-Ch));
    loglog(TM(nn(aa(i))),TD5(nn(aa(i))),'bo','markeredgecolor',0.5*[1 1 1],'markerfacecolor',CC(ik,1:3),...
        'markersize',8);
    hold on; grid on;
end
loglog([1 200],[1 200],'r');xlim([1 200]);ylim([1 200]);box on
text(30,10,['std of pdf ' num2str((10^std(log10(TM(nn)./TD5(nn)))-1)*100) '\%'],'fontsize',14)
xlabel('$2T_M$ (s)','Interpreter','Latex');ylabel('T (s)','Interpreter','Latex');
set(gca,'fontsize',14)
set(gcf, 'renderer', 'painters')
print('-dpdf','./NEWFIGS/FigS6.pdf');








%% Fig 2 of paper
figure
subplot(311)
s=scatter(log10(depth(nn)), Tpeak_stf(nn)./TD5(nn)*100,'filled');grid on;hold on;box on
    s.MarkerFaceAlpha=0.4;
    s.MarkerEdgeColor=0.5*[1 1 1];
    s.MarkerFaceColor=[173,216,230]/255;
    s.SizeData=100;
for ih=1:length(HH)-1
    [ik,~]=find(depth(nn)>=HH(ih)&depth(nn)<HH(ih+1));
    if isempty(ik);continue;end
    hh=(mean(log10(HH(ih:ih+1))));
    % bootstrap medianmediN
    medT=zeros(100,1);
    for k=1:1000
        kk=datasample((1:length(ik)),length(ik));
        medT(k)=mean(Tpeak_stf(nn(ik(kk)))./TD5(nn(ik(kk)))*100);
    end
    medB=mean(medT);stdB=std(medT);
    stdT=std(Tpeak_stf(nn(ik))./TD5(nn(ik))*100);
    plot(hh+[-0.1 0.1]/2,(medB-stdT)*[1 1],'k','linewidth',1.5);
    plot(hh+[-0.1 0.1]/2,(medB+stdT)*[1 1],'k','linewidth',1.5);
    plot(hh*[1 1],medB+[-stdT stdT],'k','linewidth',1.5);
    plot(hh+[-0.1 0.1]/2,(medB-stdB)*[1 1],'k','color',[255,99,71]/255,'linewidth',1.5);
    plot(hh+[-0.1 0.1]/2,(medB+stdB)*[1 1],'k','color',[255,99,71]/255,'linewidth',1.5);
    plot(hh*[1 1],medB+[-stdB stdB],'k','linewidth',1.5);
    plot(hh,medB,'ks','MarkerFaceColor',[255,99,71]/255,'markersize',13);
end
set(gca,'fontsize',14);ylabel('$T_M$ (\%)','Interpreter','Latex');xlabel('$\log_{10}$ H (km)','Interpreter','Latex');ylim([0 100])
title('a)','position',[0.55 110],'fontsize',20);
set(gca,'TickLabelInterpreter', 'latex');
subplot(312)
s=scatter(log10(depth(nn)), Tpeak_er(nn)./TD5(nn)*100,'filled');grid on;hold on;box on
    s.MarkerFaceAlpha=0.4;
    s.MarkerEdgeColor=0.5*[1 1 1];
    s.MarkerFaceColor=[173,216,230]/255;
    s.SizeData=100;
set(gca,'fontsize',14);ylabel('$T_R$ (\%)','Interpreter','Latex');xlabel('$\log_{10}$ H (km)','Interpreter','Latex');ylim([0 100])
title('b)','position',[0.55 110],'fontsize',20);box on
set(gca,'TickLabelInterpreter', 'latex');
subplot(313)
s=scatter(log10(depth(nn)), Tpeak_erstf(nn)./TD5(nn)*100,'filled');grid on;hold on;
    s.MarkerFaceAlpha=0.4;
    s.MarkerEdgeColor=0.5*[1 1 1];
    s.MarkerFaceColor=[173,216,230]/255;
    s.SizeData=100;
for ih=1:length(HH)-1
    [ik,~]=find(depth(nn)>=HH(ih)&depth(nn)<HH(ih+1));
    if isempty(ik);continue;end
    hh=(mean(log10(HH(ih:ih+1))));
    % bootstrap median
    medT=zeros(100,1);
    for k=1:1000
        kk=datasample((1:length(ik)),length(ik));
        medT(k)=median(Tpeak_erstf(nn(ik(kk)))./TD5(nn(ik(kk)))*100);
    end
    medB=median(medT);stdB=std(medT);
    stdT=std(Tpeak_erstf(nn(ik))./TD5(nn(ik))*100);
    plot(hh+[-0.1 0.1]/2,(medB-stdT)*[1 1],'k','linewidth',1.5);
    plot(hh+[-0.1 0.1]/2,(medB+stdT)*[1 1],'k','linewidth',1.5);
    plot(hh*[1 1],medB+[-stdT stdT],'k','linewidth',1.5);
    plot(hh+[-0.1 0.1]/2,(medB-stdB)*[1 1],'k','color',[255,99,71]/255,'linewidth',1.5);
    plot(hh+[-0.1 0.1]/2,(medB+stdB)*[1 1],'k','color',[255,99,71]/255,'linewidth',1.5);
    plot(hh*[1 1],medB+[-stdB stdB],'k','linewidth',1.5);
    plot(hh,medB,'ks','MarkerFaceColor',[255,99,71]/255,'markersize',13);
end
box on
set(gca,'fontsize',14);ylabel('$T_S$ (\%)','Interpreter','Latex');xlabel('$\log_{10}$ H (km)','Interpreter','Latex');ylim([0 100])
title('c)','position',[0.55 110],'fontsize',20);
set(gca,'TickLabelInterpreter', 'latex');box on
set(gcf, 'renderer', 'painters')
print('-dpdf','-fillpage','./NEWFIGS/Fig2a.pdf')


%% part II of Fig 2 of paper
figure
subplot(311)
[h,N]=hist(Tpeak_stf(nn)./TD5(nn)*100,50);
Bh=bar(N,h,'facecolor',[143,188,143]/255);xlim([0 100])
set(gca,'xtick',[],'xticklabel',[],'ytick',[],'yticklabel',[])
subplot(312)
[h,N]=hist(Tpeak_er(nn)./TD5(nn)*100,50);
Bh=bar(N,h,'facecolor',[143,188,143]/255);xlim([0 100])
set(gca,'xtick',[],'xticklabel',[],'ytick',[],'yticklabel',[])
subplot(313)
[h,N]=hist(Tpeak_erstf(nn)./TD5(nn)*100,50);
Bh=bar(N,h,'facecolor',[143,188,143]/255);xlim([0 100])
set(gca,'xtick',[],'xticklabel',[],'ytick',[],'yticklabel',[])
set(gcf, 'renderer', 'painters')
print('-dpdf','-fillpage','./NEWFIGS/Fig2b.pdf');


%% part II of Fig 2 of paper with Tm
figure
subplot(311)
hist(Tpeak_stf(nn)./TM(nn)*100,50);xlim([0 100])
set(gca,'fontsize',14);xlabel('$T_M / T_D$ (\%)','Interpreter','Latex');grid on
set(gca,'TickLabelInterpreter', 'latex');
title('a)','position',[10 110],'fontsize',20);
subplot(312)
hist(Tpeak_er(nn)./TM(nn)*100,50);xlim([0 100])
set(gca,'fontsize',14);xlabel('$T_R / T_D$ (\%)','Interpreter','Latex');grid on
set(gca,'TickLabelInterpreter', 'latex');
title('b)','position',[10 110],'fontsize',20);
subplot(313)
hist(Tpeak_erstf(nn)./TM(nn)*100,50);xlim([0 100])
set(gca,'fontsize',14);xlabel('$T_S / T_D $(\%)','Interpreter','Latex');grid on
set(gca,'TickLabelInterpreter', 'latex');
title('c)','position',[10 110],'fontsize',20);
set(gcf, 'renderer', 'painters')
print('-dpdf','-fillpage','./NEWFIGS/FigS6b.pdf');


%% Fig 3 of the paper
[nn1,~]=find(mw(nn)~=0&depth(nn)<=70);  
[nn2,~]=find(mw(nn)~=0&depth(nn)>=35); 
AR_MO_DP=fitlm(log10(depth(nn(nn2))),log10(AR(nn2).*mo(nn(nn2)).^(-AR_MO.Coefficients.Estimate(2))),'linear');

figure
subplot(311)
s=scatter(log10(mo(nn)), log10(AR),'filled');grid on;hold on;
    s.MarkerFaceAlpha=0.4;
    s.MarkerEdgeColor=0.5*[1 1 1];
    s.MarkerFaceColor=[173,216,230]/255;
    X=[log10(mo(nn)), ones(size(nn))];% log10(AM),1)
    y=log10(AR);
    b=AR_MO.Coefficients.Estimate ;
    b=X\y;ycalc=X*b;
    R2=1-sum((y - ycalc).^2)/sum((y - mean(y)).^2);
    plot(log10(mo(nn)),ycalc,'b','linewidth',2,'color',[128,0,128]/255);
    X=[log10(mo(nn)), ones(size(nn))];% log10(AM),1)
    y=log10(AR);
    b3=X\y;b3(1)=-0.33;
    ycalc=X*b3;
    Y = -0.25*log10(mo(nn));
    res=mean(log10(AR)-Y);
    Y = Y + res;
    hold on;plot(log10(mo(nn)),Y,'r','linewidth',2);
    xlim([18 23])
    xlabel('$\log_{10} M_0$ (Nm)','Interpreter','Latex');
    ylabel('$\log_{10}(A_R/E_R)$','Interpreter','Latex')
    text(21,-0.25,['slope ' num2str(b(1)) ' (R-squared ' num2str(AR_MO.Rsquared.Ordinary) ')'],...
        'Interpreter','Latex','fontsize',12,'color',[128,0,128]/255)
    title('a)','position',[18 0.6],'fontsize',20);box on
    set(gca,'TickLabelInterpreter', 'latex','fontsize',14);
subplot(312)
s=scatter(log10(depth(nn)), log10(AR.*mo(nn).^(-AR_MO.Coefficients.Estimate(2))),'filled');grid on;hold on;
    s.MarkerFaceAlpha=0.4;
    s.MarkerEdgeColor=0.5*[1 1 1];
    s.MarkerFaceColor=[173,216,230]/255;
    X=[ ones(size(nn2)) ,log10(depth(nn(nn2)))];% log10(AM),1)
    b=AR_MO_DP.Coefficients.Estimate;
    ycalc=X*b;
    plot(log10(depth(nn(nn2))),ycalc,'b','linewidth',2,'color',[128,0,128]/255);
    title('b)','position',[0.5 4.6],'fontsize',20); 
    xlabel('$\log_{10} H$ (km)','Interpreter','Latex');
    ylabel('$\log_{10}( A_R/E_R)_{mo}$','Interpreter','Latex')
    text(2,3.25,['slope ' num2str(b(2)) ' (R-squared ' num2str(AR_MO_DP.Rsquared.Ordinary) ')'],...
    'Interpreter','Latex','fontsize',12,'color',[128,0,128]/255);box on
    set(gca,'TickLabelInterpreter', 'latex','fontsize',14);
subplot(313)
s=scatter(Tpeak_erstf(nn)./TD5(nn)*100,log10(AS),'filled');grid on;hold on;
    s.MarkerFaceAlpha=0.4;
    s.MarkerEdgeColor=0.5*[1 1 1];
    s.MarkerFaceColor=[255,140,0]/255;
    xlabel('$T_S/ T$ (\%)','Interpreter','Latex');
    ylabel('$\log_{10}(A_S/E_S)$','Interpreter','Latex')
    title('c)','position',[0 1.3],'fontsize',20);
    set(gca,'TickLabelInterpreter', 'latex','fontsize',14);box on
print('-dpdf','-fillpage','./NEWFIGS/Fig3.pdf')

% show examples of scaled energy functions, 3 types

ik=find(AS>2.5&AS<3.9&Tpeak_erstf(nn)./TD5(nn)>.2&Tpeak_erstf(nn)./TD5(nn)<.23&depth(nn)>100);
ik2=find(AS>1.2&AS<1.7&Tpeak_erstf(nn)./TD5(nn)>.47&Tpeak_erstf(nn)./TD5(nn)<.52&depth(nn)>100);
ik3=find(Tpeak_erstf(nn)./TD5(nn)>.9&depth(nn)>300);
figure
subplot(311)
plot(tt2,erstf_sd(nn(ik),:),'r','linewidth',2);xlim([0 1])
set(gca,'Xlabel',[],'XTick',[],'Ylabel',[],'YTick',[])
subplot(312)
plot(tt2,erstf_sd(nn(ik2),:),'r','linewidth',2);xlim([0 1])
set(gca,'Xlabel',[],'XTick',[],'Ylabel',[],'YTick',[])
subplot(313)
plot(tt2,erstf_sd(nn(ik3),:),'r','linewidth',2);xlim([0 1])
set(gcf, 'renderer', 'painters');set(gca,'Xlabel',[],'XTick',[],'Ylabel',[],'YTick',[])
print('-dpdf','-fillpage','./NEWFIGS/Fig3b.pdf')

% left bottom width height
pos1=[0.07,0.25,0.3,0.65];
pos2=[0.38,0.25,0.3,0.65];
pos3=[0.69,0.25,0.3,0.65];
pos4=[0.07,0.1,0.3,0.15];
pos5=[0.38,0.1,0.3,0.15];
pos6=[0.69,0.1,0.3,0.15];

for i=1:2
    figure(100+i-1)
ax1(i)=subplot('Position',pos1,'XTickLabel',[]);
set(ax1(i),'TickLabelInterpreter', 'latex');
ax2(i)=subplot('Position',pos2,'YTick',[],'YTickLabel',[],'XTickLabel',[]);
set(ax2(i),'TickLabelInterpreter', 'latex');
ax3(i)=subplot('Position',pos3,'YTick',[],'YTickLabel',[],'XTickLabel',[]);
set(ax3(i),'TickLabelInterpreter', 'latex');
ax4(i)=subplot('Position',pos4);
set(ax4(i),'TickLabelInterpreter', 'latex');
ax5(i)=subplot('Position',pos5,'YTick',[],'YTickLabel',[]);
set(ax5(i),'TickLabelInterpreter', 'latex');
ax6(i)=subplot('Position',pos6,'YTick',[],'YTickLabel',[]);
set(ax6(i),'TickLabelInterpreter', 'latex');
end

%% Figure 4
figure(100)
subplot(ax1(1))
pcolor(tt2,1:length(aa1),stf_sd(nn(aa),:));shading flat;hold on
plot(Tpeak_stf(nn(aa))./TM(nn(aa)),1:length(aa1),'wo','markerfacecolor','w','Markersize',1);
ylabel('Quake index');set(gca,'fontsize',12);xlim([-0.05 1.1]);caxis([0 2])
title('a)','position',[0 length(aa)+10],'fontsize',14);
subplot(ax2(1))
pcolor(tt2,1:length(aa1),er_sd(nn(aa),:));shading flat;hold on
plot(Tpeak_er(nn(aa))./TM(nn(aa)),1:length(aa1),'wo','markerfacecolor','w','Markersize',1);
caxis([0 2]);xlim([-0.05 1.1]);colormap(parula);set(gca,'YTickLabel',[]);
title('b)','position',[0 length(aa)+10],'fontsize',14);
subplot(ax3(1))
pcolor(tt2,1:length(aa),erstf_sd(nn(aa),:));shading flat;hold on
plot(Tpeak_erstf(nn(aa))./TM(nn(aa)),1:length(aa1),'wo','markerfacecolor','w','Markersize',1);
caxis([0.5 1.5]);set(gca,'YTickLabel',[]);xlim([-0.05 1.1]);
title('c)','position',[0 length(aa)+10],'fontsize',14);

%% Figure S5
figure(101)
subplot(ax1(2))
pcolor(tt2,1:length(aa1),stf_sm(nn(aa),:));shading flat;hold on
plot(Tpeak_stf(nn(aa))./TM(nn(aa)),1:length(aa1),'wo','markerfacecolor','w','Markersize',1);
ylabel('Quake index');set(gca,'fontsize',12);xlim([-0.05 1.1]);caxis([0 2])
title('a)','position',[0 length(aa)+10],'fontsize',14);
subplot(ax2(2))
pcolor(tt2,1:length(aa1),er_sm(nn(aa),:));shading flat;hold on
plot(Tpeak_er(nn(aa))./TM(nn(aa)),1:length(aa1),'wo','markerfacecolor','w','Markersize',1);
caxis([0 2]);xlim([-0.05 1.1]);colormap(parula);set(gca,'YTickLabel',[]);
title('b)','position',[0 length(aa)+10],'fontsize',14);
subplot(ax3(2))
pcolor(tt2,1:length(aa),erstf_sm(nn(aa),:));shading flat;hold on
plot(Tpeak_erstf(nn(aa))./TM(nn(aa)),1:length(aa1),'wo','markerfacecolor','w','Markersize',1);
caxis([0.5 1.5]);set(gca,'YTickLabel',[]);xlim([-0.05 1.1])
title('c)','position',[0 length(aa)+10],'fontsize',14);

TT=logspace(log10(0.025),log10(0.3),100);
G1 = [ones(length(TT),1) log10(TT)'];

%% Fig Supp 2
figure(5002)
loglog(FF,xspec2(nn,:),'k','color',0.4*[1 1 1],'linewidth',0.5);grid on;hold on

figure(5003)
loglog(FF,xspec3(nn,:),'k','color',0.4*[1 1 1],'linewidth',0.5);grid on;hold on



for ih=1:length(HH2)-1
    ik=find(depth(nn)>HH2(ih)&depth(nn)<HH2(ih+1));
    disp(length(ik))
    if isempty(ik);continue;end
    
    % make better median
    med_erstf=zeros(size(tt2));med_stf=med_erstf;med_er=med_stf;med_es=med_er;
    med_erstfm=zeros(size(tt2));med_stfm=med_erstf;med_erm=med_stf;med_esm=med_er;
    for i=1:length(tt2)
        ik1=find(isnan(stf_sd(nn(ik),i))==0&isinf(stf_sd(nn(ik),i))==0&stf_sd(nn(ik),i)~=0);
        med_stf(i) = median(stf_sd(nn(ik(ik1)),i));
        ik1=find(isnan(erstf_sd(nn(ik),i))==0&isinf(er_sd(nn(ik),i))==0&er_sd(nn(ik),i)~=0);
        med_er(i) = median(er_sd(nn(ik(ik1)),i));
        ik1=find(isnan(erstf_sd(nn(ik),i))==0&isinf(erstf_sd(nn(ik),i))==0&erstf_sd(nn(ik),i)~=0);
        med_erstf(i) = median(erstf_sd(nn(ik(ik1)),i));
        
        
        ik1=find(isnan(stf_sm(nn(ik),i))==0&isinf(stf_sm(nn(ik),i))==0&stf_sm(nn(ik),i)~=0);
        med_stfm(i) = median(stf_sm(nn(ik(ik1)),i));
        ik1=find(isnan(erstf_sm(nn(ik),i))==0&isinf(er_sm(nn(ik),i))==0&er_sm(nn(ik),i)~=0);
        med_erm(i) = median(er_sm(nn(ik(ik1)),i));
        ik1=find(isnan(erstf_sm(nn(ik),i))==0&isinf(erstf_sm(nn(ik),i))==0&erstf_sm(nn(ik),i)~=0);
        med_erstfm(i) = median(erstf_sm(nn(ik(ik1)),i));
    end
    
    disp([skewness(med_stf) kurtosis(med_stf)])
    % make median of spectra
    for i=1:length(FF)
        medspec(i)=10.^(mean(log10(xspec2(nn(ik),i))));
        medspec2(i)=10.^(mean(log10(xspec3(nn(ik),i))));
    end
    
    figure(5002)
    loglog(FF,medspec,'k','color',CC3(ih,1:3),'linewidth',2);grid on
    
    figure(5003)
    loglog(FF,medspec2,'k','color',CC3(ih,1:3),'linewidth',2);grid on
    
    
    
    [~,it1]=max(med_stf);
    [~,it2]=max(med_er);
    [~,it3]=max(med_erstf);
    it=find(tt2>0&tt2<=tt2(min([it1 it2 it3])));
  
% 
%     
%     resmin1=9E9; 
%     II=find(tt2>=0&tt2<=1);   
%     for ig=1:length(agauss)
%     for isk=1:length(ske)
%         res1 = norm([med_erstf(II)-es{ig,isk}(II) ...
%             med_er(II)-sacc{ig,isk}(II) med_stf(II)-Y{ig,isk}(II)],2);
%         if res1<=resmin1
%             resmin1=res1;
%             agauss_best(ih+length(F))=agauss(ig);
%             ske_best(ih+length(F))=ske(isk);
%         end
%     end
%     end
%         disp([agauss_best(ih+length(F)) ske_best(ih+length(F))])
%         
%         
%         crap=skewedgaussian(tt3*agauss_best(ih+length(F)), ske_best(ih+length(F)));
%         ik=find(crap>=0.01*max(crap));
%         Y1=zeros(size(tt2));
%         Y1(1:min(length(ik),length(tt2)))=crap(ik(1:min(length(ik),length(tt2))));
%         Y1=Y1/sum(Y1*dt2);
%         sacc1=(diff(Y1)/dt2).^2;sacc1(end+1)=0;
%         es1=cumtrapz(sacc1.*dt2)./cumtrapz(Y1*dt2);es1=es1/es1(end);
%         sacc1=sacc1/sum(sacc1(tt2>=0&tt2<=1)*dt2);
% 
%         figure(1)
%         subplot(311)
%         plot(tt2,med_stf(:),tt2,Y1,'linewidth',2);grid on;xlim([0 1.1])
%         title(['mw ' num2str(mw(i)) 'skewness ' num2str(ske_best(i))])
%         subplot(312)
%         plot(tt2,med_er(:),tt2,sacc1,'linewidth',2);grid on;xlim([0 1.1])
%         subplot(313)
%         plot(tt2,med_erstf(:),tt2,es1,'linewidth',2);grid on;xlim([0 1.1]);legend('data','model')
% pause

    figure(100)
    subplot(ax4(1))
    plot(tt2*100,med_stf,'linewidth',2,'color',CC3(ih,1:3));grid on;hold on 
    xlabel('$t/T_D (\%)$','Interpreter','latex');xlim([-5 110]) ;ylim([0 2])
    set(ax4(1),'TickLabelInterpreter', 'latex','fontsize',12);
    
    subplot(ax5(1))
    plot(tt2*100,med_er,'linewidth',2,'color',CC3(ih,1:3));grid on;hold on 
    xlabel('$t/T_D (\%)$','Interpreter','latex');xlim([-5 110])  ;ylim([0 2])
    set(ax5(1),'TickLabelInterpreter', 'latex','fontsize',12,'YTickLabel',[]);
    
    subplot(ax6(1))
    plot(tt2*100,med_erstf,'linewidth',2,'color',CC3(ih,1:3));grid on;hold on
    xlabel('$t/T_D (\%)$','Interpreter','latex');xlim([-5 110])  ;ylim([0 2])
    set(ax6(1),'TickLabelInterpreter', 'latex','fontsize',12,'YTickLabel',[]);

    figure(101)
    subplot(ax4(2))
    plot(tt2*100,med_stfm,'linewidth',2,'color',CC3(ih,1:3));grid on;hold on 
    xlabel('$t/T_D (\%)$','Interpreter','latex');xlim([-5 110])  ;ylim([0 2.5])
    set(ax4(2),'TickLabelInterpreter', 'latex','fontsize',12);
    
    subplot(ax5(2))
    plot(tt2*100,med_erm,'linewidth',2,'color',CC3(ih,1:3));grid on;hold on 
    xlabel('$t/T_D (\%)$','Interpreter','latex');xlim([-5 110])  ;ylim([0 2.5])
    set(ax5(2),'TickLabelInterpreter', 'latex','fontsize',12,'YTickLabel',[]);
    
    
    subplot(ax6(2))
    plot(tt2*100,med_erstfm,'linewidth',2,'color',CC3(ih,1:3));grid on;hold on
    xlabel('$t/T_D (\%)$','Interpreter','latex');xlim([-5 110]);ylim([0 2.5])
    set(ax6(2),'TickLabelInterpreter', 'latex','fontsize',12,'YTickLabel',[]);
    
    
end
figure(100)
set(gcf,'PaperUnits','inches','PaperPosition',[0.1 0.1 8 11]);
print('-dpng','-r300','./NEWFIGS/Fig4.png');
figure(101)
set(gcf,'PaperUnits','inches','PaperPosition',[0.1 0.1 8 11]);
print('-dpng','-r300','./NEWFIGS/FigS7.png');

figure(5002)
set(gcf, 'renderer', 'painters')
set(gca,'fontsize',14);xlabel('Normalized frequency F" ');
ylabel('$M_0$ (Nm)','Interpreter','latex')
title(['FAS of normalized STFs ' ])
set(gcf, 'renderer', 'painters')
print('-dpdf',['./NEWFIGS/FigS8.pdf']);

figure(5003)
set(gcf, 'renderer', 'painters')
set(gca,'fontsize',14);xlabel('Normalized frequency F" ');
ylabel('$M_0$ (Nm)','Interpreter','latex')
title(['FAS of normalized STFs ' ])
set(gcf, 'renderer', 'painters')
print('-dpdf',['./NEWFIGS/Fig_resp_spec.pdf']);
