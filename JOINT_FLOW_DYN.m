% joint stramflow dynamics --- MIGLIORE VERSIONE DISPONIBILE !!!!! 

clear 
Delta_t=0.1        % (days) - temporal discretization  - può essere necessario ridurre la discretizzatzione per evitare problemi di eveni che non vengono contati
simlength=10000;   % lunghezza della simulazione (l'unità di misura dipende da quella dei parametri sotto - giorni o decimi di giorno)

% assign the parameters of the model (solution 2)

L1 =0.2   %lambda1 disjoint
L2 =0.2    %lambda2 disjoint
L12=0.02    %lambda12 joint

a1=15      %alpha1 disjoint
a2=15      %alpha2 disjoint
a1_12=30   %alpha1_12 joint
a2_12=30   %alpha2_12 joint
r=0.8      %correlation between joint jumps

k1=0.2     %recessions
k2=0.2

% genera una sequenza di eventi con frequenca L1+L2+L12 (hp: eventi joint e disjoint indipendenti)

L=L1+L2+L12;


tot=1 * simlength;         % numerosità delle realizzazzioni stocastiche degli eventi (da aumentare nel caso dasse errori)
t=zeros(tot,1);            % tempo

% crea interarrivi esponenzialmente distribuiti
dt=(rand(tot,1));
dt=-log(dt)/L;    
    

%dt=round(dt);  % discretizza--- meglio discretizzare dt o t ??? bohhhh


% time at wich all  events occur
t(1)=0;
for i=2:tot
    t(i,1)=t(i-1)+dt(i-1);   
end


t=round(t);  % discretizza--- meglio discretizzare dt o t ??? bohhhh


%crea la sequenza di tempi in cui un evento avviene (rimuove duplicati)
eventi_rimossi=0;
for i=1:(tot-1)
    if t(i)==t(i+1);
         eventi_rimossi(i)=1;  % conta gli eveni che vengono rimossi perchè concomitanti
        t(i)=0;
    else
         eventi_rimossi(i)=0;
    end
end

eventi_rimossi_percent=sum(eventi_rimossi)/length(eventi_rimossi);

t(t==0)= [];
 
% probabilità di perdere un evento per via di dt<0.5 giorni
P=1-exp(-L*0.5);
 
P1=(1-exp(-L1*0.5));
P2=(1-exp(-L2*0.5));
P12=(1-exp(-L12*0.5));

%  P1_=P*L1/L;
%  P2_=P*L2/L;
%  P12_=P*L12/L;


% altezze di pioggia disjoint esponenzialmente distribuite
h1=(rand(tot,1));
h1=-log(h1)*a1;

h2=(rand(tot,1));
h2=-log(h2)*a2;


##%% generate joint pdf of jumps - versione originale modello articolo mio
##        
##a=a2_12/a1_12 * r;  %!!!!
##rho=a1_12/a2_12;
##
##
##h=(rand(tot,1));
##hh=(rand(tot,1));
##
##h1_12=-log(h)*a1_12;
##h2_12temp=-log(hh)*a1_12;
##
##rB=(rand(tot,1));
##
##b=0;
##for i=1:length(rB)
##    
##  if rB(i)<=(a*rho);
##      rB(i)=0;
##  else
##      rB(i)=1;
##  end
##end
##
##b=rB;
##
##h2_12=a*h1_12+b.*h2_12temp;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% generate joint pdf of jumps - versione alternativa: f_2: downton distribution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
pkg load statistics

x_min=0.1
x_max=300
y_min=0.1
y_max=300
step=1

[X,Y]=meshgrid((x_min:step:x_max), (y_min:step:y_max));
x=X(1,:)';
y=Y(:,1);
x_ve=X(:);
y_ve=Y(:);


rho=r
m1=1/a1_12
m2=1/a2_12
f_2= m1.*m2./(1-rho) .*  exp (  - (m1.*X + m2.*Y) ./ (1-rho)  )  .* besseli(0,  2*sqrt(rho.*m1.*m2.*X.*Y)./(1-rho))  ;

%%sample from distribution f_2

number_samples=tot            # <----- how many samples to be drawn              
sampling_distribution= f_2;   # <----- distribution from where to sample

xy=[];
sum_distribution = nansum(sampling_distribution(:));
for i=1:length(x)
xy=[xy; [x repmat(y(i), length(x),1) ]];
end
%%generate random row number of (x,y) pair
rows = randsample(length(xy),number_samples,true, sampling_distribution(:) );  % select the distribution to use for sampling (last argument)
xy_sample = xy(rows,:);

h1_12=xy_sample(:,1);
h2_12=xy_sample(:,2);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%creates sequences of joint and disjoint effective rainfall events
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c1=zeros(simlength,1);
c2=zeros(simlength,1);

t1=zeros(length(t),1);
t2=zeros(length(t),1);
t12=zeros(length(t),1);


for i=1:length(t)
    casuale=rand;
    
    if casuale<L1/L;
        t1(i)=t(i);
        
    elseif casuale>=L1/L && casuale<(L1+L2)/L ;
        t2(i)=t(i);
        
    elseif casuale>=(L1+L2)/L;
        t12(i)=t(i);
    end
end
        
t1(t1==0)= [];      
t2(t2==0)= [];
t12(t12==0)= [];

c1(t1)=1;
c1(t12)=2;
c1=c1(1:simlength);

c2(t2)=1;
c2(t12)=2;
c2=c2(1:simlength);

c1(1,2)=0;
c2(1,2)=0;


% h1_corrected=h1*(1+P*L1/L);   % tenta di compensare il fatto che discretizzando si perdono eventi
% h2_corrected=h2*(1+P*L2/L);
% h1_12_corrected=h1_12*(1+P*L12/L);
% h2_12_corrected=h2_12*(1+P*L12/L);



% assigns the sequence of joint and disjoint events and try to correct their intensity to compensate for the events that get lost due to time discretization (especially in the dailiy discretization)
for i=1:simlength
    temp=rand;  
    if c1(i,1)==1
        c1(i,2)=h1(i) *   (a1*L1)  /  (  sum(h1(c1(:,1)==1)) /  length (c1)  ); % *   (a1*L1)  /  (  sum(h1(c1(:,1)==1)) /  length (c1)  );     %+ (rand<P1)*h1(1+round(rand*(length(h1)-1)));   %P*L1/L     % cerca di compensare il fatto che discretizzando si perdono eventi
        elseif  c1(i,1)==2;
        c1(i,2)=h1_12(i)  *   (a1_12*L12)  /  (  sum(h1_12(c1(:,1)==2)) /  length (c1)  );  %*   (a1_12*L12)  /  (  sum(h1_12(c1(:,1)==2)) /  length (c1)  );    %+ (temp<P12)*h1_12(1+round(temp*(length(h1_12)-1)));   %P*L12/L in alternativa per la probabilità di salti sovrapposti
    end
    
     if c2(i,1)==1
        c2(i,2)=h2(i) *   (a2*L2)  /  (  sum(h2(c2(:,1)==1)) /  length (c2)  ); % *   (a2*L2)  /  (  sum(h2(c2(:,1)==1)) /  length (c2)  );   %+(rand<P2)*h2(1+round(rand*(length(h2)-1)));     %P*L2/L
        elseif c2(i,1)==2;
        c2(i,2)=h2_12(i)  *   (a2_12*L12)  /  (  sum(h2_12(c2(:,1)==2)) /  length (c2)  );  % *   (a2_12*L12)  /  (  sum(h2_12(c2(:,1)==2)) /  length (c2)  );    %+ (temp<P12)*h2_12(1+round(temp*(length(h2_12)-1)));   %P*L12/L
     end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% creates streamflow timeseries 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
q1=[0;0];
q2=[0;0];

q01=0;
q02=0;
qbase1=0;
qbase2=0;
aux1=0;
aux2=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Forma originaria recessioni - conserva i volumi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=2:simlength      % recessione sempre in corso
  
    if c1(i,1)~=0    
    q01=c1(i,2)*k1+qbase1;  
    q1(i)=(q01/k1)*(1-exp(-k1));                  % forse metodo che non mantiene i volumi è meglio
    qbase1= ( q1(i)/k1)*(1-exp(-k1));
    %qbase1= q01*exp(-1*k1);      %qbase1=(q01/k1)*(exp(-0.5*k1)-exp(-1.5*k1));
    aux1=0;
    else    
    q1(i)=(q01/k1)*(exp(-(1+aux1)*k1)-exp(-(2+aux1)*k1));        %   q1(i-1)+q1(i-1)*(exp(-k1)-1);      % q1(i-1)-q1(i-1)*k1*exp(-k1*0.5);  %
    qbase1=( q1(i)/k1)*(1-exp(-k1));
    %qbase1=q01*exp(-(2+aux1)*k1);    %   qbase1=(q01/k1)*(exp(-(0.5+aux1+1)*k1)-exp(-(1.5+aux1+1)*k1));
    aux1=aux1+1;
    end
    
end

for i=2:simlength      % recessione sempre in corso
  
    if c2(i,1)~=0    
    q02=c2(i,2)*k2+qbase2;  
    q2(i)=(q02/k2)*(1-exp(-k2));                  % forse metodo che non mantiene i volumi è meglio
    qbase2= ( q2(i)/k2)*(1-exp(-k2));
    %qbase2= q02*exp(-1*k2);    % qbase2=(q02/k2)*(exp(-0.5*k2)-exp(-1.5*k2));
    aux2=0;
    else    
    q2(i)=(q02/k2)*(exp(-(1+aux2)*k2)-exp(-(2+aux2)*k2));        %   q1(i-1)+q1(i-1)*(exp(-k1)-1);      % q1(i-1)-q1(i-1)*k1*exp(-k1*0.5);  %
    qbase2=q2(i)*exp(-1*k2);
    qbase2=( q2(i)/k2)*(1-exp(-k2));
    %qbase2=q02*exp(-(2+aux2)*k2);      %  qbase2=(q02/k2)*(exp(-(0.5+aux2+1)*k2)-exp(-(1.5+aux2+1)*k2));
    aux2=aux2+1;
    end
    
end


%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Forma alternativa recessioni  Prova --- sembra andare meglio dell'alta!!!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% for i=2:simlength      
%   
%     if c1(i,1)~=0    
%     q01=c1(i,2)*k1+qbase1;  
%     q1(i)=q01*exp(-k1*0.5);                  % forse metodo che non mantiene i volumi è meglio
%     qbase1=q01*exp(-k1*1);
%     aux1=0;
%     else    
%     q1(i)=q01*exp(-k1*(0.5+1+aux1));        %   q1(i-1)+q1(i-1)*(exp(-k1)-1);      % q1(i-1)-q1(i-1)*k1*exp(-k1*0.5);  %
%     qbase1=q01*exp(-k1*(0.5+1.5+aux1));
%     aux1=aux1+1;
%     end
%     
% end
% 
% 
% for i=2:simlength      % recessione sempre in corso
%   
%     
%     if c2(i,1)~=0    
%     q02=c2(i,2)*k2+qbase2;  
%     q2(i)=q02*exp(-k2*0.5);                  % forse metodo che non mantiene i volumi è meglio
%     qbase2=q02*exp(-k2*1);
%     aux2=0;
%     else    
%     q2(i)=q02*exp(-k2*(0.5+1+aux2));        %   q1(i-1)+q1(i-1)*(exp(-k1)-1);      % q1(i-1)-q1(i-1)*k1*exp(-k1*0.5);  %
%     qbase2=q02*exp(-k2*(0.5+1.5+aux2));
%     aux2=aux2+1;
%     end
%     
% end



%%

corr_simulated=corr(q1,q2)
corr_model2=(L12*a1_12*a2_12)/sqrt((L1*a1^2+L12*a1_12^2)*(L2*a2^2+L12*a2_12^2))*0.5*(1+r)*2*sqrt(k1*k2)/(k1+k2)


Nash_Norm=1-sum((q2/mean(q2)-q1/mean(q1)).^2)/sum((q1/mean(q1)-mean(q1/mean(q1))).^2)




%ottiene serie temporale giornaliera aggrega la serie temporale - da usare se si sono scalati i parametri lambda, alpha e k

DT=10; % DT=10 se i parametri lambda e k sono stati ridotti di un fattore 10 (e alpha aumentato di un fattore 10)

q1_daily=zeros(simlength/DT-1,1);
q2_daily=zeros(simlength/DT-1,1);
day=zeros(simlength/DT-1,1);
for i=1:(simlength/DT-1)
    q1_daily(i)=mean(q1(DT*(i-1)+1:DT*(i-1)+1+DT));
    q2_daily(i)=mean(q2(DT*(i-1)+1:DT*(i-1)+1+DT));
    day(i)=(i-1)*DT+DT/2;
end
% 


%%CDF
bins_cdf=[0:0.001:max(max(q1), max(q2))];
cdf1= empirical_cdf(bins_cdf,q1);
cdf2= empirical_cdf(bins_cdf,q2);


%% remap q1 and q2 using cdf1 into q1_trans and q_2 trans (non exceedance probability)
q1_trans=zeros(length(q1),1);
q2_trans=zeros(length(q1),1);

for i=1:length(q1)
[val,idx1]=min(abs(q1(i)-bins_cdf));
[val,idx2]=min(abs(q2(i)-bins_cdf));

q1_trans(i)= cdf1(idx1);
q2_trans(i)= cdf2(idx2);
end

corr_simulated_exceedence=corr(q1_trans,q2_trans)


%% creates the binary timeseries of exceedance of a specific threshold discharge

q_threshold=[0:0.001:0.99]
corr_simulated_binary=zeros(q_threshold,1);
for i=1:length(q_threshold)
  q1_trans_binary=zeros(length(q1),1);
  q2_trans_binary=zeros(length(q2),1);

  q1_trans_binary(q1_trans> q_threshold(i)) = 1;
  q2_trans_binary(q2_trans> q_threshold(i)) = 1;

  [t] =crosstab(q1_trans_binary,q2_trans_binary);
  corr_simulated_binary(i)=corr(q1_trans_binary,q2_trans_binary);
end

scatter(q_threshold,corr_simulated_binary)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xaxis_min=1000
xaxis_max=1365
y_axis_max_1=quantile(q1,0.99)
y_axis_max_2=quantile(q2,0.99)

figure()
subplot(3,3,1)
plot(q1)
hold on
plot(q2)
xlim([xaxis_min,xaxis_max])
ylim([0,max(y_axis_max_1,y_axis_max_2)])
legend('q_1 (mm/day)','q_2 (mm/day)', 'FontSize',20)
xlabel('time (days)','FontSize',20) 
ylabel('discharge (mm/day)','FontSize',20) 

subplot(3,3,2)
hist(q1, [0.25:0.5:max(max(q1), max(q2))])
xlabel('q_1 (mm/day)','FontSize',20) 
ylabel('frequency (-)','FontSize',20) 

subplot(3,3,3)
hist(q2, [0.25:0.5:max(max(q1), max(q2))])
xlabel('q_2 (mm/day)','FontSize',20) 
ylabel('frequency (-)','FontSize',20) 

%subplot(3,3,5)
annotation('textbox',...
    [0.5 0.4 0.3 0.15],...
    'String',{['\lambda_1 =' num2str(L1)],['\lambda_2 =' num2str(L2)], ['\lambda_{12} =' num2str(L12)] , ['\alpha_1 =' num2str(a1)], ['\alpha_1^{12} =' num2str(a1_12)], ['\alpha_2 =' num2str(a2)], ['\alpha_2^{12} =' num2str(a2_12)], ['k_1 =' num2str(k1)], ['k_2=' num2str(k2)], ['r =' num2str(r)]},...
    'FontSize',14,...
    'FontName','Arial',...
    'LineStyle','-',...
    'EdgeColor',[1 1 0],...
    'LineWidth',2,...
    'BackgroundColor',[0.9  0.9 0.9],...
    'Color',[0.84 0.16 0]);

subplot(3,3,6)
plot(1-cdf1,bins_cdf)
hold on
plot(1-cdf2,bins_cdf)
xlabel('exceedance prob (-)','FontSize',20) 
ylabel('q (mm/day)','FontSize',20) 
legend('q_1','q_2', 'FontSize',20)
set(gca, 'YScale', 'log')

subplot(3,3,4)
plotyy([0:1:length(q1)-1], q1,[0:1:length(q1)-1], (1-q1_trans))
#plot(q1)
#hold on
#plot(1-q1_trans)
#line ([0 length(q1)], [1 1], "linestyle", ":", "color", "r");
legend('q_1 (mm/day)','% time exceeded (-)', 'FontSize',20) 
xlabel('time (days)','FontSize',20) 
xlim([xaxis_min,xaxis_max])
ylim([0,y_axis_max_1])


subplot(3,3,7)
plotyy([0:1:length(q2)-1], q2,[0:1:length(q2)-1], (1-q2_trans))
#plot(q2)
#hold on
#plot(1-q2_trans)
#line ([0 length(q2)], [1 1], "linestyle", ":", "color", "r");
legend('q_2 (mm/day)','% time exceeded (-)', 'FontSize',20) 
xlabel('time (days)','FontSize',20) 
xlim([xaxis_min,xaxis_max])
ylim([0,y_axis_max_2])

subplot(3,3,8)
scatter(q1,q2, '.')
xlabel('q_1 (mm/day)','FontSize',20) 
ylabel('q_2 (mm/day)','FontSize',20) 
legend(sprintf('Correlation= %f ', corr(q1,q2)),'% time exceeded (-)', 'FontSize',20)  

subplot(3,3,9)
scatter(1-q1_trans,1-q2_trans, '.')
xlabel('exceedance prob q_1 (-)','FontSize',20) 
ylabel('exceedance prob q_2 (-)','FontSize',20) 
legend(sprintf('Correlation= %f ', corr(1-q1_trans,1-q2_trans)),'% time exceeded (-)', 'FontSize',20)  


##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##%% comparison with analytical streamflow PDF
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##
binwidth=1
figure()
shape=(L1+L12)/k1
scale=(L1*a1 + L12*a1_12)/(L1+L12) *k1
media_teorica= L1*a1 + L12*a1_12
media_simulata= mean(q1)
[value bins]=hist(q1, [binwidth/2:binwidth:max(max(q1), max(q2))])
hist(q1, [binwidth/2:binwidth:max(max(q1), max(q2))])
hold on
plot(bins,simlength*binwidth*gampdf(bins,shape,scale),"r","linewidth", 2)
xlim([0,quantile(q1,0.999)])
xlabel('q (mm/day)','FontSize',20) 
legend('Frequency (-)','PDF Teorica: p(q)=c q^{\lambda/k-1}exp(-q k/\alpha)', 'FontSize',20)
