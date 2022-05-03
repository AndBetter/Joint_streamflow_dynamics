clear
pkg load statistics

x_min=0.5
x_max=300
y_min=0.5
y_max=300
step=1

[X,Y]=meshgrid((x_min:step:x_max), (y_min:step:y_max));
x=X(1,:)';
y=Y(:,1);
x_ve=X(:);
y_ve=Y(:);


% f_1: hougaard distribution - t e f: medie ; r>0; when r=1 --> corr=0; when 0<r<1 --> corr<0 (ma fa cose strane!!)

r=2
t=60
f=30

F_1= exp( -  (  (X./t).^r  +  (Y./f).^r  ).^(1./r)  );

f_1= F_1 .* ((X.*Y).^(r.-1))./(t.*f).^r  .*  (  (X./t).^r  +  (Y./f).^r  ).^(1./r - 2)  .*  ( r-1 + (  (X./t).^r  +  (Y./f).^r  ).^(1./r) )  ;

corr_teo_f_1=  1/gamma(2/r) .* (gamma(1/r)).^2 .*(gammainc( xx=1/r , a= 1, type='upper') ).^2 ./ (r * gammainc( xx=2/r , a= 1, type='upper') )  -1     %(gamma(1/r))^2/(r*gamma(2/r)) - 1

corr_teo_f_1= (gamma(1/r))^2/(r*gamma(2/r)) - 1

r_series=(0:0.01:10);
corr_vs_r=gamma(1./r_series).^2 ./(r_series.*gamma(2./r_series)) - 1;
%plot(r_series,corr_vs_r)
%f_1(f_1>0.001)=0
%surf((f_1))


% f_2: downton distribution 1/m1 e 1/m2 : medie

rho=0.7
m1=0.1
m2=0.05
f_2= m1.*m2./(1-rho) .*  exp (  - (m1.*X + m2.*Y) ./ (1-rho)  )  .* besseli(0,  2*sqrt(rho.*m1.*m2.*X.*Y)./(1-rho))  ;

%f_2(f_2>0.001)=0
surf(f_2)



##%% sample from distribution version 1
##f_ve=f_2(:);
##data=[];
##temp1=[];
##temp2=[];
##for i=1: length(x_ve)
##  temp1=repmat(x_ve(i), round(1000*f_ve(i)));
##  temp2=repmat(y_ve(i), round(1000*f_ve(i)));
##  data=[data;[temp1(:),temp2(:)]];
##end



%%sample from distribution version 2

number_samples=10000          # <----- how many samples to be drawn              
sampling_distribution= f_2;   # <----- distribution from where to sample

xy=[]
sum_distribution = nansum(sampling_distribution(:));
for i=1:length(x)
xy=[xy; [x repmat(y(i), length(x),1) ]];
end
%%generate random row number of (x,y) pair
rows = randsample(length(xy),number_samples,true, sampling_distribution(:) );  % select the distribution to use for sampling (last argument)
xy_sample = xy(rows,:);