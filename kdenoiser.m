function[out_put_signal]=kdenoiser(y,fs)

% Wavelet subband dependent thresholding for denoising of phonocardiographic signals, 
% Kritika Agrawal, Abhinash Kumar Jha, Shealini Sharma, Ayush Kumar, Vijay S Chourasia
% Signal Processing: Algorithms, Architectures, Arrangements, and Applications (SPA), 2013
% Pages:158-162, Publisher-IEEE
%% parameters
% Input   :   y                --> Noisy Signal
%             fs               --> Sampling Frequency
% Output  :   out_put_signal   --> Denoised Signal.
%% code goes here
T=1/fs;
L=size(y);
NFFT=2^nextpow2(L(1,1));
f = fs/2*linspace(0,1,NFFT/2+1);

fpass=10;
fstop=300;
apass=1;
astop=60;
fs=8000;

%----passing low pass filter to noisy signal-----


h=fdesign.lowpass('fp,fst,ap,ast',fpass,fstop,apass,astop,fs);

Hd=design(h, 'equiripple',...
    'MinOrder', 'any',...
    'StopbandShape', 'flat');
low_filter=filter(Hd,y);

z=fft(low_filter,NFFT)/L(1,1);

%----ploting amplitude responce of noisy signal after passing low pass


%----applying dwt on filtered noisy signal-----

[c p1]=wavedec(low_filter,2,'coif5');


k=0.9;
m=0.9;

%------finding threshold value for ca5---------

ca5=appcoef(c,p1,'coif5',2);

%------finding threshold value for cd1---------
mean=0;
cd1=detcoef(c,p1,1);
for i=1:length(cd1)
    mean=mean+cd1(i,1);
end
mean=mean/length(cd1);
std_deviation=0;
for i=1:length(cd1)
    std_deviation=std_deviation+(cd1(i,1)-mean).^2;
end
std_deviation=sqrt(std_deviation/(length(cd1)-1));
thr2=(k*m*std_deviation*(sqrt(2*log2(length(cd1)))));



for i=1:length(cd1)
       if abs(cd1(i,1))>=thr2
          cd1(i,1)=cd1(i,1);%sign(cd1(i,1))*(abs(cd1(i,1))-thr3);
       else
          cd1(i,1)=0;
       end
 end

l1=size(cd1);
NFFT=2^nextpow2(l1(1,1));
k1=fft(cd1,NFFT)/l1(1,1);
f1 = fs/2*linspace(0,1,NFFT/2+1);
cd2=detcoef(c,p1,2);

%----finding thr value for cd2------
mean=0;
for i=1:length(cd2)
    mean=mean+cd2(i,1);
end
mean=mean/length(cd2);
std_deviation=0;
for i=1:length(cd2)
    std_deviation=std_deviation+(cd2(i,1)-mean).^2;
end
std_deviation=sqrt(std_deviation/(length(cd2)-1));
thr3=(k*m*std_deviation*(sqrt(2*log2(length(cd2)))));

for i=1:length(cd2)
       if abs(cd2(i,1))>=thr3
          cd2(i,1)=cd2(i,1);
       else
          cd2(i,1)=0;
       end
 end




l2=size(cd2);
NFFT=2^nextpow2(l2(1,1));
k2=fft(cd2,NFFT)/l2(1,1);
f2 = fs/2*linspace(0,1,NFFT/2+1);
k=0.2;
m=0.2;
 
mean=0;
for i=1:length(ca5)
    mean=mean+ca5(i,1);
end
mean=mean/length(ca5);
std_deviation=0;
for i=1:length(ca5)
    std_deviation=std_deviation+(ca5(i,1)-mean).^2;
end
std_deviation=sqrt(std_deviation/(length(ca5)-1));
thr3=(k*m*std_deviation*(sqrt(2*log2(length(ca5)))));

for i=1:length(ca5)
       if abs(ca5(i,1))>=thr3
          ca5(i,1)=sign(ca5(i,1))*(abs(ca5(i,1))-thr3);
       else
          ca5(i,1)=0;
       end
end
 
 c=[ca5;cd1;cd2];%;cd3;cd4;cd5];
 
 out_put_signal=waverec(c,p1,'coif5');

end


