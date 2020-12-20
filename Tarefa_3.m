%% Code designed by Victor Fabre Figueiredo %%
%% Tarefa 3 da matéria Processamento de Sinais Digitais - UnB - 2020
% Realizar a filtragem de um sinal de áudio utilizando um filtro passa-baixa 
% por 3 métodos diferentes:
%   1 - Convolução no domínio do tempo;
%   2 - Overlap and add;
%   3 - Overlap and save.
% Todos os três métodos devem apresentar os mesmos resultados.
%% Read audio file and save only the first minute of it
fileName = 'Back_In_Black.mp3';
[orig_audio,Fs] = audioread(fileName);
[x,fs] = audioread(fileName, [1, 60*Fs]); %first minute of the audio
%% Convertion of stereo audio to mono audio (DO NOT CHANGE!) %%
[m, n] = size(x);

if n == 2
    y = x(:, 1) + x(:, 2);
    peakAmp = max(abs(y)); 
    y = y/peakAmp;
    peakL = max(abs(x(:, 1)));
    peakR = max(abs(x(:, 2))); 
    maxPeak = max([peakL peakR]);
    y = y*maxPeak;    
else
    y = x;
end
%% Generate or load filter
%low_pass = designfilt('lowpassfir', 'FilterOrder', 1000, 'CutoffFrequency', 500, 'SampleRate', 44100);
low_pass = load('low_pass_filter.mat');
low_pass = low_pass.low_pass;
low_pass_coeff = low_pass.Coefficients;
%low_pass_coeff = load('low_pass_coefficients.mat');
%% Do your processing here

% Transforma coluna em linha. No processo de concatenação do overlap
% and add/save, o vetor precisa ser vetor-linha e não vetor-coluna. 
y = y.';

x = y;
h = low_pass_coeff;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%1 - Convolução no domínio do tempo;
convolution=conv(x,h);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%2 - Overlap and add;
L=32;
N1=length(x);
M=length(h);

% preenche x com zeros
x=[x zeros(1,mod(-N1,L))];
N2=length(x);
h=[h zeros(1,L-1)];
% aplica fft no filtro
H=fft(h,L+M-1);
S=N2/L;
index=1:L;
X=[zeros(M-1)];

for intervalo=1:S
    % Seleciona a sequência para processar.
    xm=[x(index) zeros(1,M-1)];		
    % aplica fft
    X1=fft(xm,L+M-1);
    % realiza o produto de ffts
    Y=X1.*H;
    % cacula a fft inversa
    Y=ifft(Y);
    % Adiciona amostra em cada estágio
    Z=X((length(X)-M+2):length(X))+Y(1:M-1);	
    X=[X(1:(intervalo-1)*L) Z Y(M:M+L-1)];
    index=intervalo*L+1:(intervalo+1)*L;
end
i=1:N1+M-1;
X=X(i);
overlap_add=X;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3 - Overlap and save

x = y;
h = low_pass_coeff;

% tamanho do filtro
P = length(h);
% tamanho da entrada
N = length(x);
D = N-P+1;
x = [zeros(P-1,1)' x];
%calcula fft
H = fft(h,N);
ys = [];
for r = 0:fix(N/D)-1
    xr = x(r*D+1:r*D+N);
    Xr = fft(xr);
    Yrp = Xr .* H;
    yrp = ifft(Yrp);
    ys = [ys; yrp(P:N)];
end

figure()
subplot(3,1,1)
stem(convolution);
title('Convolução usando a função conv()')
xlabel('n');
ylabel('y(n)');

subplot(3,1,2)
stem(overlap_add);
title('Convolução usando o método Overlap Add')
xlabel('n');
ylabel('y(n)');

subplot(3,1,3)
overlap_save=ys;
stem(overlap_save);
title('Convolução usando o método Overlap Save')
xlabel('n');
ylabel('y(n)');

%% Calculo do erro
%Erro médio entre convolution e overlap_add   
diferenca = abs(overlap_add - convolution);
erro_media_conv_add= mean(diferenca); % erro médio.
display(erro_media_conv_add)

%Erro médio entre convolution e overlap_save   
overlap_save = [overlap_save zeros(1,2000)];
diferenca = abs(overlap_save - convolution);
erro_media_conv_save= mean(diferenca); % erro médio.
display(erro_media_conv_save)

%Erro médio entre overlap_add e overlap_save
diferenca = abs(overlap_add - overlap_save);
erro_media_add_save= mean(diferenca); % erro médio.
display(erro_media_add_save)

%% Listen to, plot and save your results

%sound(y,Fs);
%sound(convolution,Fs);
%sound(overlap_add,Fs);
%sound(overlap_save,Fs);

 figure
 hold on
 pspectrum(y,Fs)
 pspectrum(convolution,Fs)
 hold off

 figure
 hold on
 pspectrum(y,Fs)
 pspectrum(overlap_add,Fs)
 hold off

 figure
 hold on
 pspectrum(y,Fs)
 pspectrum(overlap_save,Fs)
 hold off


audiowrite('convolution_{Lindeberg_200051733}.wav',convolution,Fs)
audiowrite('overlap_add_{Lindeberg_200051733}.wav',overlap_add,Fs)
audiowrite('overlap_save_{Lindeberg_200051733}.wav',overlap_save,Fs)