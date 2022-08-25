
clc
clear

[Aud1, SampRate] = audioread('Audio_File_1.wav');
%Foi preciso eliminar estas amostras para que o sinal possa ser codificado
%com o codigo de Hamming (7,4)
Aud1 = Aud1(2:end-2);

%S1

%separacao das amplitudes iguais
a = unique(Aud1);
%probabilidade de cada amplitude
prob = histc(Aud1(:),a)./length(Aud1);
Tab_prob = [a , prob.*100];
%codificacao de huffman
[DICCIONARIO, compmed] = huffmandict(a, prob);
S1 = huffmanenco(Aud1, DICCIONARIO);
%histograma
figure(1)
%PQ = 2/(2^8)
histogram(Aud1, 'BinWidth', 2/(2^8));
title('Histograma das amplitudes do sinal de audio 1')
xlabel('Amplitude')
ylabel('Frequencia')
%entropia                
entrop = 0;
for i = 1:length(prob)
     if prob(i) ~= 0
        entrop = entrop + (prob(i)*log2(1/prob(i)));
     end
end
%comprimento medio do codigo -> compmed
%eficiencia
eficiencia = entrop/compmed;


%S2

% codigo hamming (7,4)
% k = 4;  
% n = 7;
 m = 3;    %m = n - k 

[MATRIZ_H,MATRIZ_G,~,~] = hammgen(m);
%Separar o fluxo em 4
u = reshape(S1,4,[]); 
u = transpose(u);
%row|col
%multiplicacao da primeira linha ate a quarta
ug1 = u(1:end,1).*MATRIZ_G(1,1:end);
ug2 = u(1:end,2).*MATRIZ_G(2,1:end);
ug3 = u(1:end,3).*MATRIZ_G(3,1:end);
ug4 = u(1:end,4).*MATRIZ_G(4,1:end);
%soma XOR das linhas para obter a mensagem u codificada
n1 = xor(ug1,ug2);
n2 = xor(n1,ug3);
ucod = xor(n2, ug4);
%juntar tudo numa linha
S2 = reshape(ucod',1,[]);

%S3
load('CHANNEL_ERRORS.mat');
CHANNEL_ERRORS = CHANNEL_ERRORS(1:end-35);
S3 = xor(CHANNEL_ERRORS,S2);

%S4
ht = transpose(MATRIZ_H);
v = reshape(S3,7,[]);
v = transpose(v);
%PADRAO DE ERRO
pe1 = v(1:end,1)-ucod(1:end,1);
pe2 = v(1:end,2)-ucod(1:end,2);
pe3 = v(1:end,3)-ucod(1:end,3);
pe4 = v(1:end,4)-ucod(1:end,4);
pe5 = v(1:end,5)-ucod(1:end,5);
pe6 = v(1:end,6)-ucod(1:end,6);
pe7 = v(1:end,7)-ucod(1:end,7);
PADRAO_ERRO = [abs(pe1), abs(pe2), abs(pe3), abs(pe4), abs(pe5), abs(pe6), abs(pe7)];
%Calculos para o sindrome v*Ht = s
vht1 = v(1:end,1).*ht(1,1:end);
vht2 = v(1:end,2).*ht(2,1:end);
vht3 = v(1:end,3).*ht(3,1:end);
vht4 = v(1:end,4).*ht(4,1:end);
vht5 = v(1:end,5).*ht(5,1:end);
vht6 = v(1:end,6).*ht(6,1:end);
vht7 = v(1:end,7).*ht(7,1:end);
n1 = xor(vht1,vht2);
n2 = xor(n1,vht3);
n3 = xor(n2,vht4);
n4 = xor(n3,vht5);
n5 = xor(n4,vht6);
SINDROME = xor(n5,vht7);
%procura as linhas do sindrome iguais as linhas em ht
[erros,ind] = ismember(SINDROME,ht,'rows');
%devolve a posicao do erro dentro da mensagem
[row,~,col_erro] = find(ind);
%Corrige o erro
linhaserro = v(row(1:end),1:end);
for i = 1:length(linhaserro)
    linhaserro(i,col_erro(i)) = not(linhaserro(i,col_erro(i)));
end
S4 = v;
S4(row(1:end),1:end) = linhaserro(1:end,1:end);
% b1 b2 b3 m1 m2 m3 m4
S4(:,1) = [];    %elimina coluna1 -> b1
S4(:,1) = [];    %elimina coluna2 -> b2
S4(:,1) = [];    %elimina coluna3 -> b3

%juntar tudo numa linha
S4 = reshape(S4',1,[]);
S4 = S4';

%S5
S4 = double(S4);
S5 = huffmandeco(S4, DICCIONARIO);
%sound(S5,SampRate);

%S6
fs = 4000000;
S5_4M = x_Audio_FS_Conversion(S5, SampRate, fs);
fc_AM = 1000000;
S6 = ammod(S5_4M, fc_AM, fs);
n = 2^nextpow2(length(S6)); 
n9 = n;
S6_espectro = fft(S6, n)./n; 
S6_espectro = fftshift(S6_espectro);
f_AM = (fs/n)*(-n/2:n/2-1); 
figure(3)
plot(f_AM, 20*log10(abs(S6_espectro))); 
title('Espectro do sinal S6')
xlabel('Frequencia (Hz)')
ylabel('Potencia (dBW)')
%DSB
%xlim([1000000-22050 1000000+22050]);
%LSB
xlim([1000000-22050 1000000]);


%S7
[Aud2, SampRate] = audioread('Audio_File_2.wav');
Aud2_4M = x_Audio_FS_Conversion(Aud2, SampRate, fs);
%Modulacao FM
fc_FM = 1.5*10^6;
fvar = 75*10^3;
S7 = fmmod(Aud2_4M, fc_FM, fs, fvar);
%Grafico do espectro do sinal FM
n = 2^nextpow2(length(S7)); % este n e igual ao n9
S7_espectro = fft(S7, n)./n; 
S7_espectro = fftshift(S7_espectro);
f_FM = (fs/n)*(-n/2:n/2-1); 
figure(4)
plot(f_FM, 20*log10(abs(S7_espectro))); 
title('Espectro do sinal S7')
xlabel('Frequencia (Hz)')
ylabel('Potencia (dBW)')
%xlim([fc_FM-194100/2 fc_FM+194100/2]);

%S8

%Filtro passa-banda AM
%Bw_AM = 2 frequencia max sinal Aud1
%BW_AM = 44.1*10^3;
BW_AM = 40*10^3;
Filtro_AM = zeros(1, length(S6_espectro));
%DSB
%Filtro_AM( f_AM < fc_AM+BW_AM/2 & f_AM > fc_AM-BW_AM/2) = 1;
%LSB
Filtro_AM( f_AM < fc_AM & f_AM > fc_AM-BW_AM/2) = 1;
%Filtro passa-banda FM
%Regra de carson -> banda transmissao = 2 * (fvar + fc)
BW_FM = 2*(fvar + 22050);
Filtro_FM = zeros(1, length(S7_espectro));
Filtro_FM(f_FM < fc_FM+BW_FM/2 & f_FM > fc_FM-BW_FM/2) = 1;
%Aplicacao dos filtros
S6_espectro_filtrado = S6_espectro.*Filtro_AM;
S7_espectro_filtrado = S7_espectro.*Filtro_FM;
%Multiplexacao em frequencia
S8_espectro = S6_espectro_filtrado + S7_espectro_filtrado;
%IFFT
S8 = real(ifft(fftshift(S8_espectro))*n9);


%S9 & S10

%Calculo da fft
n8 = 2^nextpow2(length(S8)); 
S9S10_espectro = fft(S8, n8)./n8; 
S9S10_espectro = fftshift(S9S10_espectro);
%filtragem dos sinais
S9_espectro = S9S10_espectro.*Filtro_AM;
S10_espectro = S9S10_espectro.*Filtro_FM;
%ifft sinais
S9 = real(ifft(fftshift(S9_espectro))*n8*2*2*2); 
S10 = real(ifft(fftshift(S10_espectro))*n8);



%S11 & S12

%Desmodulacao
S11_4M = amdemod(S9, fc_AM, fs);
S12_4M = fmdemod(S10, fc_FM, fs, fvar);
%x_Audio_FS_Conversion
S11 = x_Audio_FS_Conversion(S11_4M, fs, SampRate);
S12 = x_Audio_FS_Conversion(S12_4M, fs, SampRate);
%Remover o ruido final
S11 = S11(1:length(Aud1));
S12 = S12(1:length(Aud2));

%%%%%%%%%%%%%%%%%%%%%%%%%%
%Graficos auxiliares
%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Calculo da fft sinal 4MHz
% n = 2^nextpow2(length(S5_4M)); 
% S54M_espectro = fft(S5_4M, n)./n; 
% S54M_espectro = fftshift(S54M_espectro);
% f_AM4M = (fs/n)*(-n/2:n/2-1); 
% figure(2)
% plot(f_AM4M, 20*log10(abs(S54M_espectro))); 
% title('Espectro do sinal S5')
% xlabel('Frequencia (Hz)')
% ylabel('Potencia (dBW)')
    %Graficos dos sinais ja separados
    %Utilizados para ver o ruido final a ser eliminado
% figure(6)
% plot(S9);
% figure(7)
% plot(S10);
    %Graficos dos filtros
% figure(10)
% plot(f_AM, Filtro_AM);
% figure(11)
% plot(f_FM, Filtro_FM);

    %Graficos comparacao sinal original com o processado
figure(8)
subplot(2,1,1)
plot(Aud1);
subplot(2,1,2)
plot(S11);
figure(9)
subplot(2,1,1)
plot(Aud2);
subplot(2,1,2)
plot(S12);
%%%%%%%%%%%%%%%%%%%%%%%%%%
%Questoes relatorio
%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matriz de paridade P
    %P(k m)
    P = MATRIZ_G(1:end, 1:3);
%Largura de banda do sinal de audio 1
n = 2^nextpow2(length(Aud1)); 
S5_espectro = fft(S5, n)./n; 
S5_espectro = fftshift(S5_espectro);
frequencias_aud1 = (SampRate/n)*(-n/2:n/2-1); 
figure(12)
plot(frequencias_aud1, 20*log10(abs(S5_espectro))); 
title('Espectro do sinal S5 (Audio1)')
xlabel('Frequencia (Hz)')
ylabel('Potencia (dBW)')
xlim([0 22050]);
LB_Aud1 = max(frequencias_aud1); %22050Hz
%Largura de banda do sinal AM
LB_AM = 2*LB_Aud1; %44100Hz
%Largura de banda do sinal de audio 2
n = 2^nextpow2(length(Aud2));
frequencias_aud2 = (SampRate/n)*(-n/2:n/2-1); 
LB_Aud2 = max(frequencias_aud2); %22050Hz
%Largura de banda do sinal FM
LB_FM = 2*(75000 + 22050); % 194.1kHz
%Graficos dos sinais S8, S6, S7
figure(13)
subplot(2,2,3)
plot(f_AM, 20*log10(abs(S6_espectro)));
title('Espectro do sinal S6')
xlabel('Frequencia (Hz)')
ylabel('Potencia (dBW)')
xlim([0 2000000]);
subplot(2,2,4)
plot(f_FM, 20*log10(abs(S7_espectro)));
title('Espectro do sinal S7')
xlabel('Frequencia (Hz)')
ylabel('Potencia (dBW)')
xlim([0 2000000]);
subplot(2,2,[1,2])
%Sinal S8 depois de ter os filtros aplicados
plot(f_FM, 20*log10(abs(S8_espectro)));
title('Espectro do sinal S8')
xlabel('Frequencia (Hz)')
ylabel('Potencia (dBW)')
xlim([0 2000000]);
%Sinal S8 antes de ter os filtros aplicados
f = (fs/n8)*(-n8/2:n8/2-1); 
figure(5)
plot(f, 20*log10(abs(S9S10_espectro))); 
title('Espectro do sinal S8 recebido')
xlabel('Frequencia (Hz)')
ylabel('Potencia (dBW)')
%%%%%%%%%%%%%%%%%%%%%%%%%%
%ficheiros
%%%%%%%%%%%%%%%%%%%%%%%%%%
save('FT2182076_1.mat', 'S1', 'S2', 'S3', 'S4', 'S5', 'S6', 'S11', 'S12', 'DICCIONARIO', 'MATRIZ_G', 'MATRIZ_H', 'SINDROME', 'PADRAO_ERRO');
save('FT2182076_2.mat', 'S7');
save('FT2182076_3.mat', 'S8');
save('FT2182076_4.mat', 'S9');
save('FT2182076_5.mat', 'S10');