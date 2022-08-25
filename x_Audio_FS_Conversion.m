


function [Resized_audio] = x_Audio_FS_Conversion(Original_audio, Original_Sampling_frequency, Intended_Sampling_Frequency)

t = 0:1/Original_Sampling_frequency: (length(Original_audio)/Original_Sampling_frequency)-1/Original_Sampling_frequency;
FS = Intended_Sampling_Frequency;
t_extended = 0:1/FS:max(t);
Resized_audio = interp1(t,Original_audio,t_extended);
end





