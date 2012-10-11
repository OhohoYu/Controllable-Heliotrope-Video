% get imsequence. output images.
% use sequence to take a frame one at a time and creat an avi movie

% load audio
[y, Fs] = wavread('gong');
player = audioplayer(y, Fs);
% use play(player); to play the file. add timing to get the sound to play once for each line segment

[x,y,z] = size(sequence);
frames = zeros(z);
for i=1:z
	frames(i) = im2frame(sequence(:,:,i)); 
end
fps=5;
n=1;
% output the video without audio
movie2avi(frames, 'video.avi', 'fps', fps);

movie(frames,n,fps);