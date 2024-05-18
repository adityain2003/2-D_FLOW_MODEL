clc
clear
load("HEAD_ALL_IMPORT.mat")
writerObj=VideoWriter('movie_well','Uncompressed AVI')
writerObj.FrameRate=5;
open(writerObj)

for I=1:121
HEAD_I=HEAD_ALL(1+(I-1)*100:I*100,:)
subplot(1,2,1)
surf(HEAD_I);
m=min(min(HEAD_I))
el=0;
az=0;
view(az,el)
% text(100,0,'hello')
% text(10,99,1,'Head')
zlim([93 100])
xlabel('x axis (in m)')
zlabel('z axis (in m)')
subtitle('_____')
title(5*(I-1))
ax = gca;
ax.TitleHorizontalAlignment = 'right';

subplot(1,2,2)
surf(HEAD_I);
el=90;
az=00;
view(az,el)
% text(0,-10,0,'Head')
zlim([93 100])
xlabel('x axis (in m)')
ylabel('y axis (in m)')
title(100-min(min(HEAD_I)))
subtitle('_____')
ax = gca;
ax.TitleHorizontalAlignment = 'right';
frame(I)=getframe(gcf);
writeVideo(writerObj,frame(I))


end 

close(writerObj)

% movie(gcf,frame,1)

