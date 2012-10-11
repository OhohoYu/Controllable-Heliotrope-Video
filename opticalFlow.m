 function [vx,vy,warpI2]=opticalFlow(im1,im2)
 im1 = im1(1:2:end,1:2:end,:);
 im2 = im2(1:2:end,1:2:end,:);

 % code from demoflow.m
 
 alpha = 0.012;
 ratio = 0.75;
 minWidth = 20;
 nOuterFPIterations = 7;
 nInnerFPIterations = 1;
 nSORIterations = 30;

 para = [alpha,ratio,minWidth,nOuterFPIterations,nInnerFPIterations,nSORIterations];

 % end code from demoflow.m
 [vx,vy,warpI2]=Coarse2FineTwoFrames(im1,im2, para);
 
 end