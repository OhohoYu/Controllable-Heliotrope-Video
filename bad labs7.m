 
 

 
 % S = exp (-D/sigma) this is the simmilarity, given in class
 % there is a function called 'shortestpath'
 % conncomp
 % using biograph, or start with graph if using a matrix
function output = labs7 
 
 % 1. load files
 % 2. compute simmilarity
 % 3. based on a threshold, connect them in a graph of transitions
 % 4. compute optical flow along each edge
 % 5. accept user line
 % 6. compute shortest path between all edges in line
 % 7. output movie
 
 path ='C:\Users\Rob\Desktop\Study Abroad Classes\Computational Photography\Labs7\pictures';
 prefix = 'gjbLookAtTarget_';
 digits = 4;
 first = 1;
 last = 4;
 suffix = 'jpg';
 
 	load_sequence(path, prefix, first, last, digits, suffix);
	imgs=ans;
	imgs = im2double(imgs);
	userPoints = getInput(imgs(:,:,1))
    S=similarity(imgs);
	S=sparse(S);
 %   [S, pred] = graphminspantree(S, 1);
    view(biograph(S));
    
    % matrix is of the form (vx,vy,opticalflowimage,im2,im1) where im1
    % gives the list of all possible images, and im2 gives the images
    % connected to each im1
    [x,y,z] = size(imgs);
    j=0;
    for k=1:z
        % for each node, create a list of nodes it can reach, remove any
        % selfe reference, and sort the list.
        nodes = graphtraverse(S, k,'Depth', 1);
        nodes = nodes(2:end);
        nodes = sort(nodes);
        for m=1:length(nodes)
           [vx,vy,flowImages] = opticalFlow(imgs(:,:,k),imgs(:,:,nodes(m)));
		   j=j+1;
           flowNodes(k,nodes(m))=j;
           VX(j,:,:)=vx;
           VY(j,:,:)=vy;
           flowIms(:,:,j)=flowImages;
       end
    end
    
    
    
    %image 1 to desired second point
    %search tree for distance, then search for position
    
    X1=round(userPoints(1,1));
    X2=round(userPoints(2,1));
    Y1=round(userPoints(1,2));
    Y2=round(userPoints(2,2));
    
    Xdir=X2-X1;
    Ydir=Y2-Y1;
    
    nodes = graphtraverse(S, k,'Depth', 1);
    nodes = nodes(2:end);
    nodes = sort(nodes);
    
    % initialized to high values so the comparison works, acts as threshold
    CurX=100;
    CurY=100;
    % gives 
%    for i=1:z make sure to change 1's to i's
        for j=1:z
            if(1~=j && flowNodes(1,j) ~= 0)
                index=flowNodes(1,j);
                vx=VX(index,X1,Y1);
                vy=VY(index,X1,Y1); 
                Xchange=vx-Xdir;
                Ychange=vy-Ydir;
                if Xchange < CurX
                    if Ychange < CurY
                        CurX = Xchange
                        CurY = Ychange
                        curNode=j
                    end
                end
            end
        end
        [dist, path, pred]=graphshortestpath(S,1,curNode);
        
        path
        length(path)
        
        x=length(path);
        [w,y,z]=size(imgs);
        w=w/2;
        y=y/2;
        sequence = zeros(w,y,x);
        sequence(:,:,1)=imgs(1:2:end,1:2:end,1);
        vfvf='im added'
        for i=1:length(path)-1
            index=flowNodes(i,i+1);
            sequence(:,:,2)=flowIms(:,:,index);
            vfvf
        end
        fps=10;
        implay(sequence,fps);
        
%    end
    
		
			
 
 
 end 
 


 
 
 % code to read two images and compute their optical flow
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
 
function b = similarity(ims)
	[x,y,z] = size(ims);
	D=zeros(z,z);
	thresh=85;

	for i=1:z
		for j=1:z
			val =ims(:,:,i) - ims(:,:,j);
			temp=sqrt(sum(sum(val.^2)));
			if(temp<thresh)
				D(i,j)=temp;
			else
				D(i,j)=0;
			end
		end
	end
	
% return the similarity matrix	
	b=D;
 end
 
function c = getInput(im)
	f=figure;
    im = im(1:2:end,1:2:end,:);
    imshow(im,'InitialMagnification', 100);
    a = get(f,'Children');
	[x,y] = ginput;
    hold(a);
    plot(x,y);
    hold off;
    c = [x,y];
end
 

function play(path,imgs)
	fps=10;
	x=length(path);
    [w,y,z]=size(imgs);
    sequence = zeros(w,y,x);
	for i=1:x
		sequence(:,:,i) = imgs(:,:,path(i));
	 end
	implay(sequence,fps);
 end
 
