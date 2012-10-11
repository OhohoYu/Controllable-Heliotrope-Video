How to run Labs7.m

Labs7.m takes the following inputs (in the correct order):
	1. path - the location of the images to use
	2. prefix - the base name of the images
	3. first - the starting number
	4. last - the end number
	5. digits - the number of digits
	6. suffix - the file type
	7. start - the image of the first node

All are the same as the origional load_sequence provided, except the last input 'start' which is used to determine the picture at which to begin the program. This number must be in the range specified by first and last. EX: if first is 14 and last is 32, start must be in the range 1-19 inclusive. In general, start must be in the range 1 to (last - first) + 1. This is also true if first is 0.

Once the program begins, the image the user specified will be displayed, along with a crosshair. The user then left clicks on the image to designate line segments. When the user is done, enter is pressed to finish point selection.

The output movie is saved as 'video.avi' and the actual path is plotted and saved as 'actual_path.jpg'


Included files:

folder: Lego Flow - contains the pictures of the lego figure used for the advanced section

actual_path.jpg - the user picked path and the actual path of the origional pictures

video.avi - the output video of the origional pictures

lego output.avi - the output of th lego path

Labs7.pdf - commented source code

General Comments and Advanced Section Comments.pdf - comments on both sections

gong.wav - the sound file used for the advaced section

Labs7.m - the orgional source code
load_sequence.m - the altered load sequence file

Readme.txt - this readme

