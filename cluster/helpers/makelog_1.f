* This progam generates a logfile for SciAn Visualisation

	program makelog

	include 'const.inc'
	include 'timep.inc'
	include 'filep.inc'
	include 'atomp.inc'

	integer i,maxi
	character*40 fout

	call readdata

	i=2
	call getarg(i,fout)
	
	maxi=ntim/ftxyz

	open (1,file=fout)

	write (1,*) 'Version 1.2 Alpha'
	write (1,*) 'time Tue Jan  4 10:30:19 2000'
	write (1,*) 
	write (1,*) 'window SciAn'
	write (1,*) 'show recorders'
	write (1,*) 'window Recorder_Drivers'
	write (1,*) 'activate recorder@scrsave'
	write (1,*) 'close'
	write (1,*) 
	write (1,*) 'window Datasets'
	write (1,*) 'visualize dataset@vis\\.atoms'
	write (1,*) 'window Visualization_1'
	write (1,*) 'locate 193 907 62 470'
	write (1,*) 
	write (1,*) 'window Datasets'
	write (1,*) 'visobjectas dataset@TimeFile'
	write (1,*) 'window Visualize_As_1'
	write (1,*) 'set value Number_1 1'
	write (1,*) 'drag Number_Template_1'
	write (1,*) 'window Visualization_1'
	write (1,*) 'drop 588 348'
	write (1,*) 
	write (1,*) 'window Visualize_As_1'
	write (1,*) 'close'
	write (1,*) 'window Datasets'
	write (1,*) 'close'
	write (1,*) 
	write (1,*) 'window Visualization_1'
	write (1,*) 'show controls Number_2'
	write (1,*) 'window Number_2'
	write (1,*) 'set value Format_Text "Time: %6.2f ps"'
	write (1,*) 'set value Slice 1'
	write (1,*) 'set value j_Active 1'
	write (1,*) 'set value i_Active 1'
	write (1,*) 'set value Slice_0 0'
	write (1,*) 'set value Slice_1 0'
	write (1,*) 'close'
	write (1,*) 
	write (1,*) 'window Visualization_1'
	write (1,*) 'set bounds Numeric_readout_text_1 [220 430 5 35]'
	write (1,*) 'set size Numeric_readout_text_1 14'
	write (1,*) 'annotation Annotation_1 [-10 190 5 35]'
	write (1,'(A,I4,A)') ' set value Annotation_1 "',int(att),' K"'
	write (1,*) 'set size Annotation_1 14'
	write (1,*) 
	write (1,*) 'show controls Point_Cloud_1'
	write (1,*) 'window Point_Cloud_1'
	write (1,*) 'set value Dots 1'
	write (1,*) 'set value Bigger_Dots 1'
	write (1,*) 'set value Bounds 1'
	write (1,*) 'set value Draw_Outline 1'
	write (1,*) 
	write (1,*) 'set value Color 1'
*	write (1,*) 'edit palette dataset@vis\\.atoms'
*	write (1,*) 'window vis\\.atoms_Palette'
*	write (1,*) 'set value Field_Max "106"'
*	write (1,*) 'close'
	write (1,*) 'window Point_Cloud_1'
	write (1,*) 'close'
	write (1,*) 
	write (1,*) 'window Visualization_1'
	write (1,*) 'hide panel'
	write (1,*) 'locate 800 1200 600 1000'
	write (1,*) 'set rotation 0 [-0.943475 0.255815 0.210744]'
	write (1,*) 
	write (1,*) 'begin snapshot Observer_1'
	write (1,*) '  set variable LOCATION [-1.5276 1.39327 4.55642]'
	write (1,*) '  set variable ROTQUAT [0.973042 0.153898 '//
     $		    '0.169064 0.0304148]'
	write (1,*) '  set variable FOCUSDIST 5'
	write (1,*) 'end snapshot'
	write (1,*) 
	write (1,*) 'hide panel'
	write (1,*) 'locate 800 1200 600 1000'
	write (1,*) 
	write (1,*) 'show controls Clock_1'
	write (1,*) 'window Clock_1'
	write (1,*) 'set value Animation_Optimization 0'
	write (1,*) 'set value Cycle_the_clock 0'
	write (1,*) 'set value Advance_Type 3'
	write (1,*) 'set value Delta_Time_Box "1"'
	write (1,*) 'set value Frame_Rate "1"'
	write (1,*) 'set value Speed_Control 1'
	write (1,*) 
	write (1,*) 'window Visualization_1'
	write (1,*) 'begin recording 2000'
	write (1,'(A,F10.5)') '   record ',0.01+real(maxi)/30.0
	write (1,*) 'end recording'
	write (1,*) 
	write (1,*) 'window Clock_1'
	write (1,*) 'set value Speed_Control 0'
	write (1,*) 'quit'

	close (1)

	end
