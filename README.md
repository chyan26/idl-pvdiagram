# PV Diagram for IDL

 Introduction

For long time, I use IDL as my analysis tool for OIR astronomy.  IDL provided wonderful support on FIT IO, WCS mapping and plotting libraries.  The most important thing is that it allows users to expanded the function.  Recently, I used data from radio telescope to understand the gas dynamics and how they interact with molecular cloud.  PV diagram is a very important tool in radio astronomy.  However, there is no such a tool in the world of IDL.  Researchers used command line tools like GILDAS, MIRIAD and CASA to analyze the FITS cube.  However, none of these packages fits all following needs.

1. Analyzing the PV diagram of sub-region in a big radio map.
2. GUI window to interact with a user.
3. Customizing the output result in the format desired.

Therefore, I implemented this PV diagram toolkit and released the source code for IAA users.  Hope you can find this package useful.   

System Requirement

To use this program, of course, you need IDL.  Please contact your system administrator for the installation detail.  This program requires the IDL Astronomy User's Library.  You can download the source from this website.   It is also suggested that user install Cyote Library because this library provides good programs to make color plots.  

Installation

Download the IDL source code here and put this file to a directory.  Make sure this path is defined in your IDL !path variable.  


Syntax

 CALLING SEQUENCE:
       PVDIAGRAM,fits=fits,[vel_range=[min,max],scale=scale]

 INPUTS:
       FITS   - FITS cube exported from CASA or GILDAS

 OPTIONAL INPUTS:
       VEL_RANGE - velocity range for the PV diagram.
       SCALE     - scaling factor for integrated map       
       XYSTART   - positiion of begining point for PV analysis
       XYEND     - positiion of ending point for PV analysis


Example

Executing the program from prompt.  You will see a window like this.  The is the integrated intensity map of your FITS cube.  There are three buttons, Slice, Clear Line and Exit.  

 

You can use cursor to draw a line or input coordinates to draw a line.  The line will be shown like the following figure. 

When you are ready, click Slice button and the PV diagram window will show up as shown below.  This PV diagram can be saved in different formats, you can use them in your presentation or writing. 




Future Development

The design of this program follows a standard design model of software engineering.  That means the GUI, PV diagram and user inputs are all independent.   The command line takes the user inputs and starting the GUI.  GUI pass the position coordinates into a subroutine that generates the PV diagram.  A more advanced future is non-linear PV cut.   I will be working on this soon.   
