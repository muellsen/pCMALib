###########################################################
# CEC2005 Makefile
#
# 
# Updates:
# 03/16/2008
# $Log: makefile,v $
# Revision 1.1  2008/12/17 11:56:18  ofgeorg
# created makefile for CEC2005 Suit that will create the Example Program
# All Parameters that one can change are in the make.inc file (only one so far :)
#
#
# Georg Ofenbeck
# Computational Biophysics Lab, ETH Zurich, Switzerland
###########################################################

include make.inc

# -------------------End-macro-Defs---------------------------

ifeq ($(F90), ifort)
	CPP = -fpp
else
	CPP = -cpp
endif



# Here is the link step 
Example_CEC2005:CEC2005.o testCEC2005.o Example_CEC2005.o randgen.o
	 $(F90) -o Example_CEC2005 CEC2005.o testCEC2005.o Example_CEC2005.o randgen.o  

# Here are the compile steps
randgen.o:./randgen.f  
	 $(F90) $(CFLAGS) -c ./randgen.f   

CEC2005.o:./CEC2005.F90 randgen.o  
	 $(F90) $(CFLAGS) -c  $(CPP) ./CEC2005.F90 

CEC2005.mod:./CEC2005.o CEC2005.F90
	$(F90) $(CFLAGS) -c  $(CPP) CEC2005.F90

testCEC2005.o:./testCEC2005.F90 CEC2005.o CEC2005.mod randgen.o
	 $(F90) $(CFLAGS) -c ./testCEC2005.F90 

testCEC2005.mod:./testCEC2005.o testCEC2005.F90 CEC2005.mod randgen.o
	 $(F90) $(CFLAGS) -c ./testCEC2005.F90 

Example_CEC2005.o:Example_CEC2005.F90 CEC2005.o testCEC2005.o randgen.o
	 $(F90) $(CFLAGS) -c  Example_CEC2005.F90 
	 
	 
	 
	 
# This entry allows you to type " make clean " to get rid of
# all object and module files 
clean:
	rm -f  *.o *.mod *.exe
  
