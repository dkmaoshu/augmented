# Microsoft Developer Studio Project File - Name="main" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Console Application" 0x0103

CFG=main - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "main.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "main.mak" CFG="main - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "main - Win32 Release" (based on "Win32 (x86) Console Application")
!MESSAGE "main - Win32 Debug" (based on "Win32 (x86) Console Application")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
F90=df.exe
RSC=rc.exe

!IF  "$(CFG)" == "main - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "Release"
# PROP Intermediate_Dir "Release"
# PROP Target_Dir ""
# ADD BASE F90 /compile_only /nologo /warn:nofileopt
# ADD F90 /compile_only /nologo /warn:nofileopt
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /c
# ADD CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /c
# ADD BASE RSC /l 0x804 /d "NDEBUG"
# ADD RSC /l 0x804 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib /nologo /subsystem:console /machine:I386
# ADD LINK32 kernel32.lib /nologo /subsystem:console /machine:I386

!ELSEIF  "$(CFG)" == "main - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "Debug"
# PROP Intermediate_Dir "Debug"
# PROP Target_Dir ""
# ADD BASE F90 /check:bounds /compile_only /debug:full /nologo /traceback /warn:argument_checking /warn:nofileopt
# ADD F90 /check:bounds /compile_only /debug:full /nologo /traceback /warn:argument_checking /warn:nofileopt
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /GZ  /c
# ADD CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /GZ  /c
# ADD BASE RSC /l 0x804 /d "_DEBUG"
# ADD RSC /l 0x804 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib /nologo /subsystem:console /debug /machine:I386 /pdbtype:sept
# ADD LINK32 kernel32.lib /nologo /subsystem:console /incremental:no /debug /machine:I386 /pdbtype:sept

!ENDIF 

# Begin Target

# Name "main - Win32 Release"
# Name "main - Win32 Debug"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat;f90;for;f;fpp"
# Begin Source File

SOURCE=..\facility\average2.f90
NODEP_F90_AVERA=\
	".\Debug\physics.mod"\
	".\Debug\RiemannSolver.mod"\
	".\Debug\solution.mod"\
	
# End Source File
# Begin Source File

SOURCE=..\facility\boundary_conditions.f90
NODEP_F90_BOUND=\
	".\Debug\infinitive_boundary.mod"\
	".\Debug\periodic_boundary.mod"\
	".\Debug\reflective_boundary.mod"\
	
# End Source File
# Begin Source File

SOURCE=..\facility\flux.f90
NODEP_F90_FLUX_=\
	".\Debug\physics.mod"\
	".\Debug\RiemannSolver.mod"\
	".\Debug\solution.mod"\
	
# End Source File
# Begin Source File

SOURCE=..\facility\grid_parameters.f90
# End Source File
# Begin Source File

SOURCE=..\facility\in_put.f90
NODEP_F90_IN_PU=\
	".\Debug\solution.mod"\
	
# End Source File
# Begin Source File

SOURCE=..\Initial_values\boundary_conditions\infitive_boundary.f90
NODEP_F90_INFIT=\
	".\Debug\solution.mod"\
	
# End Source File
# Begin Source File

SOURCE="..\facility\main_Euler(m).f90"
NODEP_F90_MAIN_=\
	".\Debug\average.mod"\
	".\Debug\boundary_conditions.mod"\
	".\Debug\flux.mod"\
	".\Debug\in_put.mod"\
	".\Debug\out_put.mod"\
	".\Debug\output_show.mod"\
	".\Debug\solution_reconstruction.mod"\
	".\Debug\time_step.mod"\
	
# End Source File
# Begin Source File

SOURCE=..\facility\out_put.f90
NODEP_F90_OUT_P=\
	".\Debug\solution.mod"\
	
# End Source File
# Begin Source File

SOURCE=..\facility\output_4_show.f90
NODEP_F90_OUTPU=\
	".\Debug\solution.mod"\
	
# End Source File
# Begin Source File

SOURCE=..\Initial_values\boundary_conditions\periodic_boundary.f90
NODEP_F90_PERIO=\
	".\Debug\solution.mod"\
	
# End Source File
# Begin Source File

SOURCE=..\facility\physics.f90
NODEP_F90_PHYSI=\
	".\Debug\grid_and_parameters.mod"\
	
# End Source File
# Begin Source File

SOURCE=..\Initial_values\boundary_conditions\reflective_boundary.f90
NODEP_F90_REFLE=\
	".\Debug\solution.mod"\
	
# End Source File
# Begin Source File

SOURCE=..\facility\RiemannSolver.f90
NODEP_F90_RIEMA=\
	".\Debug\grid_and_parameters.mod"\
	".\Debug\physics.mod"\
	
# End Source File
# Begin Source File

SOURCE=..\facility\solution1.f90
NODEP_F90_SOLUT=\
	".\Debug\grid_and_parameters.mod"\
	
# End Source File
# Begin Source File

SOURCE=..\facility\step_reconstruction5.f90
DEP_F90_STEP_=\
	".\Debug\physics.mod"\
	".\Debug\solution.mod"\
	
# End Source File
# Begin Source File

SOURCE=..\facility\Time_step.f90
NODEP_F90_TIME_=\
	".\Debug\physics.mod"\
	".\Debug\solution.mod"\
	
# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl;fi;fd"
# End Group
# Begin Group "Resource Files"

# PROP Default_Filter "ico;cur;bmp;dlg;rc2;rct;bin;rgs;gif;jpg;jpeg;jpe"
# End Group
# End Target
# End Project
