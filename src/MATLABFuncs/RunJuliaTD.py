import numpy as np
import sys
import julia
import json
import scipy.io

def TDInterface(x,y,z,P1,P2,P3,DssVec,DdsVec,DnVec,nu,mu,DispFlag,StressFlag,HSFlag):

	#Compute something in Julia (from a module)
	j = julia.Julia()
	from julia import MyModule
	(ExxDn,EyyDn,EzzDn,ExyDn,ExzDn,EyzDn,
	ExxDss,EyyDss,EzzDss,ExyDss,ExzDss,EyzDss,
	ExxDds,EyyDds,EzzDds,ExyDds,ExzDds,EyzDds,
	UxDn,UyDn,UzDn,UxDss,UyDss,UzDss,UxDds,UyDds,UzDds)=MyModule.TD(x,y,z,P1,P2,P3,DssVec,DdsVec,DnVec,nu,mu,DispFlag,StressFlag,HSFlag)

	#Or save as .Mat file
	scipy.io.savemat('PythonOutput.mat', dict(
	ExxDn=ExxDn,EyyDn=EyyDn,EzzDn=EzzDn,ExyDn=ExyDn,ExzDn=ExzDn,EyzDn=EyzDn,
	ExxDss=ExxDss,EyyDss=EyyDss,EzzDss=EzzDss,ExyDss=ExyDss,ExzDss=ExzDss,EyzDss=EyzDss,
	ExxDds=ExxDds,EyyDds=EyyDds,EzzDds=EzzDds,ExyDds=ExyDds,ExzDds=ExzDds,EyzDds=EyzDds,
	UxDn=UxDn,UyDn=UyDn,UzDn=UzDn,UxDss=UxDss,UyDss=UyDss,UzDss=UzDss,UxDds=UxDds,UyDds=UyDds,UzDds=UzDds))