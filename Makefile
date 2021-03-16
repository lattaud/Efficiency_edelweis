ObjSuf        = o
SrcSuf        = cxx
IncSuf        = h
ExeSuf        = exe
ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     :=  $(shell root-config --libs)  -lSpectrum -lTreePlayer -lMLP -lTMVA -lMinuit -lHtml -lThread -lXMLIO -lMathMore -lMinuit2

CXX           = /usr/bin/g++ 
COMP          = g

CXXFLAGS      =  -Wall -fPIC -g -Wno-deprecated  


CXXFLAGS     += -$(COMP)
CXXFLAGS     += $(ROOTCFLAGS) 

LIBS          = $(ROOTLIBS) $(SYSLIBS)   
 


#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
 
all :  NEPAL_Template_Select.exe BuildSimuCoincArgs.exe AnalyzeSimuCoincLOOP.exe Analyze_Coinc_cleaned_buff.exe BuildSimuCoincArgs_Batch.exe   Analyze_Coinc_cleaned.exe Analyze_Coinc_cleaned_BATCH.exe

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


.SUFFIXES: .$(ObjSuf) .$(SrcSuf) .$(DllSuf) .$(ExeSuf) .C
.PRECIOUS: %.o 

%.$(ExeSuf): %.$(SrcSuf) $(objects)
	$(CXX) $(CXXFLAGS) $(objects) $< -o $@ $(LIBS)

%.$(ObjSuf): %.C %.$(IncSuf)
	$(CXX) $(CXXFLAGS) -c $< -o $@
	@echo "$@ done..."

clean:
	rm -f *.exe *.o
 
