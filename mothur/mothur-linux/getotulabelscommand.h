#ifndef Mothur_getotulabelscommand_h
#define Mothur_getotulabelscommand_h

//
//  getotulabelscommand.h
//  Mothur
//
//  Created by Sarah Westcott on 5/21/12.
//  Copyright (c) 2012 Schloss Lab. All rights reserved.
//


#include "command.hpp"
#include "inputdata.h"
#include "listvector.hpp"
#include "sharedrabundvector.h"

/**************************************************************************************************/

class GetOtuLabelsCommand : public Command {
public:
    GetOtuLabelsCommand(string);
    GetOtuLabelsCommand();
    ~GetOtuLabelsCommand(){}
    
    vector<string> setParameters();
    string getCommandName()			{ return "get.otulabels";          }
    string getCommandCategory()		{ return "OTU-Based Approaches";	} 
    
	string getHelpString();	
    string getOutputPattern(string);	
    string getCitation() { return "http://www.mothur.org/wiki/Get.otulabels"; }
    string getDescription()		{ return "Can be used with output from classify.otu, otu.association, or corr.axes to select specific otus."; }
    
    int execute(); 
    void help() { m->mothurOut(getHelpString()); }	
    
private:
    bool abort;
    string outputDir, accnosfile, constaxonomyfile, otucorrfile, corraxesfile, listfile, sharedfile, label;
    vector<string> outputNames;
    set<string> labels;
    ListVector* list;
    vector<SharedRAbundVector*> lookup;
    
    int readClassifyOtu();
    int readOtuAssociation();
    int readCorrAxes();
    int readList();
    int readShared();
    int getListVector();
    int getShared();
};

/**************************************************************************************************/






#endif
