#ifndef PARSELISTCOMMAND_H
#define PARSELISTCOMMAND_H
/*
 *  parselistcommand.h
 *  Mothur
 *
 *  Created by westcott on 2/24/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "command.hpp"
#include "groupmap.h"
#include "inputdata.h"
#include "listvector.hpp"

/***************************************************************************************/

class ParseListCommand : public Command {
	
public:
	ParseListCommand(string);
	ParseListCommand();	
	~ParseListCommand() {}
	
	vector<string> setParameters();
	string getCommandName()			{ return "parse.list";				}
	string getCommandCategory()		{ return "OTU-Based Approaches";	}
	
	string getHelpString();	
    string getOutputPattern(string);	
	string getCitation() { return "http://www.mothur.org/wiki/Parse.list"; }
	string getDescription()		{ return "parses a list file by group"; }
	
	int execute(); 
	void help() { m->mothurOut(getHelpString()); }	
	
private:
	int parse(ListVector*);
		
	ListVector* list;
	GroupMap* groupMap;
    CountTable ct;
	
	ofstream out;
	string outputDir, listfile, groupfile, label, countfile;
	set<string> labels;
	bool abort, allLines;
	vector<string> outputNames;

};

/***************************************************************************************/

#endif

