//
//  designmap.cpp
//  Mothur
//
//  Created by SarahsWork on 6/17/13.
//  Copyright (c) 2013 Schloss Lab. All rights reserved.
//

#include "designmap.h"

/************************************************************/
DesignMap::DesignMap(string file) {
    try {
        m = MothurOut::getInstance();
        defaultClass = "not found";
        read(file);
    }
	catch(exception& e) {
		m->errorOut(e, "DesignMap", "DesignMap");
		exit(1);
	}
}
/************************************************************/
int DesignMap::read(string file) {
    try {
        ifstream in;
        m->openInputFile(file, in);
        
        string headers = m->getline(in); m->gobble(in);
        vector<string> columnHeaders = m->splitWhiteSpace(headers);
        
        namesOfCategories.clear();
        indexCategoryMap.clear();
        indexNameMap.clear();
        designMap.clear();
        map<int, string> originalGroupIndexes;
        for (int i = 1; i < columnHeaders.size(); i++) {  namesOfCategories.push_back(columnHeaders[i]);  originalGroupIndexes[i-1] = columnHeaders[i]; }
        if (columnHeaders.size() > 1) { defaultClass = columnHeaders[1]; }
        
        //sort groups to keep consistent with how we store the groups in groupmap
        sort(namesOfCategories.begin(), namesOfCategories.end());
        for (int i = 0; i < namesOfCategories.size(); i++) {  indexCategoryMap[namesOfCategories[i]] = i; }
        int numCategories = namesOfCategories.size();
        
        bool error = false;
        string group;
        totalCategories.resize(numCategories);
        int count = 0;
        while (!in.eof()) {
            
            if (m->control_pressed) { break; }
            
            in >> group; m->gobble(in); 
            if (m->debug) { m->mothurOut("[DEBUG]: group = " + group + "\n"); }
            
            //if group info, then read it
            vector<string> categoryValues; categoryValues.resize(numCategories, "not found");
            for (int i = 0; i < numCategories; i++) {
                int thisIndex = indexCategoryMap[originalGroupIndexes[i]]; //find index of this category because we sort the values.
                string temp = "not found";
                in >> temp; categoryValues[thisIndex] = temp; m->gobble(in);
                
                if (m->debug) { m->mothurOut("[DEBUG]: value = " + temp + "\n"); }
                
                //do we have this value for this category already
                map<string, int>::iterator it = totalCategories[thisIndex].find(temp);
                if (it == totalCategories[thisIndex].end()) { totalCategories[thisIndex][temp] = 1; }
                else {  totalCategories[thisIndex][temp]++; }
            }
                           
            
            map<string, int>::iterator it = indexNameMap.find(group);
            if (it == indexNameMap.end()) {
                indexNameMap[group] = count;
                designMap.push_back(categoryValues);
                count++;
            }else {
                error = true;
                m->mothurOut("[ERROR]: Your design file contains more than 1 group named " + group + ", group names must be unique. Please correct."); m->mothurOutEndLine();
            }
        }
        in.close();
        
        //sanity check
        for (int i = 0; i < totalCategories.size(); i++) {
            map<string, int>::iterator it = totalCategories[i].find(namesOfCategories[i]);
            if (it != totalCategories[i].end()) { //we may have an old style design file since category name matches a value
                m->mothurOut("\n[WARNING]: Your design file has a category and value for that category named " + namesOfCategories[i] + ". Perhaps you are using an old style design file without headers? If so, please correct."); m->mothurOutEndLine();
            }
        }
        
        if (error) { m->control_pressed = true; }
        
        return 0;
    }
	catch(exception& e) {
		m->errorOut(e, "DesignMap", "read");
		exit(1);
	}
}
/************************************************************/
////groupName, returns first categories value. 
string DesignMap::get(string groupName) {
    try {
        string value = "not found";
        
        map<string, int>::iterator it2 = indexNameMap.find(groupName);
        if (it2 == indexNameMap.end()) {
            m->mothurOut("[ERROR]: group " + groupName + " is not in your design file. Please correct.\n"); m->control_pressed = true;
        }else {
            return designMap[it2->second][0];
        }
       
        return value;
    }
	catch(exception& e) {
		m->errorOut(e, "DesignMap", "get");
		exit(1);
	}
}
/************************************************************/
////categoryName, returns category values.
vector<string> DesignMap::getValues(string catName) {
    try {
        vector<string> values;
        
        map<string, int>::iterator it2 = indexCategoryMap.find(catName);
        if (it2 == indexCategoryMap.end()) {
            m->mothurOut("[ERROR]: category " + catName + " is not in your design file. Please correct.\n"); m->control_pressed = true;
        }else {
            for (map<string, int>::iterator it = totalCategories[it2->second].begin(); it != totalCategories[it2->second].end(); it++) {
                values.push_back(it->first);
            }
        }
        
        return values;
    }
	catch(exception& e) {
		m->errorOut(e, "DesignMap", "getValues");
		exit(1);
	}
}

/************************************************************/
////groupName, category returns value. example F000132, sex -> male
string DesignMap::get(string groupName, string categoryName) {
    try {
        string value = "not found";
        map<string, int>::iterator it = indexCategoryMap.find(categoryName);
        if (it == indexCategoryMap.end()) {
            m->mothurOut("[ERROR]: category " + categoryName + " is not in your design file. Please correct.\n"); m->control_pressed = true;
        }else {
            map<string, int>::iterator it2 = indexNameMap.find(groupName);
            if (it2 == indexNameMap.end()) {
                m->mothurOut("[ERROR]: group " + groupName + " is not in your design file. Please correct.\n"); m->control_pressed = true;
            }else {
                return designMap[it2->second][it->second];
            }
        }
        return value;
    }
	catch(exception& e) {
		m->errorOut(e, "DesignMap", "get");
		exit(1);
	}
}
/************************************************************/
//add group, assumes order is correct
int DesignMap::push_back(string group, vector<string> values) {
    try {
        map<string, int>::iterator it = indexNameMap.find(group);
        if (it == indexNameMap.end()) {
            if (values.size() != getNumCategories()) {  m->mothurOut("[ERROR]: Your design file has a " + toString(getNumCategories()) + " categories and " + group + " has " + toString(values.size()) + ", please correct."); m->mothurOutEndLine(); m->control_pressed = true;  return 0; }
            
            for (int i = 0; i < values.size(); i++) {
                //do we have this value for this category already
                map<string, int>::iterator it = totalCategories[i].find(values[i]);
                if (it == totalCategories[i].end()) { totalCategories[i][values[i]] = 1; }
                else {  totalCategories[i][values[i]]++; }
            }
            int count = indexNameMap.size();
            indexNameMap[group] = count;
            designMap.push_back(values);
        }else {
            m->mothurOut("[ERROR]: Your design file contains more than 1 group named " + group + ", group names must be unique. Please correct."); m->mothurOutEndLine(); m->control_pressed = true;
        }
        
        return 0;
    }
	catch(exception& e) {
		m->errorOut(e, "DesignMap", "push_back");
		exit(1);
	}
}
/************************************************************/
//set values for group, does not need to set all values. assumes group is in table already
int DesignMap::set(string group, map<string, string> values) {
    try {
        map<string, int>::iterator it = indexNameMap.find(group);
        if (it != indexNameMap.end()) {
            for (map<string, string>::iterator it2 = values.begin(); it2 != values.end(); it2++) {
                
                map<string, int>::iterator it3 = indexCategoryMap.find(it2->first); //do we have this category
                if (it3 == indexCategoryMap.end()) {
                     m->mothurOut("[ERROR]: Your design file does not contain a category called " + it2->first + ". Please correct."); m->mothurOutEndLine(); m->control_pressed = true;
                }else {
                    string oldCategory = designMap[it->second][it3->second];
                    //adjust totals for old category
                    int oldCount = totalCategories[it3->second][oldCategory];
                    if (oldCount == 1) { totalCategories[it3->second].erase(oldCategory); }
                    else {  totalCategories[it3->second][oldCategory]--; }
                    
                    designMap[it->second][it3->second] = it2->second; //reset value
                    
                    //adjust totals for new category
                    map<string, int>::iterator it4 = totalCategories[it3->second].find(it2->second);
                    if (it4 == totalCategories[it3->second].end()) { totalCategories[it3->second][it2->second] = 1; }
                    else {  totalCategories[it3->second][it2->second]++; }
                }
            }
        }else {
            m->mothurOut("[ERROR]: Your design file does not contain a group named " + group + ". Please correct."); m->mothurOutEndLine(); m->control_pressed = true;
        }
        
        return 0;
    }
	catch(exception& e) {
		m->errorOut(e, "DesignMap", "set");
		exit(1);
	}
}
/************************************************************/
//get number of groups belonging to a category or set of categories, with value or a set of values. Must have all categories and values. Example:
//  map<treatment - > early, late>, <sex -> male> would return 1. Only one group is male and from early or late.
int DesignMap::getNumUnique(map<string, vector<string> > selected) {
    try {
        int num = 0;
        
        //get indexes of categories
        vector<int> indexes;
        for (map<string, vector<string> >::iterator it = selected.begin(); it != selected.end(); it++) {
            map<string, int>::iterator it2 = indexCategoryMap.find(it->first);
            if (it2 == indexCategoryMap.end()) {
                m->mothurOut("[ERROR]: Your design file does not contain a category named " + it->first + ". Please correct."); m->mothurOutEndLine(); m->control_pressed = true; return 0;
            }else { indexes.push_back(it2->second); }
        }
        
        for (int i = 0; i < designMap.size(); i++) {
            bool hasAll = true; //innocent til proven guilty
            int count = 0;
            for (map<string, vector<string> >::iterator it = selected.begin(); it != selected.end(); it++) { //loop through each
                //check category to see if this group meets the requirements
                if (!m->inUsersGroups(designMap[i][indexes[count]], it->second)) { hasAll = false; it = selected.end(); }
                count++;
            }
            if (hasAll) { num++; }
        }
        
        return num;
    }
	catch(exception& e) {
		m->errorOut(e, "DesignMap", "getNumUnique");
		exit(1);
	}
}
/************************************************************/
//get number of groups belonging to a category or set of categories, with value or a set of values. Must have at least one categories and values. Example:
//  map<treatment - > early, late>, <sex -> male> would return 3. All three group have are either male or from early or late.
int DesignMap::getNumShared(map<string, vector<string> > selected) {
    try {
        int num = 0;
        
        //get indexes of categories
        vector<int> indexes;
        for (map<string, vector<string> >::iterator it = selected.begin(); it != selected.end(); it++) {
            map<string, int>::iterator it2 = indexCategoryMap.find(it->first);
            if (it2 == indexCategoryMap.end()) {
                m->mothurOut("[ERROR]: Your design file does not contain a category named " + it->first + ". Please correct."); m->mothurOutEndLine(); m->control_pressed = true; return 0;
            }else { indexes.push_back(it2->second); }
        }
        
        for (int i = 0; i < designMap.size(); i++) {
            bool hasAny = false; //guilty til proven innocent
            int count = 0;
            for (map<string, vector<string> >::iterator it = selected.begin(); it != selected.end(); it++) { //loop through each
                //check category to see if this group meets the requirements
                if (m->inUsersGroups(designMap[i][indexes[count]], it->second)) { hasAny = true; it = selected.end(); }
                count++;
            }
            if (hasAny) { num++; }
        }

        
        return num;
    }
	catch(exception& e) {
		m->errorOut(e, "DesignMap", "getNumShared");
		exit(1);
	}
}

/************************************************************/
//get names of groups belonging to a category or set of categories, with value or a set of values. Must have all categories and values. Example:
//  map<treatment - > early, late>, <sex -> male> would return F000132. F000132 is the only group which is male and from early or late.
vector<string> DesignMap::getNamesUnique(map<string, vector<string> > selected) {
    try {
        vector<string> names;
        
        //get indexes of categories
        vector<int> indexes;
        for (map<string, vector<string> >::iterator it = selected.begin(); it != selected.end(); it++) {
            map<string, int>::iterator it2 = indexCategoryMap.find(it->first);
            if (it2 == indexCategoryMap.end()) {
                m->mothurOut("[ERROR]: Your design file does not contain a category named " + it->first + ". Please correct."); m->mothurOutEndLine(); m->control_pressed = true; return names;
            }else { indexes.push_back(it2->second); }
        }
        
        //map int to name
        map<int, string> reverse;
        for (map<string, int>::iterator it = indexNameMap.begin(); it != indexNameMap.end(); it++) {
            reverse[it->second] = it->first;
        }
        
        for (int i = 0; i < designMap.size(); i++) {
            bool hasAll = true; //innocent til proven guilty
            int count = 0;
            for (map<string, vector<string> >::iterator it = selected.begin(); it != selected.end(); it++) { //loop through each
                //check category to see if this group meets the requirements
                if (!m->inUsersGroups(designMap[i][indexes[count]], it->second)) { hasAll = false; it = selected.end(); }
                count++;
            }
            if (hasAll) {
                map<int, string>::iterator it = reverse.find(i);
                if (it == reverse.end()) {
                    m->mothurOut("[ERROR]: should never get here, oops. Please correct."); m->mothurOutEndLine(); m->control_pressed = true; return names;
                }else { names.push_back(it->second); }
            }
        }

        
        return names;
    }
	catch(exception& e) {
		m->errorOut(e, "DesignMap", "getNamesUnique");
		exit(1);
	}
}
/************************************************************/
//get names of groups belonging to a category or set of categories, with value or a set of values. Must have all categories and values. Example:
//  map<treatment - > early, late>, <sex -> male> would return F000132. F000132 is the only group which is male and from early or late.
vector<string> DesignMap::getNamesShared(map<string, vector<string> > selected) {
    try {
        vector<string> names;
        
        //get indexes of categories
        vector<int> indexes;
        for (map<string, vector<string> >::iterator it = selected.begin(); it != selected.end(); it++) {
            map<string, int>::iterator it2 = indexCategoryMap.find(it->first);
            if (it2 == indexCategoryMap.end()) {
                m->mothurOut("[ERROR]: Your design file does not contain a category named " + it->first + ". Please correct."); m->mothurOutEndLine(); m->control_pressed = true; return names;
            }else { indexes.push_back(it2->second); }
        }
        
        //map int to name
        map<int, string> reverse;
        for (map<string, int>::iterator it = indexNameMap.begin(); it != indexNameMap.end(); it++) {
            reverse[it->second] = it->first;
        }
        
        for (int i = 0; i < designMap.size(); i++) {
            bool hasAny = false; 
            int count = 0;
            for (map<string, vector<string> >::iterator it = selected.begin(); it != selected.end(); it++) { //loop through each
                //check category to see if this group meets the requirements
                if (m->inUsersGroups(designMap[i][indexes[count]], it->second)) { hasAny = true; it = selected.end(); }
                count++;
            }
            if (hasAny) {
                map<int, string>::iterator it = reverse.find(i);
                if (it == reverse.end()) {
                    m->mothurOut("[ERROR]: should never get here, oops. Please correct."); m->mothurOutEndLine(); m->control_pressed = true; return names;
                }else { names.push_back(it->second); }
            }
        }
        
        
        return names;
    }
	catch(exception& e) {
		m->errorOut(e, "DesignMap", "getNamesShared");
		exit(1);
	}
}

/************************************************************/
//get names of groups belonging to a category or set of categories, with value or a set of values. Must have at least one categories and values. Example:
//  map<treatment - > early, late>, <sex -> male> would return F000132, F000142, F000138. All three group have are either male or from early or late.

vector<string> DesignMap::getNames(string category, string value) {
    try {
        vector<string> names;
        
        map<string, int>::iterator it = indexCategoryMap.find(category);
        if (it == indexCategoryMap.end()) {
            m->mothurOut("[ERROR]: category " + category + " is not in your design file. Please correct.\n"); m->control_pressed = true;
        }else {
            int column = it->second;
            
            //map int to name
            map<int, string> reverse;
            for (map<string, int>::iterator it2 = indexNameMap.begin(); it2 != indexNameMap.end(); it2++) {
                reverse[it2->second] = it2->first;
            }
            
            for (int i = 0; i < designMap.size(); i++) {
                if (designMap[i][column] == value) {
                    map<int, string>::iterator it2 = reverse.find(i);
                    if (it2 == reverse.end()) {
                        m->mothurOut("[ERROR]: should never get here, oops. Please correct."); m->mothurOutEndLine(); m->control_pressed = true; return names;
                    }else { names.push_back(it2->second); }
                }
            }
        }
        return names;
        
    }
	catch(exception& e) {
		m->errorOut(e, "DesignMap", "getNames");
		exit(1);
	}
}

/************************************************************/
int DesignMap::print(ofstream& out) {
    try {
       
		out << "group\t";
        for (int i = 0; i < namesOfCategories.size(); i++) { out << namesOfCategories[i] << '\t'; }
        out << endl;
        
        map<int, string> reverse; //use this to preserve order
        for (map<string, int>::iterator it = indexNameMap.begin(); it !=indexNameMap.end(); it++) { reverse[it->second] = it->first;  }
        
        for (int i = 0; i < designMap.size(); i++) {
            map<int, string>::iterator itR = reverse.find(i);
            
            if (itR != reverse.end()) { //will equal end if seqs were removed because remove just removes from indexNameMap
                out << itR->second  << '\t';
                
                for (int j = 0; j < namesOfCategories.size(); j++) {
                    out << designMap[i][j] << '\t';
                }
                out << endl;
            }
        }
        out.close();
        
        return 0;
    }
	catch(exception& e) {
		m->errorOut(e, "DesignMap", "print");
		exit(1);
	}
}
/************************************************************/
//print specific categories
int DesignMap::print(ofstream& out, vector<string> cats) {
    try {
        
		out << "group\t";
        for (int i = 0; i < namesOfCategories.size(); i++) { if (m->inUsersGroups(namesOfCategories[i], cats)) { out << namesOfCategories[i] << '\t'; } }
        out << endl;
        
        map<int, string> reverse; //use this to preserve order
        for (map<string, int>::iterator it = indexNameMap.begin(); it !=indexNameMap.end(); it++) { reverse[it->second] = it->first;  }
        
        for (int i = 0; i < designMap.size(); i++) {
            map<int, string>::iterator itR = reverse.find(i);
            
            if (itR != reverse.end()) { //will equal end if seqs were removed because remove just removes from indexNameMap
                out << itR->second  << '\t';
                
                for (int j = 0; j < namesOfCategories.size(); j++) {
                    if (m->inUsersGroups(namesOfCategories[i], cats)) {
                        out << designMap[i][j] << '\t';
                    }
                }
                out << endl;
            }
        }
        out.close();
        
        return 0;
    }
	catch(exception& e) {
		m->errorOut(e, "DesignMap", "print");
		exit(1);
	}
}

/************************************************************/

