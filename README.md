# myPhyloDB       ![myPhyloDB Logo](media/images/myPhyloDB_Logo.png)

* A user-friendly web-interface for accessing and analyzing taxonomic data from multiple projects and/or sequencing runs.
* The goal of myPhyloDB is to allow for easy comparisons and statistical analysis of microbial (i.e., fungi or bacteria) taxonomic abundance across projects, environments, and management scenarios.
* myPhyloDB accepts Mothur pre-analyzed data from any of its supported sequencing platforms. In addition, raw 454 and Illumina/MiSeq files may be uploaded and analyzed using a user customizable Mothur pipeline.


# Installers
* Installers (coming soon...) for Linux and Windows stand-alone executables can be found [here] (http://www.ars.usda.gov/services/software/software.htm).
* Executables include all dependencies, allowing myPhyloDB to be run on computers without Python, Mothur or R.

# Source code
**Dependencies:**
* CherryPy 3.6.0
* Django 1.6.5
* Django-extensions 1.4.4
* Django-registration-redux 1.2
* PypeR 1.1.2
* numPy 1.8.1
* pandas 0.14.0
* sciPy 0.14.0
* simplejson 3.5.2

**Additional Requirements (Linux):**
* R (installed in 'R/R-Linux')
* [Mothur] (http://www.mothur.org) (installed in 'mothur/mothur-linux')
* [Mothur Taxonomy Outline] (http://www.mothur.org/wiki/Taxonomy_outline) (installed in 'mothur/reference')

**Additional Requirements (Windows):**
* R-portable (installed in 'R/R-Portable')
* [Mothur] (http://www.mothur.org) (installed in 'mothur/mothur-win')
* [Mothur Taxonomy Outline] (http://www.mothur.org/wiki/Taxonomy_outline) (installed in 'mothur/reference')
