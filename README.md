# myPhyloDB       ![myPhyloDB Logo](media/images/myPhyloDB_Logo.png)
<hr>
* myPhyloDB is a user-friendly web-interface for accessing and analyzing taxonomic data from multiple projects and/or sequencing runs.
* The goal of myPhyloDB is to allow for easy comparisons and statistical analysis of microbial (i.e., fungi or bacteria) taxonomic abundance across projects, environments, and management scenarios.
* myPhyloDB accepts Mothur pre-analyzed data from any of its supported sequencing platforms. In addition, raw 454 and Illumina/MiSeq files may be uploaded and analyzed using a user customizable Mothur pipeline.
* Please visit our [website] (http://www.myphylodb.org) for additional information on how to run myPhyloDB.

<hr>

# Installers
<hr>
* Installers for the current Linux and Windows stand-alone executables can be found [here] (http://www.ars.usda.gov/services/software/download.htm?softwareid=472).
* Executables include all dependencies, allowing myPhyloDB to be run on computers without Python, Mothur or R.

<hr>

# Dependencies
<hr>
Dependencies (and the source code provided here) are only required for development purposes.  Most users should use the installers provided above.

####**Linux systems**

***Python 2.7.8***

* **Modules:** CherryPy 3.6.0, Django 1.6.5, Django-extensions 1.4.4, Django-registration-redux 1.2, PypeR 1.1.2, numPy 1.8.1, pandas 0.17.1, sciPy 0.14.0, simplejson 3.5.2

***R 3.1.2***

* **Packages:** vegan, multcomp, lsmeans, spls, mixOmics, DESeq2
* **Source:** [R for Linux] (http://lib.stat.cmu.edu/R/CRAN/)
* **Installation folder:** 'R/R-Linux'

***Mothur 1.35.1***

* **Source:** [Mothur.cen_64.zip] (https://github.com/mothur/mothur/releases/tag/v1.35.1)
* **Installation folder:** 'mothur/mothur-linux'
       
***Mothur reference files***

* **Source:** [silva.seed_v119.align] (http://www.mothur.org/wiki/Taxonomy_outline)
* **Installation folder:** 'mothur/reference/align'


* **Example:** [gg_13_5_99.pds.tax]  (http://www.mothur.org/wiki/Taxonomy_outline)
* **Installation folder:** 'mothur/reference/taxonomy'


* **Example:** [gg_13_5_99.fasta]  (http://www.mothur.org/wiki/Taxonomy_outline)
* **Installation folder:** 'mothur/reference/template'

####**Windows systems**

***Python 2.7.8***

* **Modules:** CherryPy 3.6.0, Django 1.6.5, Django-extensions 1.4.4, Django-registration-redux 1.2, PypeR 1.1.2, numPy 1.8.1, pandas 0.14.0, sciPy 0.14.0, simplejson 3.5.2

***R 3.1.2***

* **Packages:** vegan, multcomp, lsmeans, spls, mixOmics, DESeq2
* **Source:** [R Portable] (http://sourceforge.net/projects/rportable/)
* **Installation folder:** 'R/R-Portable'

***Mothur 1.35.1***

* **Source:** [Mothur.win_64.zip] (https://github.com/mothur/mothur/releases/tag/v1.35.1)
* **Installation folder:** 'mothur/mothur-win'
       
***Mothur reference files***

* **Source:** [silva.seed_v119.align] (http://www.mothur.org/wiki/Taxonomy_outline)
* **Installation folder:** 'mothur/reference/align'


* **Example:** [gg_13_5_99.pds.tax]  (http://www.mothur.org/wiki/Taxonomy_outline)
* **Installation folder:** 'mothur/reference/taxonomy'


* **Example:** [gg_13_5_99.fasta]  (http://www.mothur.org/wiki/Taxonomy_outline)
* **Installation folder:** 'mothur/reference/template'

<hr>
