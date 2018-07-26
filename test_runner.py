'''

This file should have access to all unit tests and should have most (if not all) tests run at the start
and end of every coding session, feature test, etc

Unit tests should have as few dependencies as possible (test queue code functionality without
booting into cherrypy each time for example). The idea is for each test to be as specific as it can

'''

import os
os.environ['DJANGO_SETTINGS_MODULE'] = 'myPhyloDB.settings'
import django
django.setup()

import phyloTests.test_analyses as ta
import phyloTests.test_database as td


def TestAnalyses():
    passedAll = True
    failedTests = ""
    # run tests for each analysis
    if not ta.TestAnova():
        failedTests += "Anova "
        passedAll = False
    if not ta.TestCorrPlot():
        failedTests += "CorrPlot "
        passedAll = False
    # verify all analyses passed
    if not passedAll:
        print "Failed Analyses tests:"
        print failedTests + "\n"
    return passedAll


def TestDatabase():
    passedAll = True
    failedTests = ""
    # run tests for each analysis
    if not td.TestKeggOrthology():
        failedTests += "KeggOrthology "
        passedAll = False
    # verify all analyses passed
    if not passedAll:
        print "Failed Database tests:"
        print failedTests + "\n"
    return passedAll


def runAllTests():
    passedAll = True
    failedTests = ""
    # do actual tests, set passedAll to False if any fail
    if not TestAnalyses():
        failedTests += "Analyses "
        passedAll = False
    if not TestDatabase():
        failedTests += "Database "
        passedAll = False
    # verify all tests passed
    if passedAll:
        print "Passed all tests!"
    else:
        print "Failed the following tests:"
        print failedTests
    return

# "main" method, does exactly as stated
runAllTests()
