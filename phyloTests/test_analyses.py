# analysis specific testing sets

import test_anova as ta
import test_corrplot as tcp
import time


def TestAnova():
    print "Testing Anova"
    passedAll = True
    passedTests = ""
    failedTests = ""
    # run tests for each section of anova

    '''
    actual tests here
    '''
    testStart = time.time()
    # make dummy analysis (very simple anova)
    testAnalysis = ta.createTestAnova()

    testName = "BadValidate "
    if not ta.testBadValidate():
        failedTests += testName
        passedAll = False
    else:
        passedTests += testName


    # test works but takes a while
    testName = "GoodValidate "
    if not ta.testGoodValidate(testAnalysis):
        failedTests += testName
        passedAll = False
    else:
        passedTests += testName

    testName = "BadQuery "
    if not ta.testBadQuery():
        failedTests += testName
        passedAll = False
    else:
        passedTests += testName

    testName = "GoodQuery "
    if not ta.testGoodQuery(testAnalysis):
        failedTests += testName
        passedAll = False
    else:
        passedTests += testName

    testName = "StatsOutput "
    if not ta.testStatsOutput(testAnalysis):
        failedTests += testName
        passedAll = False
    else:
        passedTests += testName

    testName = "GraphOutput "
    if not ta.testGraphOutput(testAnalysis):
        failedTests += testName
        passedAll = False
    else:
        passedTests += testName

    print "Anova testing finished. Time elapsed ", time.time()-testStart, " seconds"

    # verify all anova passed
    if not passedAll:
        print "\tPassed Anova tests:"
        print "\t" + passedTests
        print "\tFailed Anova tests:"
        print "\t" + failedTests
    else:
        print "\tPassed all Anova tests"
    return passedAll


def TestCorrPlot():
    print "Testing CorrPlot"
    passedAll = True
    passedTests = ""
    failedTests = ""
    # run tests for each section of anova

    '''
    actual tests here
    '''
    testStart = time.time()
    # make dummy analysis (very simple corrplot)
    testAnalysis = tcp.createTestCorr()

    testName = "BadValidate "
    if not tcp.testBadValidate():
        failedTests += testName
        passedAll = False
    else:
        passedTests += testName

    # test works but takes a while
    testName = "GoodValidate "
    if not tcp.testGoodValidate(testAnalysis):
        failedTests += testName
        passedAll = False
    else:
        passedTests += testName

    testName = "BadQuery "
    if not tcp.testBadQuery():
        failedTests += testName
        passedAll = False
    else:
        passedTests += testName

    testName = "GoodQuery "
    if not tcp.testGoodQuery(testAnalysis):
        failedTests += testName
        passedAll = False
    else:
        passedTests += testName

    testName = "StatsOutput "
    if not tcp.testStatsOutput(testAnalysis):
        failedTests += testName
        passedAll = False
    else:
        passedTests += testName

    testName = "GraphOutput "
    if not tcp.testGraphOutput(testAnalysis):
        failedTests += testName
        passedAll = False
    else:
        passedTests += testName

    print "Corr testing finished. Time elapsed ", time.time() - testStart, " seconds"

    # verify all anova passed
    if not passedAll:
        print "\tPassed CorrPlot tests:"
        print "\t" + passedTests
        print "\tFailed CorrPlot tests:"
        print "\t" + failedTests
    else:
        print "\tPassed all CorrPlot tests"
    return passedAll