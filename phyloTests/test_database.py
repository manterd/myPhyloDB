# database specific testing sets

import test_keggorthology as tk
import time

def TestKeggOrthology():
    print "Testing Kegg Orthology"
    passedAll = True
    passedTests = ""
    failedTests = ""
    testStart = time.time()
    # run tests for each section of anova

    '''
    actual tests here
    '''

    testName = "QuickKegg "
    if not tk.testQuickKegg():
        failedTests += testName
        passedAll = False
    else:
        passedTests += testName
    '''
    test small functions with under 5 second runtimes (static known test parameters for consistent results)
    '''

    print "Database testing finished. Time elapsed ", time.time()-testStart, " seconds"

    if not passedAll:
        print "Failed Kegg Orthology tests:"
        print failedTests + "\n"
    return passedAll