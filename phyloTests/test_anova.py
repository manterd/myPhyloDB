# specific function tests for anova analysis variations
# each positive test case (verify good output is good) also checks if time is within tolerance

import json
from django.http import HttpResponse

from functions.analysis import analysis


# dummy class intended for test case usage only
class Request:

    def __init__(self, body, user):
        self.body = body
        self.user = user

    def body(self):
        return self.body

    def user(self):
        return self.user


def testBadValidate():
    testPass = False
    # give an anova analysis object bad validation data
    # return true only if an error occurs (slightly backwards but needs to be done)
    RID = "TestingBadValidate"
    PID = 0
    stopList = []
    request = "Really, really bad json encoding. Like so bad that this string came through."
    request += " Think about it. It could happen"
    try:
        myAnalysis = analysis.Anova(request, RID, stopList, PID, debug=False)
        myAnalysis.validate()
    except Exception as e:
        testPass = True
    return testPass


def createTestAnova():
    RID = "TestingAnova"
    PID = 0
    stopList = {PID: "0"}
    myDict = {'selectAll': '2', 'transform': '0', 'taxa': '{}', 'filterData': 'no', 'metaValsCat': '{"usr_cat1" : "control","usr_cat1" : "control","usr_cat1" : "control","usr_cat1" : "control","usr_cat1" : "control","usr_cat1" : "control","usr_cat1" : "control","usr_cat1" : "control","usr_cat1" : "control","usr_cat1" : "control","usr_cat1" : "control","usr_cat1" : "control","usr_cat1" : "control","usr_cat1" : "control","usr_cat1" : "control","usr_cat1" : "treatment","usr_cat1" : "treatment","usr_cat1" : "treatment","usr_cat1" : "treatment","usr_cat1" : "treatment","usr_cat1" : "treatment","usr_cat1" : "treatment","usr_cat1" : "treatment","usr_cat1" : "treatment","usr_cat1" : "treatment","usr_cat1" : "treatment","usr_cat1" : "treatment","usr_cat1" : "treatment","usr_cat1" : "treatment"}', 'funcName': 'getCatUnivData', 'keggAll': '1', 'remUnclass': 'no', 'palette': 'Set1', 'metaValsQuant': '', 'remZeroes': 'no', 'map_taxa': 'no', 'gridVal_Y': 'None', 'gridVal_X': 'None', 'nz': '{}', 'nzAll': '1', 'treeType': '1', 'xVal': 'None', 'RID': 'TestingAnova', 'dataID': '92b3fce7a0f94cc3a44784c8feb53cab', 'kegg': '{}', 'filterMeth': '2', 'metaIDsQuant': '', 'DepVar': '0', 'metaIDsCat': '{"usr_cat1" : "dde8c08484bf4ed8b8929fa5ba6bb605","usr_cat1" : "f49b76eb74d74df7a96dc6139e7049f5","usr_cat1" : "d6916ad500394012a2aa3304fb9e8dc7","usr_cat1" : "41587496f5554d8897e0a22bfb23552d","usr_cat1" : "3abe15106eb743ad86ac550bce8b9750","usr_cat1" : "81b2cf573c7244e8afb77e5ad04f2256","usr_cat1" : "18cefacb8e8240aaa043b0d2a3a1a96f","usr_cat1" : "76915d8f32904da0b4c32273be9a57c6","usr_cat1" : "b9107dbe86bc4b28ab294d56716d0eca","usr_cat1" : "e25c98a2bc304c899d8449f2ef789f2d","usr_cat1" : "606cba425809479bafc33e745ff27142","usr_cat1" : "13d284e1ebc04f218c27410a2daadb41","usr_cat1" : "df7c87c700034b7b84dd9a590d2e92bd","usr_cat1" : "486ebf2154904a28ae461fd5ba28c1d0","usr_cat1" : "94f780071ed147e2bc80cdf84a37d32e","usr_cat1" : "98133e720df6447fb03cad8229a00275","usr_cat1" : "9f010781511f46ecafebceec7b1850d9","usr_cat1" : "9e2d953ae1f2414bb7540ec07bebeda3","usr_cat1" : "c8a939615de04987a44f7e853b633f12","usr_cat1" : "2238e24ce1f84a2db2d4ee39fdc7436e","usr_cat1" : "fb6362d4d8bf4687b0ded57983af7d13","usr_cat1" : "dc55f5a6a78b4f6f8a90544c9c81e2a6","usr_cat1" : "812f42a0f504472b838ba1cc84919d64","usr_cat1" : "a7fae1fca9634d88b2c4b9742f4f42d0","usr_cat1" : "b02f95fb584d4b11a708a95817246dd8","usr_cat1" : "e813485edf9c4db09cf50a685cf0ac33","usr_cat1" : "6907fa27097e4945af74bbb63193e23a","usr_cat1" : "0915d4670f8145cf8cc4bef98ee18209","usr_cat1" : "a9b2576e478d4fb6a1480c1a15a1b746"}', 'filterPer': '10', 'sig_only': False, 'colorVal': 'None', 'reqType': 'call', 'perZeroes': '50'}
    request = Request(json.dumps(myDict), "unit")
    myAnalysis = analysis.Anova(request, RID, stopList, PID, debug=False)
    return myAnalysis


def testGoodValidate(testAnalysis):
    testPass = False
    # give an anova analysis object good validation data
    # return true if no error occurs and function returns true
    try:
        if testAnalysis.validate() == 0:
            testPass = True
    except Exception as e:
        print "testGoodValidate error: ", e
    return testPass


def testBadQuery():
    testPass = False
    # give an anova analysis object bad query data
    # return true only if an error occurs
    RID = "TestingBadQuery"
    PID = 0
    stopList = []
    request = "Really, really bad json encoding. Like so bad that this string came through."
    request += " Think about it. It could happen"
    try:
        myAnalysis = analysis.Anova(request, RID, stopList, PID, debug=False)
        myAnalysis.query()
    except Exception as e:
        testPass = True
    return testPass


def testGoodQuery(testAnalysis):
    testPass = False
    # give an anova analysis object good query data
    # return true only if completely successful
    try:
        if testAnalysis.query() == 0:
            testPass = True
    except Exception as e:
        print "Query test error: ", e
    return testPass


def testStatsOutput(testAnalysis):
    testPass = False
    try:
        if testAnalysis.stats() == 0:
            testPass = True
    except Exception as e:
        print "Stats test error: ", e
    return testPass


def testGraphOutput(testAnalysis):
    testPass = False
    try:
        ret = testAnalysis.graph()
        if type(ret) == HttpResponse:
            testPass = True
    except Exception as e:
        print "Graph test error: ", e
    return testPass
