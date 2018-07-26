from database.views import pathTaxaJSON
from django.core.handlers.wsgi import WSGIRequest
from django.http.request import QueryDict
import json

# test kegg not test getOtuFromKoList ??? or both

# test kegg path with a level 3-4 id (ideally under 10 seconds to run this test)
def testQuickKegg():
    testPass = False
    myDict = QueryDict('key=6067cf4d7d814598a2c1b20f8372e591')
    myKeggRequest = WSGIRequest({
        'REQUEST_METHOD': 'GET',
        'wsgi.input': "heh?",
        'HTTP_X_REQUESTED_WITH': "XMLHttpRequest",
    })
    myKeggRequest.GET = myDict
    try:
        ret = pathTaxaJSON(myKeggRequest)
        retDat = json.loads(ret.content)
        if len(retDat['data']) != 0:
            testPass = True
    except Exception as errd:
        print "Error with quickKegg: ", errd
    # "6067cf4d7d814598a2c1b20f8372e591"
    #  get easy request copy from actual page, copy it here, run it
    #  Check response time and maybe number of outputs (both are ranges for variance sake)
    return testPass
