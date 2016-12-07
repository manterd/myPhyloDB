function checkDupes(list) {
            var bools = [];
            var firsts = [];
            for (var i=0; i<list.length; i++) {
                var sub = list[i].split(':');
                var found = firsts.indexOf(sub[0]);
                if (found==-1) {
                    firsts.push(sub[0]);
                    bools.push(false);
                } else {
                    bools[found] = true;
                }
            }
            return (bools.indexOf(false)==-1);
        }

function clearMeta() {
    if (dataType == 0) {
        $("#tree_metaCat").dynatree("getRoot").visit(function (node) {
            node.select(false);
        });
        $("#tree_metaQuant").dynatree("getRoot").visit(function (node) {
            node.select(false);
        });
    }
    else if (dataType == 1) {
        $("#tree_metaCat").dynatree("getRoot").visit(function (node) {
            node.select(false);
        });
    }
    else if (dataType == 2) {
        $("#tree_metaQuant").dynatree("getRoot").visit(function (node) {
            node.select(false);
        });
    }
    $("#run").css("background-color", "lightgray");
}

function clearTaxa() {
    if (treeType == 1) {
        $("#tree_taxa").dynatree("getRoot").visit(function (node) {
            node.select(false);
        });
    } else if (treeType == 2) {
        $("#tree_kegg").dynatree("getRoot").visit(function (node) {
            node.select(false);
        });
    } else if (treeType == 3) {
        $("#tree_nz").dynatree("getRoot").visit(function (node) {
            node.select(false);
        });
    }

    $("#run").css("background-color", "lightgray");
}

function makeGUID() {
            return 'xxxxxxxx-xxxx-4xxx-yxxx-xxxxxxxxxxxx'.replace(/[xy]/g, function(c) {
                var r = Math.random() * 16 | 0, v = c == 'x' ? r : (r & 0x3 | 0x8);
                return v.toString(16);
            })
        }

function startTimer() {
    refresh = setInterval(function() {
        updateStatus(RID);
    }, 1000);
    return refresh;
}
