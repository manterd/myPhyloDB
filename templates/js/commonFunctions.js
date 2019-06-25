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
    $("#tree_metaCat").dynatree("getRoot").visit(function (node) {
        node.select(false);
    });
    $("#tree_metaQuant").dynatree("getRoot").visit(function (node) {
        node.select(false);
    });
    $("#run").css("background-color", "lightgray");
}

function clearTaxa() {
    $("#tree_taxa").dynatree("getRoot").visit(function (node) {
        node.select(false);
    });
    $("#tree_kegg").dynatree("getRoot").visit(function (node) {
        node.select(false);
    });
    $("#tree_nz").dynatree("getRoot").visit(function (node) {
        node.select(false);
    });

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
        updateStatus();
    }, 1000);
    return refresh;
}

function setNorm(arg) {
    normNorm = arg
}

function menu () {
    selFiles = "{{ selFiles }}";
    normFiles = "{{ normFiles }}";
    dataID = "{{ dataID }}";
    if (normNorm == true) normFiles = "True";
    normNorm = false;

    if ((selFiles == "True") && (normFiles == "True"))  {
        $("#link0").show();
        $("#link1").show();
        $("#link2").show();
        $("#link3").show();
        $("#link4").show();
        $("#link5").show();
        $("#link6").show();
    }

    if ((selFiles == "False") && (normFiles == "True"))  {
        $("#link0").show();
        $("#link1").show();
        $("#link2").show();
        $("#link3").show();
        $("#link4").show();
        $("#link5").hide();
        $("#link6").show();
    }

    if ((selFiles == "True") && (normFiles == "False"))  {
        $("#link0").hide();
        $("#link1").hide();
        $("#link2").hide();
        $("#link3").hide();
        $("#link4").hide();
        $("#link5").show();
        $("#link6").hide();
    }
    if ((selFiles == "False") && (normFiles == "False"))  {
        $("#link0").hide();
        $("#link1").hide();
        $("#link2").hide();
        $("#link3").hide();
        $("#link4").hide();
        $("#link5").hide();
        $("#link6").hide();
    }
}

window.onload = function() {
    try {
        var ref = document.referrer.split("/");
        if (ref.length>2){
            var prevPage = ref[ref.length-2];
            if (prevPage=="history"){
                RID = localStorage.getItem("prevRID");  // this is neat in theory but presents an issue
                // when history has been loaded before, this loads old RID automatically on a new page
                if (RID != null && localStorage.getItem("fromHistory")=="true") {
                    localStorage.setItem("fromHistory", "false");
                    startTimer();
                }
            }

        }

    } catch(err) {
            // guess that didn't work, likely not here from history then
    }
};