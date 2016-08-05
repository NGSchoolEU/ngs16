// (C) Wolfgang Huber 2010-2011

// Script parameters - these are set up by R in the function 'writeReport' when copying the 
//   template for this script from arrayQualityMetrics/inst/scripts into the report.

var highlightInitial = [ false, false, false, false, false, false, false, false, true, false, false, false ];
var arrayMetadata    = [ [ "1", "GSE3325_1", "/Users/marzec01/Desktop/others/workshops/Lechu_Slowacja_bix_workshop_08_2016/workshop/data/CEL_files/GSE3325/GSM74875.CEL", "Normal", "07/01/05 12:12:28" ], [ "2", "GSE3325_2", "/Users/marzec01/Desktop/others/workshops/Lechu_Slowacja_bix_workshop_08_2016/workshop/data/CEL_files/GSE3325/GSM74876.CEL", "Normal", "07/01/05 12:36:46" ], [ "3", "GSE3325_3", "/Users/marzec01/Desktop/others/workshops/Lechu_Slowacja_bix_workshop_08_2016/workshop/data/CEL_files/GSE3325/GSM74877.CEL", "Normal", "07/01/05 11:47:58" ], [ "4", "GSE3325_4", "/Users/marzec01/Desktop/others/workshops/Lechu_Slowacja_bix_workshop_08_2016/workshop/data/CEL_files/GSE3325/GSM74879.CEL", "Normal", "07/01/05 12:24:37" ], [ "5", "GSE3325_5", "/Users/marzec01/Desktop/others/workshops/Lechu_Slowacja_bix_workshop_08_2016/workshop/data/CEL_files/GSE3325/GSM74880.CEL", "Normal", "07/01/05 12:48:56" ], [ "6", "GSE3325_6", "/Users/marzec01/Desktop/others/workshops/Lechu_Slowacja_bix_workshop_08_2016/workshop/data/CEL_files/GSE3325/GSM74881.CEL", "Tumour", "07/07/05 12:25:02" ], [ "7", "GSE3325_7", "/Users/marzec01/Desktop/others/workshops/Lechu_Slowacja_bix_workshop_08_2016/workshop/data/CEL_files/GSE3325/GSM74882.CEL", "Tumour", "07/07/05 12:12:53" ], [ "8", "GSE3325_8", "/Users/marzec01/Desktop/others/workshops/Lechu_Slowacja_bix_workshop_08_2016/workshop/data/CEL_files/GSE3325/GSM74883.CEL", "Tumour", "07/07/05 13:13:43" ], [ "9", "GSE3325_9", "/Users/marzec01/Desktop/others/workshops/Lechu_Slowacja_bix_workshop_08_2016/workshop/data/CEL_files/GSE3325/GSM74884.CEL", "Tumour", "07/07/05 11:48:25" ], [ "10", "GSE3325_10", "/Users/marzec01/Desktop/others/workshops/Lechu_Slowacja_bix_workshop_08_2016/workshop/data/CEL_files/GSE3325/GSM74885.CEL", "Tumour", "07/07/05 12:00:41" ], [ "11", "GSE3325_11", "/Users/marzec01/Desktop/others/workshops/Lechu_Slowacja_bix_workshop_08_2016/workshop/data/CEL_files/GSE3325/GSM74886.CEL", "Tumour", "07/07/05 12:49:21" ], [ "12", "GSE3325_12", "/Users/marzec01/Desktop/others/workshops/Lechu_Slowacja_bix_workshop_08_2016/workshop/data/CEL_files/GSE3325/GSM74887.CEL", "Tumour", "07/07/05 12:37:14" ] ];
var svgObjectNames   = [ "pca", "dens", "dig" ];

var cssText = ["stroke-width:1; stroke-opacity:0.4",
               "stroke-width:3; stroke-opacity:1" ];

// Global variables - these are set up below by 'reportinit'
var tables;             // array of all the associated ('tooltips') tables on the page
var checkboxes;         // the checkboxes
var ssrules;


function reportinit() 
{
 
    var a, i, status;

    /*--------find checkboxes and set them to start values------*/
    checkboxes = document.getElementsByName("ReportObjectCheckBoxes");
    if(checkboxes.length != highlightInitial.length)
	throw new Error("checkboxes.length=" + checkboxes.length + "  !=  "
                        + " highlightInitial.length="+ highlightInitial.length);
    
    /*--------find associated tables and cache their locations------*/
    tables = new Array(svgObjectNames.length);
    for(i=0; i<tables.length; i++) 
    {
        tables[i] = safeGetElementById("Tab:"+svgObjectNames[i]);
    }

    /*------- style sheet rules ---------*/
    var ss = document.styleSheets[0];
    ssrules = ss.cssRules ? ss.cssRules : ss.rules; 

    /*------- checkboxes[a] is (expected to be) of class HTMLInputElement ---*/
    for(a=0; a<checkboxes.length; a++)
    {
	checkboxes[a].checked = highlightInitial[a];
        status = checkboxes[a].checked; 
        setReportObj(a+1, status, false);
    }

}


function safeGetElementById(id)
{
    res = document.getElementById(id);
    if(res == null)
        throw new Error("Id '"+ id + "' not found.");
    return(res)
}

/*------------------------------------------------------------
   Highlighting of Report Objects 
 ---------------------------------------------------------------*/
function setReportObj(reportObjId, status, doTable)
{
    var i, j, plotObjIds, selector;

    if(doTable) {
	for(i=0; i<svgObjectNames.length; i++) {
	    showTipTable(i, reportObjId);
	} 
    }

    /* This works in Chrome 10, ssrules will be null; we use getElementsByClassName and loop over them */
    if(ssrules == null) {
	elements = document.getElementsByClassName("aqm" + reportObjId); 
	for(i=0; i<elements.length; i++) {
	    elements[i].style.cssText = cssText[0+status];
	}
    } else {
    /* This works in Firefox 4 */
	var success = false;
	i = 0; 
	/* Some of this looping could already be cached in reportInit() */
	while( (!success) & (i < ssrules.length) ) {
	    selector = ssrules[i].selectorText;  // The selector 
            if (!selector) 
		continue; // Skip @import and other nonstyle rules
            if (selector == (".aqm" + reportObjId)) {
		success = true; 
		ssrules[i].style.cssText = cssText[0+status];
	    } else {
		i++;
	    }
	}
    }

}

/*------------------------------------------------------------
   Display of the Metadata Table
  ------------------------------------------------------------*/
function showTipTable(tableIndex, reportObjId)
{
    var rows = tables[tableIndex].rows;
    var a = reportObjId - 1;

    if(rows.length != arrayMetadata[a].length)
	throw new Error("rows.length=" + rows.length+"  !=  arrayMetadata[array].length=" + arrayMetadata[a].length);

    for(i=0; i<rows.length; i++) 
 	rows[i].cells[1].innerHTML = arrayMetadata[a][i];
}

function hideTipTable(tableIndex)
{
    var rows = tables[tableIndex].rows;

    for(i=0; i<rows.length; i++) 
 	rows[i].cells[1].innerHTML = "";
}


/*------------------------------------------------------------
  From module 'name' (e.g. 'density'), find numeric index in the 
  'svgObjectNames' array.
  ------------------------------------------------------------*/
function getIndexFromName(name) 
{
    var i;
    for(i=0; i<svgObjectNames.length; i++)
        if(svgObjectNames[i] == name)
	    return i;

    throw new Error("Did not find '" + name + "'.");
}


/*------------------------------------------------------------
  SVG plot object callbacks
  ------------------------------------------------------------*/
function plotObjRespond(what, reportObjId, name)
{

    var a, i, status;

    switch(what) {
    case "show":
	i = getIndexFromName(name);
	showTipTable(i, reportObjId);
	break;
    case "hide":
	i = getIndexFromName(name);
	hideTipTable(i);
	break;
    case "click":
        a = reportObjId - 1;
	status = !checkboxes[a].checked;
	checkboxes[a].checked = status;
	setReportObj(reportObjId, status, true);
	break;
    default:
	throw new Error("Invalid 'what': "+what)
    }
}

/*------------------------------------------------------------
  checkboxes 'onchange' event
------------------------------------------------------------*/
function checkboxEvent(reportObjId)
{
    var a = reportObjId - 1;
    var status = checkboxes[a].checked;
    setReportObj(reportObjId, status, true);
}


/*------------------------------------------------------------
  toggle visibility
------------------------------------------------------------*/
function toggle(id){
  var head = safeGetElementById(id + "-h");
  var body = safeGetElementById(id + "-b");
  var hdtxt = head.innerHTML;
  var dsp;
  switch(body.style.display){
    case 'none':
      dsp = 'block';
      hdtxt = '-' + hdtxt.substr(1);
      break;
    case 'block':
      dsp = 'none';
      hdtxt = '+' + hdtxt.substr(1);
      break;
  }  
  body.style.display = dsp;
  head.innerHTML = hdtxt;
}
