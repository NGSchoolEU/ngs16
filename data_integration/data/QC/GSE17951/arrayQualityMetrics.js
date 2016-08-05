// (C) Wolfgang Huber 2010-2011

// Script parameters - these are set up by R in the function 'writeReport' when copying the 
//   template for this script from arrayQualityMetrics/inst/scripts into the report.

var highlightInitial = [ false, true, false, false, true, false, true, false, false, false, false, false, false, false, false, false, false, false, false, true, false, false ];
var arrayMetadata    = [ [ "1", "GSE17951_1", "/Users/marzec01/Desktop/others/workshops/Lechu_Slowacja_bix_workshop_08_2016/workshop/data/CEL_files/GSE17951/GSM449150.CEL", "Tumour", "07/06/07 14:12:27" ], [ "2", "GSE17951_2", "/Users/marzec01/Desktop/others/workshops/Lechu_Slowacja_bix_workshop_08_2016/workshop/data/CEL_files/GSE17951/GSM449151.CEL", "Tumour", "08/10/07 22:41:57" ], [ "3", "GSE17951_3", "/Users/marzec01/Desktop/others/workshops/Lechu_Slowacja_bix_workshop_08_2016/workshop/data/CEL_files/GSE17951/GSM449160.CEL", "Tumour", "09/18/07 12:53:03" ], [ "4", "GSE17951_4", "/Users/marzec01/Desktop/others/workshops/Lechu_Slowacja_bix_workshop_08_2016/workshop/data/CEL_files/GSE17951/GSM449176.CEL", "Tumour", "10/16/07 14:14:13" ], [ "5", "GSE17951_5", "/Users/marzec01/Desktop/others/workshops/Lechu_Slowacja_bix_workshop_08_2016/workshop/data/CEL_files/GSE17951/GSM449189.CEL", "Tumour", "11/20/07 16:44:57" ], [ "6", "GSE17951_6", "/Users/marzec01/Desktop/others/workshops/Lechu_Slowacja_bix_workshop_08_2016/workshop/data/CEL_files/GSE17951/GSM449198.CEL", "Tumour", "12/07/07 14:35:38" ], [ "7", "GSE17951_7", "/Users/marzec01/Desktop/others/workshops/Lechu_Slowacja_bix_workshop_08_2016/workshop/data/CEL_files/GSE17951/GSM449227.CEL", "Tumour", "05/07/08 14:04:51" ], [ "8", "GSE17951_8", "/Users/marzec01/Desktop/others/workshops/Lechu_Slowacja_bix_workshop_08_2016/workshop/data/CEL_files/GSE17951/GSM449229.CEL", "Tumour", "05/07/08 14:16:06" ], [ "9", "GSE17951_9", "/Users/marzec01/Desktop/others/workshops/Lechu_Slowacja_bix_workshop_08_2016/workshop/data/CEL_files/GSE17951/GSM449230.CEL", "Tumour", "05/07/08 14:38:43" ], [ "10", "GSE17951_10", "/Users/marzec01/Desktop/others/workshops/Lechu_Slowacja_bix_workshop_08_2016/workshop/data/CEL_files/GSE17951/GSM449249.CEL", "Normal", "08/05/08 14:00:01" ], [ "11", "GSE17951_11", "/Users/marzec01/Desktop/others/workshops/Lechu_Slowacja_bix_workshop_08_2016/workshop/data/CEL_files/GSE17951/GSM449250.CEL", "Normal", "08/05/08 14:11:09" ], [ "12", "GSE17951_12", "/Users/marzec01/Desktop/others/workshops/Lechu_Slowacja_bix_workshop_08_2016/workshop/data/CEL_files/GSE17951/GSM449255.CEL", "Normal", "08/19/08 12:38:35" ], [ "13", "GSE17951_13", "/Users/marzec01/Desktop/others/workshops/Lechu_Slowacja_bix_workshop_08_2016/workshop/data/CEL_files/GSE17951/GSM449257.CEL", "Normal", "08/19/08 13:00:45" ], [ "14", "GSE17951_14", "/Users/marzec01/Desktop/others/workshops/Lechu_Slowacja_bix_workshop_08_2016/workshop/data/CEL_files/GSE17951/GSM449258.CEL", "Normal", "08/19/08 13:11:53" ], [ "15", "GSE17951_15", "/Users/marzec01/Desktop/others/workshops/Lechu_Slowacja_bix_workshop_08_2016/workshop/data/CEL_files/GSE17951/GSM449259.CEL", "Normal", "08/19/08 13:23:21" ], [ "16", "GSE17951_16", "/Users/marzec01/Desktop/others/workshops/Lechu_Slowacja_bix_workshop_08_2016/workshop/data/CEL_files/GSE17951/GSM449267.CEL", "Normal", "10/13/08 14:24:50" ], [ "17", "GSE17951_17", "/Users/marzec01/Desktop/others/workshops/Lechu_Slowacja_bix_workshop_08_2016/workshop/data/CEL_files/GSE17951/GSM449269.CEL", "Normal", "10/13/08 13:40:29" ], [ "18", "GSE17951_18", "/Users/marzec01/Desktop/others/workshops/Lechu_Slowacja_bix_workshop_08_2016/workshop/data/CEL_files/GSE17951/GSM449270.CEL", "Normal", "10/13/08 14:35:55" ], [ "19", "GSE17951_19", "/Users/marzec01/Desktop/others/workshops/Lechu_Slowacja_bix_workshop_08_2016/workshop/data/CEL_files/GSE17951/GSM449286.CEL", "Normal", "03/18/09 12:28:00" ], [ "20", "GSE17951_20", "/Users/marzec01/Desktop/others/workshops/Lechu_Slowacja_bix_workshop_08_2016/workshop/data/CEL_files/GSE17951/GSM449289.CEL", "Normal", "03/18/09 13:01:38" ], [ "21", "GSE17951_21", "/Users/marzec01/Desktop/others/workshops/Lechu_Slowacja_bix_workshop_08_2016/workshop/data/CEL_files/GSE17951/GSM449292.CEL", "Normal", "04/01/09 12:52:47" ], [ "22", "GSE17951_22", "/Users/marzec01/Desktop/others/workshops/Lechu_Slowacja_bix_workshop_08_2016/workshop/data/CEL_files/GSE17951/GSM449293.CEL", "Normal", "04/01/09 13:04:00" ] ];
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
