// (C) Wolfgang Huber 2010-2011

// Script parameters - these are set up by R in the function 'writeReport' when copying the 
//   template for this script from arrayQualityMetrics/inst/scripts into the report.

var highlightInitial = [ false, false, true, true, false, true, false, false, false, false, false, false, false, false, true, false ];
var arrayMetadata    = [ [ "1", "GSE55945_1", "/Users/marzec01/Desktop/others/workshops/Lechu_Slowacja_bix_workshop_08_2016/workshop/data/CEL_files/GSE55945/GSM1348933_011508_HGU133_PLUS_2.0_MS_36D6.CEL", "Tumour", "01/16/08 12:36:10" ], [ "2", "GSE55945_2", "/Users/marzec01/Desktop/others/workshops/Lechu_Slowacja_bix_workshop_08_2016/workshop/data/CEL_files/GSE55945/GSM1348934_110607_HGU133_PLUS_2.0_MS_36C1.CEL", "Tumour", "12/06/07 14:45:27" ], [ "3", "GSE55945_3", "/Users/marzec01/Desktop/others/workshops/Lechu_Slowacja_bix_workshop_08_2016/workshop/data/CEL_files/GSE55945/GSM1348935_011508_HGU133_PLUS_2.0_MS_36D7.CEL", "Tumour", "01/16/08 12:59:55" ], [ "4", "GSE55945_4", "/Users/marzec01/Desktop/others/workshops/Lechu_Slowacja_bix_workshop_08_2016/workshop/data/CEL_files/GSE55945/GSM1348939_092707_HGU133_PLUS_2.0_NUGEN_TEST07.CEL", "Tumour", "09/27/07 16:26:48" ], [ "5", "GSE55945_5", "/Users/marzec01/Desktop/others/workshops/Lechu_Slowacja_bix_workshop_08_2016/workshop/data/CEL_files/GSE55945/GSM1348940_110607_HGU133_PLUS_2.0_MS_36C8.CEL", "Tumour", "12/06/07 16:28:16" ], [ "6", "GSE55945_6", "/Users/marzec01/Desktop/others/workshops/Lechu_Slowacja_bix_workshop_08_2016/workshop/data/CEL_files/GSE55945/GSM1348941_110807_HGU133_PLUS_2.0_MS_36A1.CEL", "Tumour", "11/08/07 15:18:24" ], [ "7", "GSE55945_7", "/Users/marzec01/Desktop/others/workshops/Lechu_Slowacja_bix_workshop_08_2016/workshop/data/CEL_files/GSE55945/GSM1348942_110607_HGU133_PLUS_2.0_MS_36C2.CEL", "Tumour", "12/06/07 15:01:43" ], [ "8", "GSE55945_8", "/Users/marzec01/Desktop/others/workshops/Lechu_Slowacja_bix_workshop_08_2016/workshop/data/CEL_files/GSE55945/GSM1348944_011508_HGU133_PLUS_2.0_MS_36C9.CEL", "Tumour", "01/15/08 14:41:11" ], [ "9", "GSE55945_9", "/Users/marzec01/Desktop/others/workshops/Lechu_Slowacja_bix_workshop_08_2016/workshop/data/CEL_files/GSE55945/GSM1348945_092707_HGU133_PLUS_2.0_NUGEN_TEST05.CEL", "Tumour", "09/27/07 15:51:10" ], [ "10", "GSE55945_10", "/Users/marzec01/Desktop/others/workshops/Lechu_Slowacja_bix_workshop_08_2016/workshop/data/CEL_files/GSE55945/GSM1348946_110607_HGU133_PLUS_2.0_MS_36C6.CEL", "Normal", "12/06/07 16:02:45" ], [ "11", "GSE55945_11", "/Users/marzec01/Desktop/others/workshops/Lechu_Slowacja_bix_workshop_08_2016/workshop/data/CEL_files/GSE55945/GSM1348947_110607_HGU133_PLUS_2.0_MS_36C7.CEL", "Normal", "12/06/07 16:15:57" ], [ "12", "GSE55945_12", "/Users/marzec01/Desktop/others/workshops/Lechu_Slowacja_bix_workshop_08_2016/workshop/data/CEL_files/GSE55945/GSM1348949_011508_HGU133_PLUS_2.0_MS_36D3.CEL", "Normal", "01/16/08 11:56:10" ], [ "13", "GSE55945_13", "/Users/marzec01/Desktop/others/workshops/Lechu_Slowacja_bix_workshop_08_2016/workshop/data/CEL_files/GSE55945/GSM1348950_011508_HGU133_PLUS_2.0_MS_36D4.CEL", "Normal", "01/16/08 12:08:25" ], [ "14", "GSE55945_14", "/Users/marzec01/Desktop/others/workshops/Lechu_Slowacja_bix_workshop_08_2016/workshop/data/CEL_files/GSE55945/GSM1348951_011508_HGU133_PLUS_2.0_MS_36D5.CEL", "Normal", "01/16/08 12:19:39" ], [ "15", "GSE55945_15", "/Users/marzec01/Desktop/others/workshops/Lechu_Slowacja_bix_workshop_08_2016/workshop/data/CEL_files/GSE55945/GSM1348952_110807_HGU133_PLUS_2.0_MS_36A4.CEL", "Normal", "11/08/07 15:41:51" ], [ "16", "GSE55945_16", "/Users/marzec01/Desktop/others/workshops/Lechu_Slowacja_bix_workshop_08_2016/workshop/data/CEL_files/GSE55945/GSM1348953_110807_HGU133_PLUS_2.0_MS_36A5.CEL", "Normal", "11/08/07 15:53:28" ] ];
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
