// =============================================================

var DEFAULT_COLOR = '#d3d3d3';

var OVERRIDES = {
  'arabidopsis thaliana': DEFAULT_COLOR, // otherwise conflicts with ara = Arachis
  'arachis ipaensis': '#aaab00',
  '*user sequences': '#990000' // TODO: use genusToColor?
}

moreBrewerColors = chroma.brewer.Set2;

genusToColor = {
  // for legume genera
  //'aes': '#......',
  'api': moreBrewerColors[0],
  'ara': '#bcbd22',
  'caj': '#ffbb78',
  'cha': moreBrewerColors[5],
  'cic': '#2ca02c',
  'gly': '#1f77b4',
  'len': '#98df8a',
  'lot': '#17becf',
  'lup': '#ff9896',
  'med': '#8c564b',
  'pha': '#e377c2',
  'pis': '#f7b6d2',
  'tri': moreBrewerColors[2],
  'USR': '#990000',
  'vic': moreBrewerColors[4],
  'vig': '#d62728'
}

genspToTaxon = {
  // legumes
  //'aesev': 'Aeschynomene evenia',
  'apiam': 'Apios americana',
  'aradu': 'Arachis duranensis',
  'arahy': 'Arachis hypogaea',
  'araip': 'Arachis ipaensis',
  'cajca': 'Cajanus cajan',
  'chafa': 'Chamaecrista fasciculata',
  'cicar': 'Cicer arietinum',
  'cicec': 'Cicer echinospermum',
  'cicre': 'Cicer reticulatum',
  'glyma': 'Glycine max',
  'lencu': 'Lens culinaris',
  'lotja': 'Lotus japonicus',
  'lupal': 'Lupinus albus',
  'lupan': 'Lupinus angustifolius',
  'medsa': 'Medicago sativa',
  'medtr': 'Medicago truncatula',
  'phalu': 'Phaseolus lunatus',
  'phavu': 'Phaseolus vulgaris',
  'pissa': 'Pisum sativum',
  'tripr': 'Trifolium pratense',
  'trire': 'Trifolium repens',
  'vicfa': 'Vicia faba',
  'vigan': 'Vigna angularis',
  'vigra': 'Vigna radiata',
  'vigun': 'Vigna unguiculata',
  // non-legumes
  'arath': 'Arabidopsis thaliana',
  'cucsa': 'Cucumis sativus',
  'prupe': 'Prunus persica',
  'solly': 'Solanum lycopersicum',
  'vitvi': 'Vitis vinifera'
}
taxonToGensp = function(taxon) {
  if (taxon == '*User sequences') return 'USR';
  var parts = taxon.toLowerCase().split(' ');
  return parts[0].substring(0, 3) + parts[1].substring(0, 2);
}

function fnv32a(str, hashSize) {
  /* a consistent hashing algorithm
     https://gist.github.com/vaiorabbit/5657561
     http://isthe.com/chongo/tech/comp/fnv/#xor-fold
  */
  var FNV1_32A_INIT = 0x811c9dc5;
  var hval = FNV1_32A_INIT;
  for (var i = 0; i < str.length; i++) {
    hval ^= str.charCodeAt(i);
    hval += (hval << 1) + (hval << 4) + (hval << 7) + (hval << 8) + (hval << 24);
  }
  return (hval >>> 0) % hashSize;
}

function taxonToColor(taxon) {
  var colorCache = {};
  var LIGHTNESS_FACTOR = 1;
  var MIN_LIGHTNESS = 0.3;

  taxon = taxon.toLowerCase();
  if (OVERRIDES === undefined) {
    OVERRIDES = {};
  }
  if (OVERRIDES[taxon] !== undefined) {
    return OVERRIDES[taxon];
  }
  if (colorCache[taxon] !== undefined) {
    return colorCache[taxon];
  }
  var color = null;
  var parts = taxon.split(' ');
  var genus = parts[0];
  var gen = genus.substring(0, 3);
  var species = parts[1];

  if (gen in genusToColor) {
    var genusColor = genusToColor[gen];
    var hue = chroma(genusColor).hsl()[0];
    var lightness = MIN_LIGHTNESS + (fnv32a(species, 1000) / 1000) * (1 - 2*MIN_LIGHTNESS);
    color = chroma(hue, 1, lightness*LIGHTNESS_FACTOR, 'hsl').hex();
  } else {
    color = DEFAULT_COLOR;
  }
  colorCache[taxon] = color;
  return color;
}

// =============================================================

var legumeMine = new intermine.Service({ root: 'https://mines.legumeinfo.org/legumemine/service/' });

const LINKOUTS_BASE_URL = 'https://linkouts.services.legumeinfo.org';

function lineWithLink(href, text) {
  return '<p><a href="' + href + '" target=_blank>' + text + '</a></p>';
}

var newickModule = require('biojs-io-newick');

const TIMEOUT_DURATION_MSEC = 500;
const LABEL_HEIGHT_PIXELS = 10;
const AXIS_TICKS = 12;
const AXIS_SAMPLE_PIXELS = 30;

// =============================================================

// Send a tree (in Newick format) to the phylotree chart
shinyjs.setPhylotree = function(args) {
  sessionStorage.clear();
  sessionStorage.setItem("newickTree", args[0]);
  var elementId = args[1];
  if (args.length > 2) {
    sessionStorage.setItem("highlightedProteins", args[2])
  }
  drawPhylotree(elementId);

  // this tree also defines the taxa for the Taxa chart
  setTaxa("taxa");
}

// Update the phylotree chart layout (vertical or radial)
function getPhylotreeWidth() {
  return window.innerWidth - 200;
}

function updateLayout(layoutType, width) {
  var layout = (layoutType === "Radial" ? tnt.tree.layout.radial() : tnt.tree.layout.vertical()); // otherwise "Vertical" or null
  gPhylotree.layout(layout.width(width));
}

shinyjs.setPhylotreeLayout = function(args) {
  var layout = args[0];
  // var elementId = args[1]; // unused

  sessionStorage.setItem("treeLayout", layout);
  updateLayout(layout, getPhylotreeWidth());
  gPhylotree.update();
  setTimeout(() => { updateXAxis(); }, TIMEOUT_DURATION_MSEC);
  setTimeout(highlightSequences, TIMEOUT_DURATION_MSEC);
}

// Phylotree can be global, as long as there is only one per page
gPhylotree = null;
gSelectedNode = null;
gFocusedNode = null;

// Render the phylotree
function drawPhylotree(elementId) {
  // read the (Newick) tree
  var newick = sessionStorage.getItem("newickTree");
  if (newick == null) return;
  gPhylotree = tnt.tree().data(tnt.tree.parse_newick(newick));

  // filter by taxa
  var taxa = sessionStorage.getItem("taxa");
  if (taxa !== null) {
    taxa = JSON.parse(taxa);
    var rootNode = (gFocusedNode == null ? gPhylotree.root() : gFocusedNode);
    var leaves = rootNode.find_all(function(node) {
      if (!node.is_leaf()) return false;
      var gensp = node.node_name().substring(0, 5);
      if (gensp.startsWith("USR")) gensp = "USR";
      return taxa.includes(gensp);
    }, true);
    var subtree = gPhylotree.root().subtree(leaves, true);
    gPhylotree.data(subtree.data());
  }

  // set tree layout
  var layout = sessionStorage.getItem("treeLayout");
  var width = getPhylotreeWidth();
  updateLayout(layout, width);

  // set node colors and labels
  var ssn = sessionStorage.getItem("showSingletonNodes");
  var showSingletonNodes = (ssn == undefined || ssn == 'true');
  gPhylotree.node_display(tnt.tree.node_display.circle()
    .size(function(node) {
      return (!isRootNode(node) && isSingletonNode(node) && !showSingletonNodes) ? 0 : 6;
    })
    .fill(function(node) {
      var gen = node.node_name().substring(0, 3)
      if (gen == 'USR') return(genusToColor[gen])
      var gensp = node.node_name().substring(0, 5)
      if (gensp in genspToTaxon) return(taxonToColor(genspToTaxon[gensp]))
      else if (node.is_leaf() && !node.is_collapsed()) return(DEFAULT_COLOR);
      else if (isRootNode(node)) return 'black';
      else if (isSingletonNode(node)) return 'gray';
      return 'white'; // internal nodes
    })
  );
  gPhylotree.label(tnt.tree.label.text()
    .height(16)
    .text(function(node) {
      if (node.is_collapsed()) return(node.get_all_leaves(true).length + ' leaf nodes');
      if (!node.is_leaf()) return('');
      return node.node_name();
    })
  );

  var treeDiv = document.getElementById(elementId);
  treeDiv.innerHTML = ""; // to clear it
  gPhylotree(treeDiv);

  gPhylotree.on('click', node => onTreeNodeClick(node));

  drawDistanceScale(width);

  // wrap in timeout in order to wait to render the tree before computing label bounding boxes
  setTimeout(highlightSequences, TIMEOUT_DURATION_MSEC);
}

function isRootNode(node) {
  return node.root_dist() == 0;
}
function isSingletonNode(node) {
  if (gFocusedNode == null) return false;
  var nn = node.find_all(function(n) { return n.id() === gFocusedNode.id(); });
  if (nn.length !== 1) return false;
  var children = node.children(true);
  return children != undefined && children.length === 1;
}
function isFocusedNode(node) {
  if (gFocusedNode == null) return false;
  return node.id() == gFocusedNode.id();
}

function drawDistanceScale(width) {
  // use TnT Tree API to calculate intergenic distance over some screen distance.
  $("#phylotreeDistanceScale").empty(); // to clear it
  d3.select("#phylotreeDistanceScale")
    .attr('width', width)
    .attr('height', 40)
    .append('g')
      .attr('transform', 'translate(20, 20)')
      .attr('class', 'x axis')
    .call(getXAxis(gPhylotree, width)
  );
  setXAxisAttributes();
}

function nontrivialMatch(str1, str2) {
  if (str1 == '' || str2 == '') return false;
  return str1.includes(str2);
}

// Highlight if (a) user sequence, (b) specified in URL, (c) selected in MSA
function highlightSequences() {
  var highlightedProteins = sessionStorage.getItem("highlightedProteins")
  var msaSelectedProteins = msaSelected();
  d3.selectAll('text.tnt_tree_label')
    .filter((d) => d.name.toLowerCase().startsWith('usr') || nontrivialMatch(highlightedProteins, d.name) || nontrivialMatch(msaSelectedProteins, d.name))
    .each(function(d) {
      d.bbox = this.getBBox();
      if (d.bbox.x < -1) d.bbox.x += d.bbox.width; // to correctly highlight sequences on the left side of the radial layout
    });
  var top = Infinity;
  d3.selectAll('g.tnt_tree_node')
    .filter((d) => d.name.toLowerCase().startsWith('usr') || nontrivialMatch(highlightedProteins, d.name) || nontrivialMatch(msaSelectedProteins, d.name))
    .insert('svg:rect', ':first-child')
    .attr('x', (d) => {
      if (d.textAnchor === 'end') {
        // textAnchor was set dynamically for the radial layout
        return d.bbox.x + d.bbox.width + LABEL_HEIGHT_PIXELS;
      }
      return d.bbox.x + LABEL_HEIGHT_PIXELS;
    })
    .attr('y', (d) => d.bbox.y/2)
    .attr('width', (d) => d.bbox.width + 2)
    .attr('height', (d) => d.bbox.height + 1)
    .attr('class', 'hilite-node')
    .each(function() {
      var offset = $(this).offset();
      if (offset.top > 0 && offset.top < top) {
        top = offset.top;
      }
    });
}

function getXAxis(tree, width) {
  var distance = tree.scale_bar(AXIS_SAMPLE_PIXELS, 'pixel');
  var scale = d3.scale.linear()
    .domain([0, distance * width/AXIS_SAMPLE_PIXELS ])
    .range([0, width]);
  var axis = d3.svg.axis()
    .scale(scale)
    .ticks(AXIS_TICKS)
    .orient('bottom');
  return axis;
}
function updateXAxis() {
  d3.selectAll('#phylotreeDistanceScale g.x.axis')
    .call(getXAxis(gPhylotree, getPhylotreeWidth()));
  setXAxisAttributes();
}
function setXAxisAttributes() {
  d3.select("#phylotreeDistanceScale").select(".domain")
    .attr('fill', 'none')
    .attr('stroke', 'black');
  d3.select("#phylotreeDistanceScale").selectAll("g.tick > line")
    .attr('stroke', 'black');
  d3.selectAll("#phylotreeDistanceScale").selectAll("g.tick > text")
    .style('font-size', '10px');
}

// User clicked on link to sequence -> scroll to its phylotree node
function onScrollToHilite(seqName) {
  var node = gPhylotree.root().find_node(node => {
    return node.node_name() === seqName;
  });
  var nodeId = node.property('_id');
  var selector = '#tnt_tree_node_phylotree_' + nodeId;
  var offset = $(selector).offset();
  $('html,body').scrollTop(offset.top - 100);
}

function showDialog(title, content) {
  $('<div id="dlg1">' + content + '</div>').dialog({
    title: title,
    width: 512,
    modal: true
  });
}

// User clicked a phylotree node (internal or leaf)
function onTreeNodeClick(node) {
  const isCollapsed = node.is_collapsed();
  if (!isCollapsed && node.is_leaf()) {
    // Leaf nodes: post dialog with linkouts
    const node_fullname = node.node_name();
    var node_name = node_fullname;
    if (node_name.startsWith('USR.')) node_name = node_name.substring(4);
    var gensp = node_name.substring(0, 5);
    if (gensp in genspToTaxon) {
      // Get corresponding gene identifier from LegumeMine
      var query = {
        "from": "Protein",
        "select": [
          "genes.primaryIdentifier",
        ],
        "where": [
          {
            "path": "primaryIdentifier",
            "op": "=",
            "value": node_name,
            "code": "A"
          }
        ]
      };
      legumeMine.rows(query).then(function(rows) {
        if (rows.length == 0) {
          showDialog(node_name, 'Non-legume, no linkouts available.');
          return;
        }
        var gene_id = rows[0][0];
        var url = LINKOUTS_BASE_URL + '/gene_linkouts?genes=' + gene_id;
        $.getJSON(url, function(data) {
          var content = '';
          if (data.length > 0) {
            $.each(data, function(d, obj) {
              content += lineWithLink(obj.href, obj.text);
            });
          } else {
            content += 'No linkouts found.';
          }
/*
          // for the organism and feature links (unused)
          var taxon = genspToTaxon[gensp];
          taxon = taxon.split(' ');
          var genus = taxon[0];
          var species = taxon[1];
          content += lineWithLink(LINKOUTS_BASE_URL + '/organism/' + genus + '/' + species,
            'View organism: ' + genus + ' ' + species); // + ' (' + common_name + ')'
          content += lineWithLink(LINKOUTS_BASE_URL + '/node/' + node.data()._id,
            'View feature: ' + node_name);
*/
          showDialog(gene_id, content);
        });
      });
    } else {
      // phylogram leaf nodes should always have a gensp, so we may never reach here, but post a message if we do
      showDialog(node_fullname, 'No linkouts available.');
    }
  } else {
    // Internal node actions
    gSelectedNode = node;
    var internalNodeActions = '';
    var toggleText = (isCollapsed ? 'Expand' : 'Collapse');
    if (gFocusedNode == null || !(isSingletonNode(node) || isRootNode(node))) {
      internalNodeActions += '<p><a onclick="doCollapseExpand();">' + toggleText + ' subtree at this node</a></p>';
    }
    if (!(isFocusedNode(node) || isCollapsed || isSingletonNode(node) || isRootNode(node))) {
      internalNodeActions += '<p><a onclick="doFocus();">Focus on subtree at this node</a></p>';
    }
    internalNodeActions += '<p><a onclick="doExport();">Export subtree in Newick format</a></p>';

    var node_names = gSelectedNode.get_all_leaves(true).map(n => n.node_name());
    var query = {
      "from": "Protein",
      "select": [
        "genes.primaryIdentifier",
      ],
      "where": [
        {
          "path": "primaryIdentifier",
          "op": "ONE OF",
          "values": node_names,
          "code": "A"
        }
      ]
    };
    legumeMine.rows(query).then(function(rows) {
      var genes = rows.join(",");
      if (genes.length > 0) {
        var url = 'https://gcv.legumeinfo.org/gcv2/gene;lis=' + genes;
        internalNodeActions += lineWithLink(url, 'GCV Multiple Alignment');
      }
      showDialog('Interior node', internalNodeActions);
    });
  }
}

function doCollapseExpand() {
  $('#dlg1').dialog('destroy');

  gSelectedNode.toggle();
  gPhylotree.update();
  gSelectedNode = null;

  var rootNode = gPhylotree.root();
  var subtree = rootNode.subtree(rootNode.get_all_leaves(false), true);
  var newick = newickModule.parse_json(subtree.data());
  sessionStorage.setItem("subtree", newick);
  sessionStorage.setItem("taxa", JSON.stringify(Object.keys(countTaxa(newick)).map(taxonToGensp)));
  setTaxa("taxa");
  // no need to call drawPhylotree() for collapse/expand
  drawMSAView();

  setTimeout(() => { updateXAxis(); }, TIMEOUT_DURATION_MSEC);
  highlightSequences();
}

function doFocus() {
  $('#dlg1').dialog('destroy');

  Shiny.setInputValue('toggleSubtree', true);
  gFocusedNode = gSelectedNode; // to maintain the focused subtree
  gSelectedNode = null;

  // expand all collapsed nodes on focus, so as not to lose any
  var collapsedNodes = gPhylotree.root().find_all(function(node) { return node.is_collapsed(); }, true);
  var collapsedNodeIds = collapsedNodes.map(function(node) { return node.id(); } );
  collapsedNodes.forEach(function(node) { node.toggle(); });
  gPhylotree.update();

  var subtree = gPhylotree.root().subtree(gFocusedNode.get_all_leaves(true), true);
  var newick = newickModule.parse_json(subtree.data());
  sessionStorage.setItem("subtree", newick);
  sessionStorage.setItem("taxa", JSON.stringify(Object.keys(countTaxa(newick)).map(taxonToGensp)));

  setTaxa("taxa");
  drawPhylotree("phylotree");
  drawMSAView();
}

shinyjs.showSingletonNodes = function(args) {
  sessionStorage.setItem("showSingletonNodes", args[0]);
  // var elementId = args[1]; // unused

  // add css classes to singleton nodes
  d3.selectAll('#phylotree .inner.tnt_tree_node')
  .classed('singleton-node', (d) =>  {
    d.singleton = d.children && d.children.length === 1;
    return d.singleton;
  })
  .filter((d) => d.singleton)
  .classed('singleton-node-visible', args[0]);
}

shinyjs.clearSubtreeFocus = function(args) {
  var elementId = args[0];
  gPhylotree = null;
  gSelectedNode = null;
  gFocusedNode = null;
  // show singleton nodes, to prevent them from disappearing (and while toggleSubtree is still true)
  Shiny.setInputValue("showSingletonNodes", true);
  $("#showSingletonNodes").prop("checked", true); // and set its checkbox
  Shiny.setInputValue("toggleSubtree", false);
  sessionStorage.removeItem("subtree");
  sessionStorage.removeItem("taxa");

  setTaxa("taxa");
  drawPhylotree(elementId);
  drawMSAView();
}

function doExport() {
  $('#dlg1').dialog('destroy');

  var content = newickModule.parse_json(gSelectedNode.data());
  gSelectedNode = null;
  $('<div id="dlg2" style="font-size: 12px">' + content + '</div>').dialog({
    title: 'Subtree in Newick format',
    width: 600,
    height: 200,
    modal: true,
    buttons: [
      {
        text: 'Copy to Clipboard',
        click: function() {
          var textArea = document.createElement("textarea");
          textArea.value = content;
          textArea.style.position = "fixed";
          textArea.style.left = "-999999px";
          textArea.style.top = "-999999px";
          this.appendChild(textArea); // make it a child of the dialog to preserve modality
          textArea.select();
          return new Promise((res, rej) => {
            document.execCommand('copy') ? res() : rej();
            textArea.remove();
          });
        }
      }
    ]
  });
}

// =============================================================

// Set the Taxa and Legend chart data (Newick tree and selected taxa) from sessionStorage
function setTaxa(elementId) {
  var elementTag = "#" + elementId;

  // read the (Newick) tree
  var newick = sessionStorage.getItem("subtree");
  if (newick == null) newick = sessionStorage.getItem("newickTree");
  var count = countTaxa(newick);

  // read currently selected taxa
  var taxa = sessionStorage.getItem("taxa");
  if (taxa !== null) taxa = JSON.parse(taxa);

  // transform to objects expected by nvd3 chart
  var data = _.map(Object.keys(count).sort(), k => {
    d = {
      label: k,
      value: count[k],
      color: taxonToColor(k),
      disabled: !(taxa === null || taxa.includes(taxonToGensp(k)))
    };
    return d;
  });
  nv.addGraph(function() {
    var chart = nv.models.pieChart()
      .x(function(d) { return d.label })
      .y(function(d) { return d.value })
      .color(function(d) { return d.color })
      .showLabels(true)
      .labelThreshold(0.05)
      .labelType("percent")
      .donut(true)
      .donutRatio(0.35)
      .valueFormat(d3.format("d"))
      .legendPosition("right")
      .duration(50)
    ;
    d3.select(elementTag)
      .datum(data)
      .transition().duration(50)
      .call(chart)
    ;
    nv.utils.windowResize(chart.update);

    chart.dispatch.on('stateChange', evt => {
      // alert("stateChange: " + JSON.stringify(evt));
      handleTaxaSelection(elementTag);
    });
    chart.dispatch.on('renderEnd', () => {
      doResize();
    });

    return chart;
  });
}

// Tabulate the tree's leaf nodes by taxon
function countTaxa(tree) {
  var root = tnt.tree.parse_newick(tree);
  var leaves = [];
  function recParse(node) {
    if (node.children === undefined) {
      if (isNaN(node.name)) leaves.push(node.name);
    } else {
      for (var childNode of node.children) recParse(childNode);
    }
  }
  recParse(root);

  var count = {};
  for (var leaf of leaves) {
    gen = leaf.substring(0, 3);
    if (gen == 'USR') {
      taxon = '*User sequences';
    } else {
      gensp = leaf.substring(0, 5);
      if (gensp in genspToTaxon) {
        taxon = genspToTaxon[gensp];
      } else {
        taxon = gensp;
      }
    }
    if (taxon in count) {
      ++count[taxon];
    } else {
      count[taxon] = 1;
    }
  }
  return count;
}

// User clicked on the legend to toggle a taxon
handleTaxaSelection = function(elementTag) {
  // first determine the selected taxa (in gensp format)
  taxa = [];
  var legend = d3.select(elementTag);
  var data = legend.datum();
  data.forEach(function(series) {
    if (series.disabled === undefined || !series.disabled) {
      taxa.push(taxonToGensp(series.label));
    }
  });
  if (sessionStorage.getItem("subtree") == null && taxa.length == data.length) {
    sessionStorage.removeItem("taxa");
  } else {
    sessionStorage.setItem("taxa", JSON.stringify(taxa));
  }

  // Redraw phylotree and MSA viewer, filtered by selected taxa
  // (TODO: should really notify all registered event listeners)
  drawPhylotree("phylotree");
  drawMSAView();
}

// Reset all taxa to enabled/visible
shinyjs.resetTaxa = function(args) {
  elementId = args[0];
  elementTag = "#" + elementId;
  var legend = d3.select(elementTag);
  var data = legend.datum();
  data.forEach(function(series) { series.disabled = false });
  sessionStorage.removeItem("taxa");

  handleTaxaSelection(elementTag);
  doResize();
}

// Force a resize, to redraw the chart
function doResize() {
  window.dispatchEvent(new Event("resize"));
}

// Export phylotree (including its distance scale) to SVG file
shinyjs.takePhylotreeSnapshot = function(args) {
  // Distance scale
  var svgAxis = document.getElementById("phylotreeDistanceScale");
  var h1 = parseInt(svgAxis.getAttribute("height"));
  // Phylotree
  const tnt_st_id = "tnt_st_" + args[0]; // id of the svg element's child
  var svgPhylotree = document.getElementById(tnt_st_id).parentElement;
  var h2 = parseInt(svgPhylotree.getAttribute("height"));

  var serializer = new XMLSerializer();
  // Wrap the distance scale and phylotree in an overall SVG element
  var source = "<svg id='combined'>";
    source += serializer.serializeToString(svgAxis);
    source += serializer.serializeToString(svgPhylotree);
  source += "</svg>";

  // Shift svgPhylotree down by h1
  source = source.replace("height=\"" + h2 + "\" fill=\"none\"", "height=\"" + (h2 + h1) + "\" fill=\"none\"");
  source = source.replace("translate(20,20)", "translate(20," + (20 + h1) + ")");

  source = source.replace(/(\w+)?:?xlink=/g, 'xmlns:xlink='); // Fix root xlink without namespace
  source = source.replace(/ns\d+:href/g, 'xlink:href'); // Safari NS namespace fix.
  if (!source.match(/^<svg[^>]+xmlns="http\:\/\/www\.w3\.org\/2000\/svg"/)) {
    source = source.replace(/^<svg/, '<svg xmlns="http://www.w3.org/2000/svg"');
  }

  var preface = '<?xml version="1.0" standalone="no"?>\r\n';
  var svgBlob = new Blob([preface, source], { type: "image/svg+xml;charset=utf-8" });
  var svgUrl = URL.createObjectURL(svgBlob);
  var downloadLink = document.createElement("a");
  downloadLink.href = svgUrl;
  downloadLink.download = "phylotree.svg";
  document.body.appendChild(downloadLink);
  downloadLink.click();
  document.body.removeChild(downloadLink);
}

// =============================================================

// MSA view can be global, as long as there is only one per page
gMsaView = null;

// Send multiple sequence alignment to the MSA view
shinyjs.setMSA = function(args) {
  // Replace special characters in sequence names to prevent FASTA parser from treating them as delimiters
  var seqs = msa.io.fasta.parse(args[0].replaceAll("|", ".").replaceAll(":", "."));
  sessionStorage.setItem('msaSeqs', JSON.stringify(seqs));
  var elementId = args[1];
  var msaElement = document.getElementById(elementId);
  msaElement.innerHTML = "";

  // TODO: persist settings across sessions, possibly in localStorage?
  gMsaView = msa({
    el: msaElement,
    bootstrapMenu: true,
    seqs: seqs,
    vis: {
      overviewbox: false,
      labelId: false
    },
    zoomer: {
      labelNameLength: 150,
      labelFontsize: 9,
      autoResize: true // only works for width
    }
  });
  // to effect autoResize for height (= show more or less rows when manually resizing)
  gMsaView.g.zoomer.autoHeight();
  $('#msa').height(225);
  // to redraw phylotree on MSA selection
  gMsaView.g.on("row:click", function(data) {
    drawPhylotree("phylotree");
  });
  gMsaView.g.selcol.on("reset", function() {
    drawPhylotree("phylotree");
  });

  drawMSAView();
}

// Render the MSA view, filtered by selected taxa and subtree
// (always call drawPhylotree() before drawMSAView() to ensure that gPhylotree is set)
function drawMSAView() {
  var taxa = sessionStorage.getItem("taxa");
  if (taxa !== null) taxa = JSON.parse(taxa);
  var rootNode = (gFocusedNode == null ? gPhylotree.root() : gFocusedNode);
  var visibleNodes = rootNode.get_all_leaves().map(function(node) { return node.node_name(); });

  // Reload the sequences and set their 'hidden' attribute according to current visibility,
  // then call seqs.reset() which is much faster than set("hidden", !visible) for each one
  var seqs = JSON.parse(sessionStorage.getItem('msaSeqs'));
  for (seq of seqs) {
    var seqName = seq["name"];
    var gensp = seqName.substring(0, 5);
    if (gensp.startsWith("USR")) gensp = "USR";

    var visible = (taxa === null || taxa.includes(gensp));
    var isAncestor = (gFocusedNode != null && gFocusedNode.find_node_by_name(seqName) != null);
    visible &&= ((gFocusedNode === null || isAncestor) && visibleNodes.includes(seqName));
    seq["hidden"] = !visible;
  }
  gMsaView.seqs.reset(seqs);

  gMsaView.render();
}

// Return the selected rows (sequence names), as a comma-separated string
function msaSelected() {
  var ss = '';
  for (var i = 0; i < gMsaView.g.selcol.length; i++) {
    var seqId = gMsaView.g.selcol.at(i).get('seqId');
    var seqName = gMsaView.seqs._byId[seqId].get('name');
    if (i > 0) ss += ',';
    ss += seqName;
  }
  return ss;
}

// =============================================================

