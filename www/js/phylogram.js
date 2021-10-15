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
  for ( var i = 0; i < str.length; ++i ) {
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

const LINKOUTS_BASE_URL = 'https://legumeinfo.org';

function lineWithLink(href, text) {
  return '<p><a href="' + href + '" target=_blank>' + text + '</a></p>';
}

var newickModule = require('biojs-io-newick');

// =============================================================

// Send a tree (in Newick format) to the phylotree chart
shinyjs.setPhylotree = function(args) {
  sessionStorage.setItem("newickTree", args[0]);
  var elementId = args[1];
  drawPhylotree(elementId);

  // this tree also defines the taxa for the Taxa chart
  setTaxa("taxa");
}

// Update the phylotree chart layout (vertical or radial)
shinyjs.setPhylotreeLayout = function(args) {
  sessionStorage.setItem("treeLayout", args[0]);
  var elementId = args[1];
  drawPhylotree(elementId);
}

// Phylotree can be global, as long as there is only one per page
gPhylotree = null;
gSelectedNode = null;

// Render the phylotree
function drawPhylotree(elementId) {
  // read the (Newick) tree
  var newick = sessionStorage.getItem("newickTree");
  gPhylotree = tnt.tree().data(tnt.tree.parse_newick(newick));

  // filter by taxa
  var taxa = sessionStorage.getItem("taxa");
  if (taxa !== null) {
    taxa = JSON.parse(taxa);
    var leaves = gPhylotree.root().find_all(function(node) {
      if (!node.is_leaf()) return false;
      var gensp = node.node_name().substring(0, 5);
      if (gensp.startsWith("USR")) gensp = "USR";
      return taxa.includes(gensp);
    });
    var subtree = gPhylotree.root().subtree(leaves);
    gPhylotree.data(subtree.data());
  }

  // set tree layout
  var width = window.innerWidth - 200;
  var layout = sessionStorage.getItem("treeLayout");
  if (layout === "Radial") {
    gPhylotree.layout(tnt.tree.layout.radial().width(width));
  } else {
    // "Vertical" or null
    gPhylotree.layout(tnt.tree.layout.vertical().width(width));
  }

  // set node colors and labels
  gPhylotree.node_display(tnt.tree.node_display.circle()
    .size(6)
    .fill(function(node) {
      var gen = node.node_name().substring(0, 3)
      if (gen == 'USR') return(genusToColor[gen])
      var gensp = node.node_name().substring(0, 5)
      if (gensp in genspToTaxon) return(taxonToColor(genspToTaxon[gensp]))
      else if (node.is_leaf() && !node.is_collapsed()) return(DEFAULT_COLOR);
      else if (isRootNode(node)) return 'black';
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
  setTimeout(highlightUserSequences, 1000);
}

function isRootNode(node) {
  return node.root_dist() == 0;
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
}

function highlightUserSequences() {
  d3.selectAll('text.tnt_tree_label')
    .filter((d) => d.name.toLowerCase().startsWith('usr'))
    .each(function(d) {
      d.bbox = this.getBBox();
    });
  var top = Infinity;
  d3.selectAll('g.tnt_tree_node')
    .filter((d) => d.name.toLowerCase().startsWith('usr'))
    .insert('svg:rect', ':first-child')
    .attr('x', (d) => {
      if (d.textAnchor === 'end') {
        // textAnchor was set dynamically for the radial layout
        return d.bbox.x + d.bbox.width + 10;
      }
      return d.bbox.x + 10;
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
  var distance = tree.scale_bar(30, 'pixel');
  var scale = d3.scale.linear()
    .domain([0, distance * width/30 ])
    .range([0, width]);
  var axis = d3.svg.axis()
    .scale(scale)
    .ticks(12)
    .orient('bottom');
  return axis;
}
function updateXAxis() {
  d3.selectAll('#phylotreeDistanceScale g.x.axis')
    .call(getXAxis(gPhylotree, window.innerWidth - 200));
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
    var node_name = node.node_name();
    var gensp = node_name.substring(0, 5);
    if (gensp in genspToTaxon) {
      var taxon = genspToTaxon[gensp];
      taxon = taxon.split(' ');
      var genus = taxon[0];
      var species = taxon[1];
      var url = LINKOUTS_BASE_URL + '/phylotree_links/' + genus + '/' + species + '/' + node_name + '/json';
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
        content += lineWithLink(LINKOUTS_BASE_URL + '/organism/' + genus + '/' + species,
          'View organism: ' + genus + ' ' + species); // + ' (' + common_name + ')'
        content += lineWithLink(LINKOUTS_BASE_URL + '/node/' + node.data()._id,
          'View feature: ' + node_name);
*/
        showDialog(node_name, content);
      });
    } else {
      showDialog(node_name, '<p>User sequence, no linkouts.</p>');
    }
  } else {
    // Internal node actions
    gSelectedNode = node;
    var internalNodeActions = '';
    var toggleText = (isCollapsed ? 'Expand' : 'Collapse');
    internalNodeActions += '<p><a onclick="doCollapseExpand();">' + toggleText + ' subtree at this node</a></p>';
    if (!isCollapsed) internalNodeActions += '<p><strike>Focus on subtree at this node</strike></p>';
    internalNodeActions += '<p><a onclick="doExport();">Export subtree in Newick format</a></p>';
    var url = LINKOUTS_BASE_URL + '/famreps_links';
    var query = 'famreps=' + gSelectedNode.get_all_leaves(true).map(n => n.node_name()).join(',');
    // add internal node linkouts, if any
    $.post({
      url: url,
      data: query,
      dataType: 'json'
    })
    .done(function(data) {
      $.each(data, function(d, obj) {
        internalNodeActions += lineWithLink(obj.href, obj.text);
      });
    })
    .always(function() {
      showDialog('Interior node', internalNodeActions);
    });
  }
}

function doCollapseExpand() {
  $('#dlg1').dialog('destroy');

  gSelectedNode.toggle();
  gPhylotree.update();
  setTimeout(() => { updateXAxis(); }, 500);
  highlightUserSequences();

  gSelectedNode = null;
}

function doExport() {
  $('#dlg1').dialog('destroy');

  var content = newickModule.parse_json(gSelectedNode.data());
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
      },
      {
        text: 'Dismiss',
        click: function() {
          $(this).dialog('close');
        }
      }
    ]
  });

  gSelectedNode = null;
}

// =============================================================

// Set the Taxa and Legend chart data (Newick tree and selected taxa) from sessionStorage
function setTaxa(elementId) {
  var elementTag = "#" + elementId;

  // read the (Newick) tree
  var newick = sessionStorage.getItem("newickTree");
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
      leaves.push(node.name);
    } else {
      for (i in node.children) recParse(node.children[i]);
    }
  }
  recParse(root);

  var count = {};
  for (l in leaves) {
    leaf = leaves[l];
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
      count[taxon] = count[taxon] + 1;
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
  if (taxa.length == data.length) {
    sessionStorage.removeItem("taxa");
  } else {
    sessionStorage.setItem("taxa", JSON.stringify(taxa));
  }

  // Redraw MSA viewer and phylotree, filtered by selected taxa
  // (TODO: should really notify all registered event listeners)
  drawMSAView();
  drawPhylotree("phylotree");
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

// =============================================================

// MSA view can be global, as long as there is only one per page
gMsaView = null;

// Send multiple sequence alignment to the MSA view
shinyjs.setMSA = function(args) {
  // Replace special characters in sequence names to prevent FASTA parser from treating them as delimiters
  var seqs = msa.io.fasta.parse(args[0].replaceAll("|", ".").replaceAll(":", "."));
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

  drawMSAView();
}

// Render the MSA view, filtered by selected taxa
function drawMSAView() {
  var taxa = sessionStorage.getItem("taxa");
  if (taxa !== null) taxa = JSON.parse(taxa);
  for (var s = 0; s < gMsaView.seqs.length; s++) {
    var gensp = gMsaView.seqs.at(s).get('name').substring(0, 5);
    if (gensp.startsWith("USR")) gensp = "USR";
    gMsaView.seqs.at(s).set('hidden', !(taxa === null || taxa.includes(gensp)));
  }

  gMsaView.render();
}

// =============================================================

