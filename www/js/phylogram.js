// =============================================================

var DEFAULT_COLOR = '#d3d3d3';

var OVERRIDES = {
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

// TODO: standardize function names
shinyjs.displayPhylotree = function(args) {
  var newick = args[0];
  localStorage.setItem("newickTree", newick);
  var elementId = args[1];
  drawPhylotree(elementId);
}

shinyjs.changePhylotreeLayout = function(args) {
  var useVerticalLayout = args[0];
  localStorage.setItem("verticalLayout", useVerticalLayout);
  var elementId = args[1];
  drawPhylotree(elementId);
}

function filterPhylotree(taxa, elementId) {
  localStorage.setItem("taxa", taxa);
  drawPhylotree(elementId);
}

function drawPhylotree(elementId) {
  // read the (Newick) tree
  var newick = localStorage.getItem("newickTree");
  var tree = tnt.tree().data(tnt.tree.parse_newick(newick));

  // filter by taxa
  var taxa = localStorage.getItem("taxa");
  if (taxa !== null) {
    var leaves = tree.root().find_all(function(node) {
      if (!node.is_leaf()) return false;
      var gensp = node.node_name().substring(0, 5);
      if (gensp.startsWith("USR")) gensp = "USR";
      return taxa.includes(gensp);
    });
    var subtree = tree.root().subtree(leaves);
    tree.data(subtree.data());
  }

  // set tree layout
  var width = window.innerWidth - 200;
  var useVerticalLayout = (localStorage.getItem("verticalLayout") == "true");
  if (useVerticalLayout === null) {
    useVerticalLayout = true;
    localStorage.setItem("verticalLayout", useVerticalLayout);
  }
  if (useVerticalLayout) {
    tree.layout(tnt.tree.layout.vertical().width(width));
  } else {
    tree.layout(tnt.tree.layout.radial().width(width));
  }

  // set node colors and labels
  tree.node_display(tnt.tree.node_display.circle()
    .size(6)
    .fill(function(node) {
      var gen = node.node_name().substring(0, 3)
      if (gen == 'USR') return(genusToColor[gen])
      var gensp = node.node_name().substring(0, 5)
      if (gensp in genspToTaxon) return(taxonToColor(genspToTaxon[gensp]))
      else if (node.is_leaf()) return(DEFAULT_COLOR);
      else if (node.root_dist() == 0) return('black');
      return 'white';
    })
  );
  tree.label(tnt.tree.label.text()
    .height(16)
    .text(function(node) {
      if (!node.is_leaf()) return('');
      return node.node_name();
    })
  );

  var treeDiv = document.getElementById(elementId);
  treeDiv.innerHTML = ""; // to clear it
  tree(treeDiv);

  // use TnT Tree API to calculate intergenic distance over some screen distance.
  $("#phylotreeDistanceScale").empty(); // to clear it
  d3.select("#phylotreeDistanceScale")
    .attr('width', width)
    .attr('height', 40)
    .append('g')
    .attr('transform', 'translate(20, 20)')
    .attr('class', 'x axis')
    .call(getXAxis(tree, width)
  );
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

// =============================================================

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

shinyjs.displayTaxaView = function(args) {
  var taxa = args[0];
  var elementId = args[1];
  var elementTag = "#" + elementId;
  var tree = args[2];

  // transform to objects expected by nvd3 chart
  var x = countTaxa(tree);
  var data = _.map(Object.keys(x).sort(), k => {
    return { label: k, value: x[k], color: taxonToColor(k) };
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

    return chart;
  });

  doResize(); // does not always work here?
}

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

  // filter MSA viewer and phylotree by selected taxa
  // (TODO: should really notify all registered event listeners)
  filterMSAView(taxa, "msa");
  filterPhylotree(taxa, "phylotree");
}

shinyjs.resetTaxaView = function(args) {
  elementId = args[0];
  elementTag = "#" + elementId;
  var legend = d3.select(elementTag);
  var data = legend.datum();
  data.forEach(function(series) { series.disabled = false });

  handleTaxaSelection(elementTag);
  doResize();
}

// force a resize, to redraw the chart
function doResize() {
  window.dispatchEvent(new Event("resize"));
}

// =============================================================

function filterMSAView(taxa, elementId) {
  document.getElementById(elementId).innerHTML = "";
  var seqs = JSON.parse(localStorage.getItem("msaSeqs"));
  seqs = seqs.filter(function(seq) {
    var gensp = seq.name.substring(0, 5);
    if (gensp.startsWith("USR")) gensp = "USR";
    return taxa.includes(gensp);
  });
  var m = msa({
    el: document.getElementById(elementId),
    //bootstrapMenu: true, // TODO: get buttons working
    seqs: seqs,
    vis: {
      overviewbox: false,
      labelId: false
    },
    zoomer: {
      labelNameLength: 150,
      labelFontsize: 9,
      autoResize: true
    }
  });
  m.render();
}

shinyjs.displayMSAView = function(args) {
  var seqs = msa.io.fasta.parse(args[0]);
  localStorage.setItem("msaSeqs", JSON.stringify(seqs));
  var elementId = args[1];

  var m = msa({
    el: document.getElementById(elementId),
    //bootstrapMenu: true, // TODO: get buttons working
    seqs: seqs,
    vis: {
      overviewbox: false,
      labelId: false
    },
    zoomer: {
      labelNameLength: 150,
      labelFontsize: 9,
      autoResize: true
    }
  });
  m.render();
}

// =============================================================

